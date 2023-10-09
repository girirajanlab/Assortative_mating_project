import pandas as pd
import numpy as np

# Get all SSC proband and parent samples
samp_df=pd.read_csv('SSC data/Phenotype Data Set 9/Individuals_by_Distribution_v15.csv')

df=samp_df[samp_df['SSC ID'].str.contains('.p1')].copy()
df['FID']=df['SSC ID'].str.split('.', expand=True)[0]
df['proband']=df['SSC ID']
df['Mother']=df.FID+'.mo'
df['Father']=df.FID+'.fa'

# Confirm that both parents are in collection
df.loc[~df.Mother.isin(samp_df['SSC ID'].to_list()), 'Mother']=np.nan
df.loc[~df.Father.isin(samp_df['SSC ID'].to_list()), 'Father']=np.nan
df=df[(~df.Mother.isnull()) & (~df.Father.isnull())]

# Add additional family information
df['spousepair']=df.Mother+'.'+df.Father

df=df[['FID', 'Mother', 'Father', 'spousepair', 'proband']]
print(df.shape)

# Gather quantitative phenotypes from parents and children
# Phenotypes and locations are listed in Analysis_files/Warrier_ASD_phenotypes.xlsx

# Proband data
pro_df=pd.read_excel('Analysis_files/Linear_model_phenotypes.xlsx', sheet_name='Proband_features')
# Save names for later
pro_phenos=[]
for idx, row in pro_df.iterrows():
    filename=row['SSC_filename']
    if 'Analysis_files' not in filename:
        filename=dropbox+filename
    mdf=pd.read_csv(filename)
    mdf[row['Measure']]=mdf[row['SSC_colname']].apply(pd.to_numeric, errors='coerce')
    mdf['proband']=mdf.individual
    df=pd.merge(df, mdf[['proband', row['Measure']]], on='proband', how='outer')
    pro_phenos.append(row['Measure'])

# Parent data
par_df=pd.read_excel('Analysis_files/Linear_model_phenotypes.xlsx', sheet_name='Parent_features')
for idx, row in par_df.iterrows():
    for p in ['Mother', 'Father']:
        filename='SSC data/Phenotype Data Set 9/'+p+' Data/'+row['SSC_filename']
        mdf=pd.read_csv(filename)
        mdf[p+'_'+row['Measure']]=mdf[row['SSC_colname']]
        mdf['proband']=mdf.individual.str.split('.', expand=True)[0]+'.p1'
        df=pd.merge(df, mdf[['proband', p+'_'+row['Measure']]], on='proband', how='outer')

# Calculate mean parental SRS and BAPQ
df['Biparental mean SRS']=(df.Mother_SRS + df.Father_SRS)/2
df['Biparental mean BAPQ']=(df.Mother_BAPQ+df.Father_BAPQ)/2

# Add rare variant burden
counts=pd.read_csv('Analysis_files/SSC_burden_table.csv')
counts.index=counts.Sample.to_list()
burden_dict=counts.burden_nonsynonymous.to_dict()
df['SNV burden']=df.proband.map(burden_dict)

# Annotate for presence of mutation in Tier S SFARI gene in probands
vars = pd.read_csv('analysis_files/SSC_variants.csv')
gene_annos = pd.read_csv("Gene_annotations.csv")
tier_s = gene_annos[gene_annos.SFARI_gene_score=='S']['gene_id']
missense_sfari_samps=list(vars[(vars.Mut_type=='missense') & (vars.Gene_id_.isin(tier_s))]['Sample'].unique())
lof_sfari_samps=list(vars[(vars.Mut_type=='lof') & (vars.Gene_id_.isin(tier_s))]['Sample'].unique())
splice_sfari_samps=list(vars[(vars.Mut_type=='splice') & (vars.Gene_id_.isin(tier_s))]['Sample'].unique())
df['Tier S SFARI SNV']=np.nan
df.loc[(~df['SNV burden'].isnull()), 'Tier S SFARI SNV']=0
df.loc[df.proband.isin(splice_sfari_samps), 'Tier S SFARI SNV']='Splice'
df.loc[df.proband.isin(missense_sfari_samps), 'Tier S SFARI SNV']='Missense'
df.loc[df.proband.isin(lof_sfari_samps), 'Tier S SFARI SNV']='LOF'

# Add total length del and dup from large CNVs in probands
cnvs=pd.read_csv('SSC_Sanders_CNVs.csv')
cnvs_done=pd.read_excel('Analysis_files/Sanders_TableS1.xlsx')
cnvs_samps=cnvs_done[cnvs_done['CNV_ThisManuscriptSsc']!='0']['Proband'].to_list()
del_dict={}
dup_dict={}
for p in df.proband.to_list():
    if p not in cnvs_samps:
        continue
    del_len=sum(cnvs[(cnvs.patientID==p) & (cnvs['Del/Dup']=='Del')]['Size'].to_list())
    del_dict[p]=del_len
    dup_len=sum(cnvs[(cnvs.patientID==p) & (cnvs['Del/Dup']=='Dup')]['Size'].to_list())
    dup_dict[p]=dup_len

df['Del (bp)']=df.proband.map(del_dict)
df['Dup (bp)']=df.proband.map(dup_dict)

# Add autism PRS
prs = pd.read_csv('Analysis_files/SummaryGeneticData_REACH_SSC_SPARK.SBayesR.20210915.csv')
prs=prs[(prs.Cohort=='SSC') & (~prs.IID.str.contains('REACH'))]
sample_map=pd.read_csv('Analysis_files/nygc_sfari_id_sample_map.csv')
prs.index=prs.IID.map(dict(zip(sample_map[' Sample ID'].to_list(), sample_map['SFARI ID'].to_list())))
df['Autism PRS']=df.proband.map(prs['PS ASD'].to_dict())
df['Mother_PRS']=df.proband.map(prs['PS ASD Mother'].to_dict())
df['Father_PRS']=df.proband.map(prs['PS ASD Father'].to_dict())

# Also add in proband sex and age
desc=pd.read_csv("SSC data/Phenotype Data Set 9/Proband Data/ssc_core_descriptive.csv")
desc.index=desc.individual.to_list()
df['Sex']=df.proband.map(desc.sex.to_dict())
df['Age']=df.proband.map(desc.age_at_ados.to_dict())

# Add in ancestry PCs
pcs=pd.read_csv('Analysis_files/EIGENSTRAT_SSC_ancestry_pcs.csv')
df=pd.merge(df, pcs, left_on='proband', right_on='IID', how='outer')

# Add in parental kinship
kin=pd.read_csv('Analysis_files/SSC_king.kin', delim_whitespace=True)
kin['spousepair1']=kin.ID1+'.'+kin.ID2
kin['spousepair2']=kin.ID2+'.'+kin.ID1
kinship_dict=dict(zip(kin.spousepair1.to_list()+kin.spousepair2.to_list(), kin.Kinship.to_list()+kin.Kinship.to_list()))
df['kinship_coeff']=df.spousepair.map(kinship_dict)

# Annotate families with the anlyses they can be used in
# Linear models
# Need all inputs for probands (all 10 genetic pcs, proband sex,  proband age, autism PRS, SFARI mutation, SNV burden, bp del, bp dup, parental mean srs, parental mean bapq) and at least one response factor
df['linear_models']='.'
for f in pro_phenos:
    df.loc[(~df[f].isnull()), 'linear_models']='X'
df.loc[(df.PC1.isnull()) | (df.Sex.isnull()) | (df.Age.isnull()) | (df['Autism PRS'].isnull()) | (df['Tier S SFARI SNV'].isnull()) | (df['SNV burden'].isnull()) | (df['Del (bp)'].isnull()) | (df['Dup (bp)'].isnull()) | (df['Biparental mean SRS'].isnull()) | (df['Biparental mean BAPQ'].isnull()), 'linear_models']=np.nan
df.replace('.', np.nan, inplace=True)

# Parent phenotype comparisons
eig=pd.read_csv('Analysis_files/SSC_snv_table.csv')
df['Parent_phenotype_analysis']=np.nan
df.loc[((~df['Biparental mean SRS'].isnull()) | (~df['Biparental mean BAPQ'].isnull())) & (df.Mother.isin(eig.Sample.to_list())) & (df.Father.isin(eig.Sample.to_list())), 'Parent_phenotype_analysis']='X'

# Kinship and deleterious SNV comparisons
df['kinship_SFARI_analysis']=np.nan
df.loc[(~df.kinship_coeff.isnull()) & (~df['Tier S SFARI SNV'].isnull()) & (df.kinship_coeff>-0.1), 'kinship_SFARI_analysis']='X'

# Also annotate probands that can be used for ROH analysis
df['ROH_proband']=False
df.loc[(~df.kinship_coeff.isnull()), 'ROH_proband']=True

# Save to file
df.to_csv('Analysis_files/2_SSC_data.csv', index=False)