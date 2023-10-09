#!/bin/python
import pandas as pd
import numpy as np

# Gather families from SPARK cohort
fam=pd.read_csv('SPARK data/SPARK_Collection_version5/roles_index.csv', low_memory=False)
fam=fam[['subject_sp_id', 'family_sf_id', 'sex', 'proband', 'biomother_sp_id', 'biofather_sp_id']]
fam.columns=['IID', 'FID', 'sex', 'is_proband', 'Mother', 'Father']

# Restrict to probands and their parents
df=fam[(fam.is_proband) & ((~fam.Mother.isnull()) | (~fam.Father.isnull()))]

df=df[['FID', 'Mother', 'Father', 'IID', 'sex']]
df.columns=['FID', 'Mother', 'Father', 'proband', 'proband_sex']
df['spousepair']=np.nan
df.loc[(~df.Mother.isnull()) & (~df.Father.isnull()), 'spousepair']=df.Mother+'.'+df.Father

# Also make a table listing all individuals separately
indiv_tab=pd.DataFrame({'FID':df.FID.to_list()*3, 'IID':df.Mother.to_list()+df.Father.to_list()+df.proband.to_list(), 'sex':['Female']*df.shape[0]+['Male']*df.shape[0]+df.proband_sex.to_list(),
                        'role':['mother']*df.shape[0]+['father']*df.shape[0]+['proband']*df.shape[0]})
indiv_tab=indiv_tab[~indiv_tab.IID.isnull()]
indiv_tab['spousepair']=indiv_tab.IID.map(dict(zip(df.Mother.to_list()+df.Father.to_list()+df.proband.to_list(), df.spousepair.to_list()*3)))
indiv_tab.sort_values(by='FID', inplace=True)

# Gather phenotype data

# We want the following phenotypes:
# behav_adhd	  ADHD (Attention Deficit-Hyperactivity Disorder) or ADD
# mood_anx	 Anxiety disorder, such as panic, phobia, agoraphobia, or generalized anxiety disorder (GAD) except for social anxiety 
# mood_bipol	 Bipolar (Manic-Depressive) Disorder
# mood_dep	 Depression or dysthymia
# mood_ocd	 Obsessive-Compulsive Disorder
# pers_dis	 Personality Disorder
# sleep_dx	Sleep Disorder or sleep problem diagnosed by a professional

phenos = ['age_at_eval_months', 'behav_adhd', 'mood_anx', 'mood_bipol', 'mood_dep', 'mood_ocd', 'pers_dis', 'sleep_dx']

med_screen = pd.read_csv("SPARK data/SPARK_Collection_version5/basic_medical_screening.csv", low_memory=False)
med_screen.fillna(0, inplace=True)
med_screen = med_screen[med_screen.subject_sp_id.isin(indiv_tab.IID.to_list())]

# Add phenotype information
for i in phenos:
    indiv_tab[i]=indiv_tab.IID.map(dict(zip(med_screen.subject_sp_id.to_list(), med_screen[i].to_list())))

indiv_tab['phenotype']=np.nan
for i in phenos[1:]:
    indiv_tab.loc[~(indiv_tab[i].isnull()), 'phenotype'] = 'X'

# Remove parents from individual table if they do not have any phenotype information and are not part of a complete parent pair
indiv_tab=indiv_tab[(indiv_tab.role=='proband') | (~indiv_tab.phenotype.isnull()) | (~indiv_tab.spousepair.isnull())]
# Remove probands if they do not have a parent in the table
indiv_tab['Mother']=indiv_tab.IID.map(dict(zip(df.proband.to_list(), df.Mother.to_list())))
indiv_tab['Father']=indiv_tab.IID.map(dict(zip(df.proband.to_list(), df.Father.to_list())))
indiv_tab=indiv_tab[(indiv_tab.role!='proband') | (indiv_tab.Mother.isin(indiv_tab.IID.to_list())) | (indiv_tab.Father.isin(indiv_tab.IID.to_list()))]

# Save phenotype data as separate table condensed by family
pro_df = indiv_tab[indiv_tab.role=='proband'].copy()
m_df=indiv_tab[indiv_tab.role=='mother'].copy()
f_df=indiv_tab[indiv_tab.role=='father'].copy()

pm = pd.merge(pro_df[['FID', 'sex']+phenos], m_df[['FID']+phenos], on='FID', how='outer', suffixes=['_proband', '_mother'])
pf=pd.merge(pro_df[['FID']+phenos], f_df[['FID']+phenos], on='FID', how='outer', suffixes=['_proband', '_father'])

pmf=pd.merge(pm, pf, on=['FID']+[i+'_proband' for i in phenos], how='outer')

pmf['proband_phenotype']=np.nan
pmf['mother_phenotype']=np.nan
pmf['father_phenotype']=np.nan
for p in phenos:
    for g in ['proband', 'mother', 'father']:
        pmf.loc[~pmf[p+'_'+g].isnull(), g+'_phenotype']='X'

# Ensure all families in file can be used for parent-child comparisons or spouse-spouse comparisons
pmf['parent_child']=np.nan
pmf.loc[(pmf.proband_phenotype=='X') & ((pmf.mother_phenotype=='X') | (pmf.father_phenotype=='X')), 'parent_child']='X'

pmf['spouse']=np.nan
pmf.loc[((pmf.mother_phenotype=='X') & (pmf.father_phenotype=='X')), 'spouse']='X'

pmf=pmf[(pmf.parent_child=='X') | (pmf.spouse=='X')]
pmf.to_csv('Analysis_files/1_family_phenotypes.csv', index=False)

# Restrict family table to families remaining in individual table
df=df[(df.Mother.isin(indiv_tab.IID.to_list())) | (df.Father.isin(indiv_tab.IID.to_list()))]

# Annotate families with whether they will be used in phenotypic analysis (both parents have phenotype data or parent and child have phenotype data)
pheno_dict=dict(zip(indiv_tab.IID.to_list(), indiv_tab.phenotype.to_list()))
df['Mother_phenotype']=df.Mother.map(pheno_dict)
df['Father_phenotype']=df.Father.map(pheno_dict)
df['proband_phenotype']=df.proband.map(pheno_dict)

df['Parent_child_pheno']=np.nan
df.loc[(~df.proband_phenotype.isnull()) & ((~df.Mother_phenotype.isnull()) | (~df.Father_phenotype.isnull())), 'Parent_child_pheno']='X'

df['Spouse_pheno']=np.nan
df.loc[(~df.Mother_phenotype.isnull()) & (~df.Father_phenotype.isnull()), 'Spouse_pheno']='X'

# Add kinship information to family table
kin=pd.read_csv('Analysis_files/SPARK_king.kin', delim_whitespace=True)
kin=kin[['FID', 'ID1', 'ID2', 'Kinship']]
# Remove kinship coefficients <= -0.1
kin=kin[kin.Kinship>-0.1]
kin['spousepair1']=kin.ID1+'.'+kin.ID2
kin['spousepair2']=kin.ID2+'.'+kin.ID1
kin_dict=dict(zip(kin.spousepair1.to_list()+kin.spousepair2.to_list(), kin.Kinship.to_list()+kin.Kinship.to_list()))
df['kinship_coeff']=df.spousepair.map(kin_dict)

# Remove probands whose cannot be used in phenotype analysis (parent and child both have phenotype data) or ROH-kinship analysis (parents have kinship coeff)
kin_spouses=df[~df.kinship_coeff.isnull()]['spousepair'].to_list()
pheno_fams=df[df.Parent_child_pheno=='X']['FID'].to_list()
indiv_tab=indiv_tab[(indiv_tab.role!='proband') | (indiv_tab.FID.isin(pheno_fams)) | (indiv_tab.spousepair.isin(kin_spouses))]

# Annotate probands to use in ROH analysis
indiv_tab['ROH_proband']=False
indiv_tab.loc[(indiv_tab.role=='proband') & (indiv_tab.spousepair.isin(kin_spouses)), 'ROH_proband']=True

# Save tables to file
indiv_tab.to_csv('Analysis_files/1_individual_table.csv', index = False)
df.to_csv('Analysis_files/1_family_table.csv', index=False)