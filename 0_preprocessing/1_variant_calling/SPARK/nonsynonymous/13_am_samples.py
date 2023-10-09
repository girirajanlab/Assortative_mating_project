import pandas as pd
import numpy as np

# Get just samples needed for assortative mating analysis
# Use FAM files to get pairs and children
fam1=pd.read_csv('analysis_files/SPARK_DATA1.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
fam2=pd.read_csv('analysis_files/SPARK_DATA2.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])

fam=pd.concat([fam1, fam2])

# Get sex information
sex_dict=dict(zip(fam.IID.to_list(), fam.Sex.to_list()))

# Keep only complete parent pairs
print(fam.shape)
fam=fam[(fam.Mother!='0') & (fam.Father!='0')]
print(fam.shape)

# Create a table with all spouse pairs from SPARK
spouse_tab=fam[['FID', 'Mother', 'Father']].copy()
# Remove duplicate rows
spouse_tab.drop_duplicates(inplace=True)

spouse_tab['spousepair']=spouse_tab.Mother+'.'+spouse_tab.Father

# Get SNV burden
snv=pd.read_csv('analysis_files/12_sample_burden.csv')
burden_dict=dict(zip(snv.Sample.to_list(), snv.burden.to_list()))
spouse_tab['female_burden']=spouse_tab.Mother.map(burden_dict)
spouse_tab['male_burden']=spouse_tab.Father.map(burden_dict)

spouse_tab['mean_burden']=(spouse_tab.female_burden+spouse_tab.male_burden)/2

spouse_tab['spouse_SNV']='X'
spouse_tab.loc[(spouse_tab.female_burden.isnull()) | (spouse_tab.male_burden.isnull()), 'spouse_SNV']=np.nan
print(spouse_tab.spouse_SNV.value_counts())

# Make table with burden and kinship values
kin=pd.read_csv('analysis_files/SPARK_spouse.txt', sep='\t')

kin['order1']=kin.ID1+'.'+kin.ID2
kin['order2']=kin.ID2+'.'+kin.ID1

print(kin[~((kin.order1.isin(spouse_tab.spousepair.to_list())) | (kin.order2.isin(spouse_tab.spousepair.to_list())))])

kinship_dict=dict(zip(kin.order1.to_list()+kin.order2.to_list(), kin.Kinship.to_list()+kin.Kinship.to_list()))
spouse_tab['kinship_coeff']=spouse_tab.spousepair.map(kinship_dict)

spouse_tab['kinship']='X'
spouse_tab.loc[spouse_tab.kinship_coeff.isnull(), 'kinship']=np.nan
print(spouse_tab.kinship.value_counts())

# Add children
fam['spousepair']=fam.Mother+'.'+fam.Father
child_dict={}
for s in list(fam.spousepair.unique()):
	child_dict[s]=';'.join(list(fam[fam.spousepair==s]['IID'].unique()))
spouse_tab['child_ID']=spouse_tab.spousepair.map(child_dict)

# Add proband
roles=pd.read_csv('analysis_files/roles_index.csv')
print(roles)
probands=roles[roles.proband]['subject_sp_id'].to_list()
proband_only=fam[fam.IID.isin(probands)]
print(proband_only.spousepair.value_counts())
proband_dict={}
for s in list(proband_only.spousepair.unique()):
	proband_dict[s]=';'.join(list(proband_only[proband_only.spousepair==s]['IID'].unique()))
spouse_tab['proband']=spouse_tab.spousepair.map(proband_dict)

# Add proband burden
spouse_tab['proband_burden']=spouse_tab.proband.map(burden_dict)

# Save to file
spouse_tab.to_csv('SPARK_AM_spouse_table.csv', index=False)
