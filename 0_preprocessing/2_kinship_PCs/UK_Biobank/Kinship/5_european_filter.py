import pandas as pd

# Filter spouses for only those of European ancestry
# Do this using Data Field 22006 from the ukb30075.csv file
# This field is "Genetic ethnic group" and
# "Indicates samples who self-identified as 'White British' according to Field 21000 and have very similar genetic ancestry based on a principal components analysis of the genotypes."

samps=pd.read_csv('plink_files/4_filter/UKB_QC_het.valid.sample', sep='\t')

pheno_df = pd.read_csv('download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid', '22006-0.0'])
pheno_df.columns=['IID', 'European']

df=pd.merge(samps, pheno_df, on='IID', how='left')
df=df[df.European==1]

# Check for complete spouse pairs
spouses=pd.read_csv('list_files/1_spouses.csv')
df=pd.merge(df, spouses, on='IID', how='left')

spousepairs=df.Couple.value_counts()
comp_pairs=spousepairs[spousepairs==2].index.to_list()
df=df[df.Couple.isin(comp_pairs)]

# SAve needed samples as FAM
df[['FID', 'IID']].to_csv('list_files/5_european_pairs.fam', sep=' ', index=False)

# Also make a new file for updating IDs
df[['FID', 'IID', 'Couple', 'IID']].to_csv('list_files/5_update_fid.list', sep=' ', index=False)

