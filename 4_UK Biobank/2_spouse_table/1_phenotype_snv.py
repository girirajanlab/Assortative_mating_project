import pandas as pd
import numpy as np

# Filter spouses to remove same-sex pairs, related pairs, and pairs with the same parental age of death
df=pd.read_csv('../../0_preprocessing/2_kinship_PCs/UK_Biobank/Kinship/list_files/1_spouses.csv')
df=df[(df.detail.isnull()) | (df.detail=='Missing genotype data in pair after QC')]
df.IID=df.IID.astype(int)

# Add in phenotype data
# 20544 Mental health problems ever diagnosed by a professional

# Add phenotype information
# Get columns of relevant phenotypes
file = open('download/ukb30075.csv', 'r', encoding='unicode_escape')
for line in file:
        header = line.rstrip().split(',')
        break
file.close()

cols = [i.replace('"', '') for i in header if '20544-' in i]

# 31 is sex, 22006 is European ancestry
pheno_df = pd.read_csv('download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid', '31-0.0', '22006-0.0']+cols)
pheno_df.eid=pheno_df.eid.astype(int)
pheno_df = pheno_df[pheno_df.eid.isin(df.IID.to_list())]
pheno_df.dropna(how='all', axis=1, inplace=True)

# Add in sex
pheno_df['Sex']=pheno_df['31-0.0'].map({0:'female', 1:'male'})
df['sex']=df.IID.map(dict(zip(pheno_df.eid.to_list(), pheno_df.Sex.to_list())))

# Add in European ancestry
df['European']=df.IID.map(dict(zip(pheno_df.eid.to_list(), pheno_df['22006-0.0'].to_list())))

# 20544 uses data-coding 1401
data_coding = pd.read_csv('Analysis_files/coding1401.tsv', sep = '\t', index_col = 0)
pheno_df['diagnoses'] = '.'
for col in cols:
        pheno_df[col] = pheno_df[col].map(data_coding.meaning.to_dict())
        pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses=='.'), 'diagnoses'] = ''
        pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses!=''), 'diagnoses'] = pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses!=''), 'diagnoses'] + ';' + pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses!=''), col]
        pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses==''), 'diagnoses'] = pheno_df.loc[(~pheno_df[col].isnull()) & (pheno_df.diagnoses==''), col]

meanings=data_coding.meaning.to_list()
for m in meanings:
	pheno_df[m]=0
	# The () in some meanings causes issues for pandas - annotate these separately
	if 'group A' in m:
		pheno_df.loc[pheno_df.diagnoses.str.contains('group A'), 'Prefer not to answer (group A)']=1
	elif 'group B' in m:
		pheno_df.loc[pheno_df.diagnoses.str.contains('group B'), 'Prefer not to answer (group B)']=1
	elif 'OCD' in m:
		pheno_df.loc[pheno_df.diagnoses.str.contains('OCD'), 'Obsessive compulsive disorder (OCD)']=1
	elif 'Any other phobia' in m:
		pheno_df.loc[pheno_df.diagnoses.str.contains('Any other phobia'), 'Any other phobia (eg disabling fear of heights or spiders)']=1
	elif 'ADHD' in m:
		pheno_df.loc[pheno_df.diagnoses.str.contains('ADHD'), 'Attention deficit or attention deficit and hyperactivity disorder (ADD/ADHD)']=1
	else:
		pheno_df.loc[pheno_df.diagnoses.str.contains(m), m]=1

# Add NAs based on "prefer not to answer"
groupa=['Depression', 'Mania, hypomania, bipolar or manic-depression', 'Anxiety, nerves or generalized anxiety disorder', 'Social anxiety or social phobia',
	'Agoraphobia', 'Any other phobia (eg disabling fear of heights or spiders)', 'Panic attacks', 'Obsessive compulsive disorder (OCD)']
groupb=['Anorexia nervosa', 'Bulimia nervosa', 'Psychological over-eating or binge-eating', 'Schizophrenia', 'Any other type of psychosis or psychotic illness', 'A personality disorder',
	"Autism, Asperger's or autistic spectrum disorder", 'Attention deficit or attention deficit and hyperactivity disorder (ADD/ADHD)']
for g in groupa:
	pheno_df.loc[pheno_df['Prefer not to answer (group A)']==1, g]=np.nan
for g in groupb:
	pheno_df.loc[pheno_df['Prefer not to answer (group B)']==1, g]=np.nan

# Add phenotypes of interest to spouse file
pheno_df['IID']=pheno_df.eid
df=pd.merge(df, pheno_df[['IID', 'Anxiety, nerves or generalized anxiety disorder', 'Mania, hypomania, bipolar or manic-depression', 'Depression', 'A personality disorder', 'Schizophrenia']], on='IID', how='left')

# Add SNV burden
snv=pd.read_csv('../../0_preprocessing/1_variant_calling/UK_Biobank/synonymous/tables/17_burden_table.csv')
snv['IID'] = snv.Sample.astype(int)
df=pd.merge(df, snv[['IID', 'burden_nonsynonymous', 'burden_synonymous']], on='IID', how='left')

# Save to file
print(df)
df.columns=['psuedoid', 'Couple', 'detail', 'IID', 'sex', 'European', 'Anxiety', 'BPD_mania', 'Depression', 'Personality disorder', 'Schizophrenia', 'nonsynonymous_burden', 'synonymous_burden']
df=df[['IID', 'Couple', 'detail', 'sex', 'European', 'Anxiety', 'BPD_mania', 'Depression', 'Personality disorder', 'Schizophrenia', 'nonsynonymous_burden', 'synonymous_burden']]
df.to_csv('Analysis_files/1_phenotype_snv.csv', index=False)
