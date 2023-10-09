import pandas as pd
import numpy as np

# Reorganize data by spouses
df=pd.read_csv('Analysis_files/1_phenotype_snv.csv')
df['phenotype_data']='X'
df.loc[(df.Anxiety.isnull()) & (df.BPD_mania.isnull()) & (df.Depression.isnull()) & (df['Personality disorder'].isnull()) & (df.Schizophrenia.isnull()), 'phenotype_data']=np.nan

m_df=df[df.sex=='male']
f_df=df[df.sex=='female']

spouse_df=pd.merge(m_df, f_df, on=['Couple', 'detail'], how='outer', suffixes=['_male', '_female'])
print(spouse_df)
print(spouse_df.columns)

spouse_df=spouse_df[['Couple', 'detail', 'IID_male', 'IID_female', 'European_male', 'European_female',
			'Anxiety_male', 'BPD_mania_male', 'Depression_male', 'Personality disorder_male', 'Schizophrenia_male',
			'nonsynonymous_burden_male', 'synonymous_burden_male', 'phenotype_data_male',
			'Anxiety_female', 'BPD_mania_female', 'Depression_female', 'Personality disorder_female', 'Schizophrenia_female',
			'nonsynonymous_burden_female', 'synonymous_burden_female', 'phenotype_data_female']]
spouse_df.columns=['Couple', 'detail', 'Male', 'Female', 'European_male', 'European_female',
			'Anxiety_male', 'BPD_mania_male', 'Depression_male', 'Personality disorder_male', 'Schizophrenia_male',
			'nonsynonymous_burden_male', 'synonymous_burden_male', 'phenotype_data_male',
			'Anxiety_female', 'BPD_mania_female', 'Depression_female', 'Personality disorder_female', 'Schizophrenia_female',
			'nonsynonymous_burden_female', 'synonymous_burden_female', 'phenotype_data_female']

# Add kinship
kin=pd.read_csv('../../0_preprocessing/2_kinship_PCs/UK_Biobank/Kinship/king.kin', sep='\t')
kin.index=kin.FID

spouse_df['kinship_coeff']=spouse_df.Couple.map(kin.Kinship.to_dict())

# Annotate pairs with use for anlyses
spouse_df['phenotype_analysis']=np.nan
spouse_df.loc[(spouse_df.phenotype_data_male=='X') & (spouse_df.phenotype_data_female=='X'), 'phenotype_analysis']='X'

spouse_df['burden_analysis']=np.nan
spouse_df.loc[(~spouse_df.nonsynonymous_burden_male.isnull()) & (~spouse_df.nonsynonymous_burden_female.isnull()), 'burden_analysis']='X'

spouse_df['kinship_burden_analysis']=np.nan
spouse_df.loc[(spouse_df.burden_analysis=='X') & (~spouse_df.kinship_coeff.isnull()), 'kinship_burden_analysis']='X'

spouse_df['european_pair']=np.nan
spouse_df.loc[(spouse_df.European_male==1) & (spouse_df.European_female==1), 'european_pair']='X'

for c in ['phenotype_analysis', 'burden_analysis', 'kinship_burden_analysis']:
	print(spouse_df[c].value_counts())

# Save to file
spouse_df.to_csv('UKB_family_table.csv')
