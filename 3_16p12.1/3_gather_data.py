import pandas as pd
import numpy as np

# Compile genetic data for spouse pairs
df=pd.read_csv('Analysis_files/1_spouse_phenotypes.csv')

cohort=pd.read_excel('Analysis_files/Table_S1_Cohort_description.xlsx')
cohort.index=cohort.Sample.to_list()

df['Mother']=df.spousepair.str.split('.', expand=True)[0]
df['Father']=df.spousepair.str.split('.', expand=True)[1]
df['FID']=df.Mother.map(cohort.Family.to_dict())
df=df[['FID', 'Mother', 'Father', 'spousepair', 'spouse_phenotypes']]

# Save to file
print(df)
df.to_csv('Analysis_files/3_genetic_data.csv', index=False)

