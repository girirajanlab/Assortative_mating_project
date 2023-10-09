#!/bin/python
import pandas as pd

# Use the SPARK samples file from Dropbox to get the list of proband samples we need
df = pd.read_csv('Analysis_files/1_individual_table.csv')
print(df)
df=df[df.ROH_proband]
print(df)

# Reformat as FAM file
df['Sex']=df.sex.map({'Male':1, 'Female':2})
df=df[['FID', 'IID', 'Father', 'Mother', 'Sex']]
df['Phenotype']=0

df.to_csv('Analysis_files/1_roh_samples.fam', sep = ' ', index = False, header = False)

