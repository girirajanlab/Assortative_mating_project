#!/bin/python
import pandas as pd
from os.path import exists

# Merge all calls into single file
# Remove calls that are:
# 1. < 50kb, or
# 2. Contain < 5 SNPs

samps=pd.read_csv('Analysis_files/1_samples.csv')

df=pd.DataFrame()
for s in samps.IID.to_list():
	if exists('tables/6_merge_individual/'+s+'.txt'):
		sdf=pd.read_csv('tables/6_merge_individual/'+s+'.txt', sep='\t')
		df=pd.concat([df, sdf])

df = df[(df.Length>=50000) & (df.NumSNP>=5)]

# Save to file
df.to_csv('tables/7_size_filtered.txt', sep='\t', index=False)
