#!/bin/python
import pandas as pd

# Merge cohorts into single file
# Remove calls that are:
# 1. < 50kb, or
# 2. Contain < 5 SNPs

mv1=pd.read_csv('tables/1Mv1_merged.txt', sep='\t')
mv3=pd.read_csv('tables/1Mv3_merged.txt', sep='\t')
omni=pd.read_csv('tables/Omni2.5_merged.txt', sep='\t')

out_df=pd.concat([mv1, mv3, omni])
out_df = out_df[(out_df.Length>=50000) & (out_df.NumSNP>=5)]

# Save to file
out_df.to_csv('tables/7_size_filtered.txt', sep='\t', index=False)
