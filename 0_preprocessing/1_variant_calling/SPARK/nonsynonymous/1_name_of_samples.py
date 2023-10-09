#!/bin/python3
import pandas as pd
import sys
import os

spark_mastertable=sys.argv[1]
df = pd.read_csv(spark_mastertable, sep='\t')

df['Sample'] = df['spid']

df['VCF_filename'] = df['wes_gvcf']
df['Last_digit'] = df['Sample'].apply(lambda s: s[-1])

df = df[['Sample', 'Last_digit', 'VCF_filename', 'father', 'mother']]

df.to_csv('samples.csv', index=False)
