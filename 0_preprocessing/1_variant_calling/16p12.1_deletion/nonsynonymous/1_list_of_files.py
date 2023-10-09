#!/bin/python3
import pandas as pd
import os

filenames = os.listdir('final_calls/annovar/')
filenames = [s for s in filenames1 if s.endswith('_input.vcf')]
filenames = [f'final_calls/annovar/{s}' for s in filenames1]

df = pd.DataFrame(index=filenames)
df['Filename'] = df.index.to_series()
df['Sample'] = df['Filename'].apply(lambda s: s.split('/')[-1].split('_')[0])

df.to_csv('samples.csv', index=False)
