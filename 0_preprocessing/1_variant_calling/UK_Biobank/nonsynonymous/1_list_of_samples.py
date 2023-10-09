#!/bin/python3
import pandas as pd
from pandarallel import pandarallel

pandarallel.initialize(10)

root_dir = '.'
filename = 'ukb48799_23151.bulk'
with open(filename, 'r') as f:
	bulk_file_lines = f.readlines()

bulk_file_lines = [s.strip() for s in bulk_file_lines]

df = pd.DataFrame(index=bulk_file_lines)
df['Bulk_file_line'] = bulk_file_lines
df['Sample'] = df['Bulk_file_line'].apply(lambda s: s.split(' ')[0])
df['Last_two_digits'] = df['Sample'].apply(lambda s: s[-2:])
df['VCF_filename'] = '/data5/UK_Biobank/release_2021_oct/48799/' + df['Last_two_digits'] + '/' + df['Sample'] + '_23151_0_0.gz'

filename = f'{root_dir}/samples.csv'
df.to_csv(filename, index=False)
