#!/bin/python3
import pandas as pd
import subprocess
import multiprocessing as mp

df = pd.read_csv('samples.csv')
df.sort_values('Last_two_digits', inplace=True)

startfile_directory='intermediate_files/'
output_dir='vcfs/by_sample/'

for i in range(100):
	command = f'mkdir -p {output_dir}/{i:02d}'
	# print(command)
	subprocess.run(command, shell=True)

def process(i):
	sample_name = df.at[i, 'Sample']
	invcf = f'{startfile_directory}/{str(sample_name)[-2:]}/{sample_name}/{sample_name}.1_filter1.vcf'
	outvcf = f'{output_dir}/{str(sample_name)[-2:]}/{sample_name}.vcf.gz'

	command = f'cat {invcf} | bgzip > {outvcf} ; tabix -p vcf {outvcf}'
	subprocess.run(command, shell=True)

print(df)
num_samples = df.shape[0]

pool = mp.Pool(100)
pool.map(process, [s for s in range(num_samples)])

