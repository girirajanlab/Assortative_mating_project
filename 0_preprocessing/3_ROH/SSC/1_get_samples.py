#!/bin/python
import pandas as pd

# As this analysis is a companion to the assortative mating analysis,
# we only want to calculate ROH for children of spouse pairs used in that analysis

# Use the SSC spousepairs file from Dropbox to get the list of proband samples we need
spouse = pd.read_csv('Analysis_files/ssc_spouse.txt', sep = '\t')

df = pd.DataFrame(spouse.FID1.str.split('.', expand = True)[0].to_list(), columns = ['FID'])
df['ID'] = df.FID+'.p1'
df['father'] = df.FID+'.fa'
df['mother'] = df.FID+'.mo'
df['sex'] = 0
df['phenotype'] = 1

df.to_csv('Analysis_files/1_roh_samples.fam', sep = ' ', index = False, header = False)

