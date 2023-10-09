import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42

df=pd.read_csv('SSC_AM_family_table.csv')
df=df[df.kinship_SFARI_analysis=='X']

# Compare kinship values of probands with LOF mutations in Tier S SFARI genes to those without SFARI mutations
df['label']='No Tier S SFARI SNV'
df.loc[df['Tier S SFARI SNV']!='0.0', 'label']='Tier S SFARI SNV'

f = plt.figure(figsize = (2, 6))
sns.violinplot(data=df, x='kinship_SFARI_analysis', y='kinship_coeff', hue='label', hue_order=['Tier S SFARI SNV', 'No Tier S SFARI SNV'], palette=['#788B97', '#CC5C6D'], inner='quartiles', cut=0)
plt.savefig('Figures/5_kinship_sfari.pdf')

x=df[df.label=='Tier S SFARI SNV']['kinship_coeff'].to_numpy()
y=df[df.label=='No Tier S SFARI SNV']['kinship_coeff'].to_numpy()
t, p = stats.ttest_ind(x, y, alternative='less')

with open('Result_tables/5_kinship_sfari.csv', 'w') as outfile:
    outfile.write(','.join(['Test', 'Group1', 'Group2', 'statistic', 'p', 'samplesize_group1', 'samplesize_group2'])+'\n')
    outfile.write(','.join(['One tailed t test', 'Probands with SNVs in Tier S SFARI genes', 'Probands without SNVs in Tier S SFARI genes', str(t), str(p), str(len(x)), str(len(y))])+'\n')
