import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Check the correlation of the euclidean distance on the first two ancestry PCs in SPARK and SSC pairs and their kinship coefficient

cohorts=['SPARK', 'SSC']

df=pd.read_csv('Analysis_files/2_all_spouse_table.csv')
df=df[df.cohort.isin(cohorts)]

# Group PCs by spouses
df=df[['Sample', 'PC1', 'PC2', 'spousepair', 'Sex', 'cohort']]
df.sort_values(by='spousepair', inplace=True)

mdf=df[df.Sex=='Male']
fdf=df[df.Sex=='Female']

mdf=mdf[mdf.spousepair.isin(fdf.spousepair.to_list())]
fdf=fdf[fdf.spousepair.isin(mdf.spousepair.to_list())]

# Calculate Euclidean distance between spouses on the first two PCs
dist=np.linalg.norm(mdf[['PC1', 'PC2']].values - fdf[['PC1', 'PC2']].values, axis=1)

df=pd.merge(fdf, mdf, on=['spousepair', 'cohort'], suffixes=['_female', '_male'])
df['euc_dist']=dist

# Annotate with kinship
spark=pd.read_csv('../1_SPARK/SPARK_family_table.csv')
# Rename spousepairs in SPARK to match PC dataframe
spark['testpair']=spark.Father+'.'+spark.Mother
spark.loc[(~spark.spousepair.isin(df.spousepair.to_list())), 'spousepair']=spark.loc[(~spark.spousepair.isin(df.spousepair.to_list())), 'testpair']
ssc=pd.read_csv('../2_SSC/SSC_AM_family_table.csv')
# Rename spousepairs in SSC to match PC dataframe
ssc['spousepair']=ssc.FID.astype(str)+'.fa.'+ssc.FID.astype(str)+'.mo'
kin=pd.concat([spark[['spousepair', 'kinship_coeff']], ssc[['spousepair', 'kinship_coeff']]])
kin=kin[kin.kinship_coeff>-0.1]

df=pd.merge(df, kin, on='spousepair', how='inner')
print(df)

colors = ['#4c2a85', '#2F9C95']

pdf=PdfPages('Figures/4_PC_kinship_correlations.pdf')
stat_lst=[]
fig, axs = plt.subplots(nrows = 2, figsize = (3, 5.5))
for i in range(2):
	ax = axs[i]

	sub_df = df[df.cohort==cohorts[i]]
	sns.scatterplot(data = sub_df, x = 'euc_dist', y = 'kinship_coeff', color = colors[i], ax = ax, size = 'cohort', sizes=[30], alpha=0.5)
	sns.regplot(data = sub_df, x = 'euc_dist', y = 'kinship_coeff', color = 'k', scatter = False, ax = ax)
	ax.set_title(cohorts[i])
	
	# Do stats
	# Pearson R
	x = sub_df['euc_dist'].to_numpy()
	y = sub_df['kinship_coeff'].to_numpy()
	res = stats.pearsonr(x, y)
	r=res.statistic
	p=res.pvalue
	ci=res.confidence_interval(confidence_level=0.95)
	stat_lst.append([cohorts[i], 'PC euclidean distance vs. kinship', 'Pearson R', len(x), r, ci.low, ci.high, p])
	
	# Annotate stats
	s='R=%.3f\np=%.2E' % (r, p)
	ax.text(0.5, 0.75, s, transform=ax.transAxes)
	
plt.tight_layout()
pdf.savefig()
plt.close()
pdf.close()

stat_df = pd.DataFrame(stat_lst, columns = ['cohort', 'comparison', 'test', 'sample_size', 'test_stat', 'ci_low', 'ci_high', 'p'])
print(stat_df)
stat_df.to_csv('Result_tables/4_PC_kinship_correlations_SPARK_SSC.csv', index=False)