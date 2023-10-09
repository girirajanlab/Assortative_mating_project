import pandas as pd

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Examine correlations in PCs between spouses in SSC and SPARK

# Concat PCs from all cohorts
pc_16p=pd.read_csv('Analysis_files/16p12_snv_table.csv')
pc_16p['cohort']='16p12.1 deletion'

spouse_16p=pd.read_csv('../3_16p12.1/16p12.1_family_table.csv')
sample=spouse_16p.Mother.to_list()+spouse_16p.Father.to_list()
partner=spouse_16p.Father.to_list()+spouse_16p.Mother.to_list()
sexes=['Female']*spouse_16p.shape[0]+['Male']*spouse_16p.shape[0]

spark=pd.read_csv('Analysis_files/SPARK_snv_table.csv')
spark['cohort']='SPARK'

spouse_spark=pd.read_csv('../1_SPARK/SPARK_family_table.csv')
sample+=spouse_spark.Mother.astype(str).to_list()+spouse_spark.Father.astype(str).to_list()
partner+=spouse_spark.Father.astype(str).to_list()+spouse_spark.Mother.astype(str).to_list()
sexes+=['Female']*spouse_spark.shape[0]+['Male']*spouse_spark.shape[0]

ssc=pd.read_csv('Analysis_files/SSC_snv_table.csv')
ssc['cohort']='SSC'

spouse_ssc=pd.read_csv('../2_SSC/SSC_AM_family_table.csv')
sample+=spouse_ssc.Mother.astype(str).to_list()+spouse_ssc.Father.astype(str).to_list()
partner+=spouse_ssc.Father.astype(str).to_list()+spouse_ssc.Mother.astype(str).to_list()
sexes+=['Female']*spouse_ssc.shape[0]+['Male']*spouse_ssc.shape[0]

ukb=pd.read_csv('Analysis_files/UKB_snv_table.csv')
ukb['cohort']='UK Biobank'

spouse_ukb=pd.read_csv('../4_UK Biobank/2_spouse_table/UKB_family_table.csv')
sample+=spouse_ukb.Female.astype(str).to_list()+spouse_ukb.Male.astype(str).to_list()
partner+=spouse_ukb.Male.astype(str).to_list()+spouse_ukb.Female.astype(str).to_list()
sexes+=['Female']*spouse_ukb.shape[0]+['Male']*spouse_ukb.shape[0]

# Merge
df=pd.concat([pc_16p, spark])
df=pd.concat([df, ssc])
df=pd.concat([df, ukb])

spouse_dict=dict(zip(sample, partner))
df['spouse']=df.Sample.astype(str).map(spouse_dict)

spousepairs=dict(zip(sample, ['.'.join(sorted([i, spouse_dict[i]])) for i in sample]))
df['spousepair']=df.Sample.astype(str).map(spousepairs)

df['Sex']=df.Sample.astype(str).map(dict(zip(sample, sexes)))

# Save this file for other analyses
df.to_csv('Analysis_files/2_all_spouse_table.csv', index=False)

# Group PCs by spouses
df=df[['Sample', 'PC1', 'PC2', 'spousepair', 'Sex', 'cohort']]

mdf=df[df.Sex=='Male']
fdf=df[df.Sex=='Female']

df=pd.merge(fdf, mdf, on=['spousepair', 'cohort'], suffixes=['_female', '_male'])
df=df[(~df.Sample_female.isnull()) & (~df.Sample_male.isnull())]

cohorts=['SPARK', 'SSC']
colors = ['#4c2a85', '#2F9C95']

pdf=PdfPages('Figures/2_PC_correlations_SPARK_SSC.pdf')
stat_lst=[]
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (6, 5.5))
col=0
for pc in ['PC1', 'PC2']:
	for i in range(2):
		ax = axs[i, col]

		sub_df = df[df.cohort==cohorts[i]]
		sns.scatterplot(data = sub_df, x = pc+'_male', y = pc+'_female', color = colors[i], ax = ax, size = 'cohort', sizes=[30], alpha=0.5)
		sns.regplot(data = sub_df, x = pc+'_male', y = pc+'_female', color = 'k', scatter = False, ax = ax)
		ax.set_title(cohorts[i]+' '+pc)
		lox, hix=ax.get_xlim()
		loy, hiy=ax.get_ylim()
		hi=max([hix, hiy])
		lo=min([lox, loy])
		ax.set_xlim(lo, hi)
		ax.set_ylim(lo, hi)
		
		# Do stats
		# Pearson R
		x = sub_df[pc+'_male'].to_numpy()
		y = sub_df[pc+'_female'].to_numpy()
		res = stats.pearsonr(x, y)
		r=res.statistic
		p=res.pvalue
		ci=res.confidence_interval(confidence_level=0.95)
		stat_lst.append([cohorts[i], pc, 'Pearson R', len(x), r, ci.low, ci.high, p])
	col+=1
	
plt.tight_layout()
pdf.savefig()
plt.close()
pdf.close()

stat_df = pd.DataFrame(stat_lst, columns = ['cohort', 'comparison', 'test', 'sample_size', 'test_stat', 'ci_low', 'ci_high', 'p'])
print(stat_df)
stat_df.to_csv('Result_tables/2_PC_correlations_SPARK_SSC.csv', index=False)