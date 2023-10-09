import pandas as pd
import statsmodels.api as sm

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Examine burden correlations between spouses in all cohorts
df=pd.read_csv('Analysis_files/2_all_spouse_table.csv')

# Regress out 10 genetic PCs and synonymous burden from rare variant burden by cohort
cohorts=['16p12.1 deletion', 'SPARK', 'SSC', 'UK Biobank']
resids=pd.DataFrame()
for cohort in cohorts:
	pred_cols=['synonymous_burden']+ ['PC'+str(i) for i in range(1, 11)]
	
	subdf=df[df.cohort==cohort].copy()
	for pc in ['variant_burden']+pred_cols:
		subdf[pc]=(subdf[pc]-subdf[pc].mean())/subdf[pc].std()

	y=subdf['variant_burden'].to_numpy()
	x=subdf[pred_cols].to_numpy()
	X=sm.add_constant(x)
	mod=sm.OLS(y, X)
	res=mod.fit()
	
	subdf['residuals']=res.resid
	
	resids=pd.concat([resids, subdf[['Sample', 'residuals']]])

df=pd.merge(df, resids, on='Sample')

# Group data by spouses
df=df[['Sample', 'residuals', 'spousepair', 'Sex', 'cohort']]

mdf=df[df.Sex=='Male']
fdf=df[df.Sex=='Female']

df=pd.merge(fdf, mdf, on=['spousepair', 'cohort'], suffixes=['_female', '_male'])
df=df[(~df.Sample_female.isnull()) & (~df.Sample_male.isnull())]

colors = ['k', '#4c2a85', '#2F9C95', '#C75146']

pdf=PdfPages('Figures/3_burden_correlations.pdf')
stat_lst=[]
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (6, 6), sharex=True, sharey=True)
for i in range(4):
	if i < 2:
		ax = axs[0, i]
	else:
		ax = axs[1, i-2]

	sub_df = df[df.cohort==cohorts[i]]
	if i==0:
		print(sub_df)
	sns.scatterplot(data = sub_df, x = 'residuals_male', y = 'residuals_female', color = colors[i], ax = ax)
	sns.regplot(data = sub_df, x = 'residuals_male', y = 'residuals_female', color = 'k', scatter = False, ax = ax)
	ax.set_title(cohorts[i])
	ax.set_xlim(-4, 8)
	ax.set_ylim(-4, 8)
	
	# Do stats
	# Pearson R
	x = sub_df['residuals_male'].to_numpy()
	y = sub_df['residuals_female'].to_numpy()
	res = stats.pearsonr(x, y)
	r=res.statistic
	p=res.pvalue
	ci=res.confidence_interval(confidence_level=0.95)
	stat_lst.append([cohorts[i], 'residual burden', 'Pearson R', len(x), r, ci.low, ci.high, p])
	
	# Annotate stats
	s='R=%.3f\np=%.3f' % (r, p)
	ax.text(-2.5, 5, s)
	
plt.tight_layout()
pdf.savefig()
plt.close()

# Also check correlations in PRS in SSC samples
ssc=pd.read_csv('../2_SSC/SSC_AM_family_table.csv')
prs_dict=dict(zip(ssc.Mother.to_list()+ssc.Father.to_list(), ssc.Mother_PRS.to_list()+ssc.Father_PRS.to_list()))
df=pd.read_csv('Analysis_files/2_all_spouse_table.csv')
df=df[df.cohort=='SSC']
df['PRS']=df.Sample.map(prs_dict)
df=df[~df.PRS.isnull()]

# Regress out 10 PCs
pred_cols=['PC'+str(i) for i in range(1, 11)]
	
for pc in ['PRS']+pred_cols:
	df[pc]=(df[pc]-df[pc].mean())/df[pc].std()

y=df['PRS'].to_numpy()
x=df[pred_cols].to_numpy()
X=sm.add_constant(x)
mod=sm.OLS(y, X)
res=mod.fit()

df['residuals']=res.resid

# Group data by spouses
df=df[['Sample', 'residuals', 'spousepair', 'Sex']]

mdf=df[df.Sex=='Male']
fdf=df[df.Sex=='Female']

df=pd.merge(fdf, mdf, on=['spousepair'], suffixes=['_female', '_male'])
df=df[(~df.residuals_female.isnull()) & (~df.residuals_male.isnull())]

fig, ax = plt.subplots(figsize=(5,5))
sns.scatterplot(data = df, x = 'residuals_male', y = 'residuals_female', color = colors[2], ax = ax)
sns.regplot(data = df, x = 'residuals_male', y = 'residuals_female', color = 'k', scatter = False, ax = ax)
plt.title('SSC PRS')

# Do stats
# Pearson R
x = df['residuals_male'].to_numpy()
y = df['residuals_female'].to_numpy()
res = stats.pearsonr(x, y)
r=res.statistic
p=res.pvalue
ci=res.confidence_interval(confidence_level=0.95)
stat_lst.append(['SSC', 'residual PRS', 'Pearson R', len(x), r, ci.low, ci.high, p])

# Annotate stats
s='R=%.3f\np=%.3f' % (r, p)
plt.text(-3, 3, s)
	
pdf.savefig()
plt.close()

pdf.close()

stat_df = pd.DataFrame(stat_lst, columns = ['cohort', 'comparison', 'test', 'sample_size', 'test_stat', 'ci_low', 'ci_high', 'p'])
print(stat_df)
stat_df.to_csv('Result_tables/3_burden_correlations.csv', index=False)