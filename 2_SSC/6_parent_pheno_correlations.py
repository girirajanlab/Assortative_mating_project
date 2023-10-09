import pandas as pd

import scipy.stats as stats
import statsmodels.api as sm

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Check correlations of parental SRS and BAPQ scores
df=pd.read_csv('SSC_AM_family_table.csv')

pdf=PdfPages('Figures/6_parent_quantitative_correlations.pdf')

# Regress out top 10 PCs
modf=df[['Mother', 'Mother_SRS', 'Mother_BAPQ']]
fadf=df[['Father', 'Father_SRS', 'Father_BAPQ']]
modf.columns=['Sample', 'SRS', 'BAPQ']
fadf.columns=modf.columns
df2=pd.concat([modf, fadf])

reg=pd.read_csv('Analysis_files/SSC_snv_table.csv')
df2=pd.merge(df2, reg, on='Sample', how='inner')

pred_cols=['PC'+str(i) for i in range(1, 11)]
for pc in pred_cols:
	df2=df2[~df2[pc].isnull()]
	df2[pc]=(df2[pc]-df2[pc].mean())/df2[pc].std()

for t in ['SRS', 'BAPQ']:
	testdf=df2[~df2[t].isnull()]
	testdf[t]==(testdf[pc]-testdf[pc].mean())/testdf[pc].std()
	
	y=testdf[t].to_numpy()
	x=testdf[pred_cols].to_numpy()
	X=sm.add_constant(x)
	
	mod=sm.OLS(y, X)
	res=mod.fit()
	print(res.summary())
	testdf[t+'resid']=res.resid
	df2[t+'resid']=df2.Sample.map(dict(zip(testdf.Sample.to_list(), testdf[t+'resid'].to_list())))

df2['Family']=df2.Sample.str.split('.', expand=True)[0].astype(int)

residdf=pd.merge(df2[df2.Sample.str.contains('.mo')][['Sample', 'Family', 'SRSresid', 'BAPQresid']], df2[df2.Sample.str.contains('.fa')][['Sample', 'Family', 'SRSresid', 'BAPQresid']], on='Family', how='inner', suffixes=['_mother', '_father'])

stat_lst=[]
fig, axs = plt.subplots(ncols=2, figsize=(9, 4))
for i, t in enumerate(['SRSresid', 'BAPQresid']):
	plotdf=residdf[(~residdf[t+'_mother'].isnull()) & (~residdf[t+'_father'].isnull())]
	res=stats.pearsonr(plotdf[t+'_mother'].to_numpy(), plotdf[t+'_father'].to_numpy())
	ci=res.confidence_interval(confidence_level=0.95)
	stat_lst.append([t+' - 10 PCs', res.statistic, ci.low, ci.high, res.pvalue, plotdf.shape[0]])
	
	# Plot correlations
	sns.scatterplot(plotdf, y=t+'_mother', x=t+'_father', color='#2F9C95', ax=axs[i])
	sns.regplot(plotdf, y=t+'_mother', x=t+'_father', color='k', scatter=False, ax=axs[i])
	
	# Annotate values
	txt='R=%.2f\np=%.2E' % (res.statistic, res.pvalue)
	axs[i].annotate(txt, xy=(0.02, 0.85), xycoords='axes fraction')
	
	axs[i].set_title(t+' - 10 PCs')
pdf.savefig()
plt.close()

pdf.close()

stat_df=pd.DataFrame(stat_lst, columns=['phenotype', 'Pearson_R', 'CI_lower', 'CI_upper', 'p', 'n'])
print(stat_df)
stat_df.to_csv('Result_tables/6_parent_quantitative_correlations.csv', index=False)
	