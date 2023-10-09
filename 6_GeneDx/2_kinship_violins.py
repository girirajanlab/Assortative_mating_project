import pandas as pd

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# Make violin plots of raw kinship values in each group
df = pd.read_csv('Analysis_files/1_GeneDx_formatted_data.csv')
df['split'] = 1
print(df)

stat_lst=[]

pdf = PdfPages("Figures/2_kinship_violins.pdf")

# Syndromic vs. Variably Expressive
f1 = plt.figure(figsize = (3, 6))
sns.violinplot(y='Parental Kinship Values', x='split', hue='Expression', data = df, hue_order=['syndromic', 'variably expressive'], palette="Set2", inner="quartiles", cut=0)
plt.title('Syndromic vs. Variably-expressive')
plt.xlabel(None)
pdf.savefig(f1, bbox_inches='tight')
plt.close()

# Inherited vs. De novo
f1 = plt.figure(figsize = (3, 6))
sns.violinplot(y='Parental Kinship Values', x='split', hue='Inherited', data = df, hue_order = ['de novo', 'inherited'], inner="quartiles", cut=0)
plt.title('De novo vs. Inherited')
plt.xlabel(None)
pdf.savefig(f1, bbox_inches='tight')
plt.close()

# DEL vs. DUP
f1 = plt.figure(figsize = (3, 6))
sns.violinplot(y='Parental Kinship Values', x='split', hue='Del/Dup', data = df, hue_order = ['del', 'dup'], palette=["#3470BC", "#A9373B"], inner="quartiles", cut=0)
plt.title('Del vs. Dup')
plt.xlabel(None)
pdf.savefig(f1, bbox_inches='tight')
plt.close()

pdf.close()

# Statistics
def ttest(var, g1, g2):
    x=df[df[var]==g1]['Parental Kinship Values'].to_numpy()
    y=df[df[var]==g2]['Parental Kinship Values'].to_numpy()
    
    res=stats.ttest_ind(x, y, alternative='less')
    
    return([var, g1, g2, 'One tailed t-test', res.statistic, res.pvalue, len(x), len(y)])

stat_lst=[]
stat_lst.append(ttest('Expression', 'syndromic', 'variably expressive'))
stat_lst.append(ttest('Inherited', 'de novo', 'inherited'))
stat_lst.append(ttest('Del/Dup', 'del', 'dup'))

stat_df=pd.DataFrame(stat_lst, columns=['Variable', 'Group1', 'Group2', 'Test', 'statistic', 'p', 'samplesize1', 'samplesize2'])
print(stat_df)

stat_df.to_csv('Result_tables/2_kinship_violins.csv', index=False)

