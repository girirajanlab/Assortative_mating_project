import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

import scipy.stats as stats

# Assess the changes in SNV burden over generations
cohort=pd.read_excel('Analysis_files/Table_S1_Cohort_description.xlsx')

# Get probands, parents, and grandparents (parents of carrier parent)
pro_pars=cohort[cohort.Relationship=='Proband']['Mother'].to_list()+cohort[cohort.Relationship=='Proband']['Father'].to_list()
cpar_pars=cohort[(cohort.Sample.isin(pro_pars)) & (cohort['16p12.1 deletion carrier']=='Carrier')]['Mother'].to_list()+cohort[(cohort.Sample.isin(pro_pars)) & (cohort['16p12.1 deletion carrier']=='Carrier')]['Father'].to_list()
cohort=cohort[(cohort.Relationship=='Proband') | (cohort.Sample.isin(pro_pars+cpar_pars))]

# Annotate with generation
cohort['generation']=3
cohort.loc[cohort.Relationship.isin(['Mother', 'Father']), 'generation']=2
cohort.loc[cohort.Relationship.isin(['Grandmother', 'Grandfather']), 'generation']=1
print(cohort.generation.value_counts())

cohort=cohort[~cohort['SNV burden'].isnull()]

print(cohort[cohort.Sample.isin(cpar_pars)])

# Plot
sns.boxplot(data=cohort, y='SNV burden', x='generation', hue='16p12.1 deletion carrier', hue_order=['Carrier', 'Noncarrier'], palette=['#EFAF49', '#007A8D'], fliersize=0)
sns.stripplot(data=cohort, y='SNV burden', x='generation', hue='16p12.1 deletion carrier', hue_order=['Carrier', 'Noncarrier'], dodge=True, palette=['k', 'k'], marker='$\circ$', s=10, jitter=0.3)
plt.savefig('Figures/4_generation_burden.pdf')

# Perform t-tests between carriers in each generation
stat_lst=[]
for gen1 in [1, 2]:
    for gen2 in [2, 3]:
        if gen1==gen2:
            continue
        g1=cohort[(cohort['16p12.1 deletion carrier']=='Carrier') & (cohort.generation==gen1)]['SNV burden'].to_numpy()
        g2=cohort[(cohort['16p12.1 deletion carrier']=='Carrier') & (cohort.generation==gen2)]['SNV burden'].to_numpy()
        res=stats.ttest_ind(g1, g2, alternative='less')
        stat_lst.append([gen1, gen2, 'One tailed t-test', res.statistic, res.pvalue, len(g1), len(g2)])
stat_df=pd.DataFrame(stat_lst, columns=['Generation1', 'Generation2', 'Test', 'statistic', 'p', 'samplesize1', 'samplesize2'])

# Save to file
print(stat_df)
stat_df.to_csv('Result_tables/4_generation_burden.csv', index=False)

# Annotate pairs with whether they were used in this analysis
df=pd.read_csv('Analysis_files/3_genetic_data.csv')
df['generation_burden_analysis']=np.nan
gen2=cohort[cohort.generation==2]['Sample'].to_list()
df.loc[(df.Mother.isin(gen2)) & (df.Father.isin(gen2)), 'generation_burden_analysis']='P'
gen1=cohort[cohort.generation==1]['Sample'].to_list()
df.loc[(df.Mother.isin(gen1)) & (df.Father.isin(gen1)), 'generation_burden_analysis']='GP'

# Also annotate samples with whether they can be used in the burden correlations
df['burden_correlation']=np.nan
df['spouse_config']=df.Father+'.'+df.Mother
corr_table=pd.read_csv('../5_MultiCohort/Analysis_files/2_all_spouse_table.csv')
used_pairs=corr_table.spousepair.value_counts()
used_pairs=used_pairs[used_pairs==2]
df.loc[df.spouse_config.isin(used_pairs.index.to_list()), 'burden_correlation']='X'

df['Used']=np.nan
df.loc[(df.spouse_phenotypes=='X') | (df.generation_burden_analysis.isin(['P', 'GP'])) | (df.burden_correlation=='X'), 'Used']='X'
print(df)
print(df.Used.value_counts())
df.sort_values(by='FID', inplace=True)

# Save
df.to_csv('16p12.1_family_table.csv', index=False)