import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Plot kinship-ROH correlations in SPARK and SSC
spark=pd.read_csv('../1_SPARK/SPARK_family_table.csv')
spark['cohort']='SPARK'
spark['Sex']=spark['proband_sex']
ssc=pd.read_csv('../2_SSC/SSC_AM_family_table.csv')
ssc['cohort']='SSC'
df=pd.concat([spark[['cohort', 'spousepair', 'kinship_coeff', 'Sex', 'ROH_NSEG']], ssc[['cohort', 'spousepair', 'kinship_coeff', 'Sex', 'ROH_NSEG']]])

df['kinship_ROH']=np.nan
df.loc[(~df.kinship_coeff.isnull()) & (~df.ROH_NSEG.isnull()), 'kinship_ROH']='X'
df=df[df.kinship_ROH=='X']

df['Sex']=df.Sex.map({'Male':'Male', 'Female':'Female', 'male':'Male', 'female':'Female'})

print(df)

# Plot
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6.5, 3), sharey=True)
stat_lst=[]
for i, cohort in enumerate(['SSC', 'SPARK']):
    subdf=df[df.cohort==cohort]
    sns.scatterplot(data = subdf, y = 'kinship_coeff', x = 'ROH_NSEG', hue = 'Sex', palette = ['#FF8360', '#A1C5EE'], hue_order = ['Female', 'Male'], size = 'Sex', sizes=[30, 30], alpha = 0.5, ax=axs[i])
    # Do stats and add best-fit line
    colors=['#FF8360', '#A1C5EE']
    sexes=['Female', 'Male']
    for j in range(2):
        sns.regplot(data = subdf[subdf.Sex==sexes[j]], y = 'kinship_coeff', x = 'ROH_NSEG', color=colors[j], scatter=False, ax=axs[i])
        y = subdf[(subdf.Sex==sexes[j])]['kinship_coeff'].to_numpy()
        x = subdf[(subdf.Sex==sexes[j])]['ROH_NSEG'].to_numpy()
        res = stats.pearsonr(x, y)
        r=res.statistic
        p=res.pvalue
        ci=res.confidence_interval(confidence_level=0.95)
        stat_lst.append([sexes[j], cohort, 'ROH NSEG', r, ci.low, ci.high, p, subdf[(subdf.Sex==sexes[j])].shape[0]])
    axs[i].set_title(cohort)

plt.savefig('Figures/1_kinship_roh.pdf')

stat_df = pd.DataFrame(stat_lst, columns = ['Sex', 'Cohort', 'Attribute', 'Pearson_R', 'ci_low', 'ci_high', 'p', 'n'])
stat_df.to_csv('Result_tables/1_kinship_roh.csv', index=False)