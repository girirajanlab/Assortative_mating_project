import pandas as pd

from statsmodels.stats.proportion import proportions_ztest

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns

# Plot the proportion of children developing phenotypes based on whether 0, 1, or both of their parents have the phenotype
# Also separate by sex
phenos=pd.read_csv('Analysis_files/1_family_phenotypes.csv')
phenos=phenos[(phenos.parent_child=='X') & (phenos.spouse=='X')]

# Skip personality disorder due to low sample size in probands
phenotypes=['behav_adhd', 'mood_anx', 'mood_bipol', 'mood_dep', 'mood_ocd', 'sleep_dx']
for p in phenotypes:
    phenos[p+'_parent']=phenos[[p+'_mother', p+'_father']].sum(axis=1)

# Get proportion of children with phenotypes in each group
dat=[]
for sex in ['Male', 'Female']:
    for np in [0, 1, 2]:
        for p in phenotypes:
            subdf=phenos[(phenos.sex==sex) & (phenos[p+'_parent']==np) & (~phenos[p+'_proband'].isnull())]
            w=sum(subdf[p+'_proband'].to_list())
            perc=(w*100)/subdf.shape[0]
            
            dat.append([sex, np, p, int(w), subdf.shape[0], perc])

props=pd.DataFrame(dat, columns=['Sex', 'number_of_parents', 'phenotype', 'number_with', 'total', 'percent'])
# Save
props.to_csv('Result_tables/6_phenotype_proportions_by_assortment.csv', index=False)

# Two sample proportion tests
stats=[]
for p in phenotypes:
    for sex in ['Male', 'Female']:
        for g1 in [0, 1, 2]:
            for g2 in [0, 1, 2]:
                if g1>=g2:
                    continue
                c1=props[(props.phenotype==p) & (props.number_of_parents==g1) & (props.Sex==sex)]['number_with'].to_list()[0]
                c2=props[(props.phenotype==p) & (props.number_of_parents==g2) & (props.Sex==sex)]['number_with'].to_list()[0]
                t1=props[(props.phenotype==p) & (props.number_of_parents==g1) & (props.Sex==sex)]['total'].to_list()[0]
                t2=props[(props.phenotype==p) & (props.number_of_parents==g2) & (props.Sex==sex)]['total'].to_list()[0]
                res=proportions_ztest(count=[c1, c2], nobs=[t1, t2], alternative='smaller')
                stats.append([p, sex, g1, g2, c1, t1, c2, t2, res[0], res[1]])
stats_df=pd.DataFrame(stats, columns=['phenotype', 'sex', 'group1', 'group2', 'n_with_group1', 'total_group1', 'n_with_group2', 'total_group2', 'statistic', 'p'])
stats_df['bonferroni_p']=stats_df.p*stats_df.shape[0]
stats_df.loc[stats_df.bonferroni_p>1, 'bonferroni_p']=1
# Save to file
stats_df.to_csv('Result_tables/6_proportion_tests.csv', index=False)
stats_df['star']='n.s.'
stats_df.loc[stats_df.bonferroni_p<=0.05, 'star']='*'
stats_df.loc[stats_df.bonferroni_p<=0.01, 'star']='**'
stats_df.loc[stats_df.bonferroni_p<=0.001, 'star']='***'

# Make barplots showing data
# Plot each sex on a separate graph
fig, axs=plt.subplots(ncols=2, figsize=(15, 6), sharey=True)
props['x']=props.phenotype.map({'behav_adhd':0, 'mood_anx':1, 'mood_bipol':2, 'mood_dep':3, 'mood_ocd':4, 'sleep_dx':5})
for i, sex in enumerate(['Male', 'Female']):
    subdf=props[props.Sex==sex]
    sns.barplot(subdf, x='x', y='percent', hue='number_of_parents', palette='viridis_r', ax=axs[i])
    
    # Annotate signifiance on plot
    for x, p in enumerate(phenotypes):
        max_val=max(props[(props.phenotype==p) & (props.Sex==sex)]['percent'])
        axs[i].plot([-0.3+x, 0+x], [max_val+2, max_val+2], color='k', lw=1)
        axs[i].plot([0+x, 0.3+x], [max_val+8, max_val+8], color='k', lw=1)
        axs[i].plot([-0.3+x, 0.3+x], [max_val+5, max_val+5], color='k', lw=1)
        
        axs[i].text(-0.15+x, max_val+3, stats_df[(stats_df.phenotype==p) & (stats_df.sex==sex) & (stats_df.group1==0) & (stats_df.group2==1)]['star'].to_list()[0], ha='center')
        axs[i].text(0.15+x, max_val+9, stats_df[(stats_df.phenotype==p) & (stats_df.sex==sex) & (stats_df.group1==1) & (stats_df.group2==2)]['star'].to_list()[0], ha='center')
        axs[i].text(0+x, max_val+6, stats_df[(stats_df.phenotype==p) & (stats_df.sex==sex) & (stats_df.group1==0) & (stats_df.group2==2)]['star'].to_list()[0], ha='center')
    
    # Clean up
    axs[i].set_xticks([0, 1, 2, 3, 4, 5], ['ADHD', 'Anxiety', 'Bipolar\ndisorder', 'Depression', 'OCD', 'Sleep\ndisorder'])
    axs[i].set_xlabel('')
    axs[i].set_title(sex)
    if i==0:
        axs[i].set_ylabel('Percentage of probands with phenotype')
    else:
        axs[i].set_ylabel('')
    
plt.suptitle('Proband phenotypes by parent assortment')
plt.tight_layout()
plt.savefig('Figures/6_phenotype_by_assortment.pdf')
