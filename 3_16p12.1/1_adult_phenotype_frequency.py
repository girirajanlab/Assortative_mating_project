import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42

# Cohort information
df = pd.read_excel('Analysis_files/Table_S1_Cohort_description.xlsx')
df=df[df.Relationship!='Proband']

# For phenotype frequency, restrict to complete spouse pairs
df=df[~df.Spouse.isnull()]
c_dict=dict(zip(df.Sample.to_list(), df['16p12.1 deletion carrier'].to_list()))
df['spouse_status']=df.Spouse.map(c_dict)
freqdf=df[((df['16p12.1 deletion carrier']=='Carrier') & (df.spouse_status=='Noncarrier')) | ((df['16p12.1 deletion carrier']=='Noncarrier') & (df.spouse_status=='Carrier')) | ((df['16p12.1 deletion carrier']=='Carrier') & (df.spouse_status.isnull()))].copy()
freqdf.reset_index(drop=True, inplace=True)
freqdf.to_csv('Analysis_files/1_phenotype_samples.csv', index=False)

# Get phenotype frequencies in adults
plot_df=pd.DataFrame({'phenotype':['Seizures', 'Schizophrenic Features', 'Depression', 'Anxiety', 'Addiction']*2,
                        'carrier':['Carrier']*5+['Noncarrier']*5, 'Count':[0]*10, 'Freq':[0]*10, 'Sample_size':[0]*10})
colors=['#EFAF49', '#007A8D']

for pheno in ['Seizures', 'Schizophrenic Features', 'Depression', 'Anxiety', 'Addiction']:
    for c in ['Carrier', 'Noncarrier']:
        subdf=freqdf[(freqdf['16p12.1 deletion carrier']==c) & (~freqdf[pheno].isnull())].copy()
        plot_df.loc[(plot_df.phenotype==pheno) & (plot_df.carrier==c), 'Count']=subdf[subdf[pheno]>0].shape[0]
        plot_df.loc[(plot_df.phenotype==pheno) & (plot_df.carrier==c), 'Freq']=subdf[subdf[pheno]>0].shape[0]*100/subdf.shape[0]
        plot_df.loc[(plot_df.phenotype==pheno) & (plot_df.carrier==c), 'Sample_size']=subdf.shape[0]

# Plot carriers and noncarriers
sns.barplot(data=plot_df, x='phenotype', y='Freq', hue='carrier', hue_order=['Carrier', 'Noncarrier'], palette=['#EFAF49', '#007A8D'])
plt.savefig('Figures/1_phenotype_bars.pdf')

# Save counts to file
plot_df.to_csv('Result_tables/1_phenotype_frequency.csv', index=False)

# Reformat data for spouse correlations
# Duplicate samples with multiple spouses
for r in df[df.Spouse.str.contains('/')].index.to_list():
    newrow=df.loc[r].to_frame().T
    newrow['Spouse']=newrow.Spouse.str.split('/', expand=True)[1]
    df=pd.concat([df, newrow])
df.reset_index(drop=True, inplace=True)
df.Spouse=df.Spouse.str.split('/', expand=True)[0]
df['spousepair']=df.Sample+'.'+df.Spouse
df.loc[df.Sex=='Male', 'spousepair']=df.Spouse+'.'+df.Sample

# Combine spouse data together in a single line
cols=['Family', 'spousepair', 'Sample', '16p12.1 deletion carrier', 'Seizures', 'Schizophrenic Features', 'Depression', 'Anxiety', 'Addiction']
spouse_df=pd.merge(df[df.Sex=='Male'][cols], df[df.Sex=='Female'][cols], on=['Family', 'spousepair'], suffixes=['_male', '_female'])
print(spouse_df[(spouse_df['16p12.1 deletion carrier_male']=='Carrier') & (spouse_df['16p12.1 deletion carrier_female']=='Carrier')])
spouse_df['Carrier']='None'
spouse_df.loc[spouse_df['16p12.1 deletion carrier_male']=='Carrier', 'Carrier']='Male'
spouse_df.loc[spouse_df['16p12.1 deletion carrier_female']=='Carrier', 'Carrier']='Female'

cols=[]
for c in ['Carrier', 'Noncarrier']:	
	for p in ['Seizures', 'Schizophrenic Features', 'Depression', 'Anxiety', 'Addiction']:
		spouse_df[c+'.'+p]=np.nan
		spouse_df.loc[spouse_df['16p12.1 deletion carrier_male']==c, c+'.'+p]=spouse_df.loc[spouse_df['16p12.1 deletion carrier_male']==c, p+'_male']
		spouse_df.loc[spouse_df['16p12.1 deletion carrier_female']==c, c+'.'+p]=spouse_df.loc[spouse_df['16p12.1 deletion carrier_female']==c, p+'_female']
		if c=='Noncarrier':
			spouse_df.loc[spouse_df.Carrier=='None', c+'.'+p]=np.nan
		cols.append(c+'.'+p)

# Annotate pairs used in phenotpyic analyses
spouse_df['carrier_phenotype']=np.nan
for p in [i for i in cols if 'Carrier.' in i]:
    spouse_df.loc[~spouse_df[p].isnull(), 'carrier_phenotype']='X'
spouse_df['noncarrier_phenotype']=np.nan
for p in [i for i in cols if 'Noncarrier.' in i]:
    spouse_df.loc[~spouse_df[p].isnull(), 'noncarrier_phenotype']='X'
spouse_df['spouse_phenotypes']=np.nan
spouse_df.loc[(spouse_df.carrier_phenotype=='X') & (spouse_df.noncarrier_phenotype=='X'), 'spouse_phenotypes']='X'

# Save to file
spouse_df[['Family', 'spousepair', 'Sample_male', 'Sample_female', 'Carrier']+cols+['carrier_phenotype', 'noncarrier_phenotype', 'spouse_phenotypes']].to_csv('Analysis_files/1_spouse_phenotypes.csv', index=False)
