import pandas as pd
import numpy as np

# Add ROH information to family dataframe
df=pd.read_csv('Analysis_files/2_SSC_data.csv')

roh=pd.read_csv('Analysis_files/ROH_summary.csv')
roh.index=roh.Sample.to_list()

df['ROH_NSEG']=df.proband.map(roh.NSEG.to_dict())

# Annotate families used in ROH analysis
df['kinship_ROH_analysis']=np.nan
df.loc[(~df.ROH_NSEG.isnull()) & (~df.kinship_coeff.isnull()) & (df.kinship_coeff>-0.1), 'kinship_ROH_analysis']='X'

# Also annotate families used in burden correlation analyses
corr_table=pd.read_csv('../5_MultiCohort/Analysis_files/2_all_spouse_table.csv')
used_pairs=corr_table.spousepair.value_counts()
used_pairs=used_pairs[used_pairs==2]
df['testpair']=df.Father+'.'+df.Mother
df.loc[((df.spousepair.isin(used_pairs.index.to_list())) | (df.testpair.isin(used_pairs.index.to_list()))), 'burden_correlation']='X'

# Annotate families with those used in any analysis
df['Used']=np.nan
df.loc[(df.linear_models=='X') | (df.Parent_phenotype_analysis=='X') | (df.kinship_ROH_analysis=='X') | (df.kinship_SFARI_analysis=='X') | (df.burden_correlation=='X'), 'Used']='X'
print(df.Used.value_counts())

# Also annotate type of family (trio or spouse pair)
df['type']=''
df.loc[(df.Parent_phenotype_analysis=='X') | (df.burden_correlation=='X'), 'type']='spouse pair'
df.loc[(df.linear_models=='X') | (df.kinship_ROH_analysis=='X') | (df.kinship_SFARI_analysis=='X'), 'type']='trio'

print(df.type.value_counts())

# Save
df.to_csv('SSC_AM_family_table.csv', index=False)
