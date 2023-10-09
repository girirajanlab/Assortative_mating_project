import pandas as pd
import numpy as np

fam_df=pd.read_csv('Analysis_files/1_family_table.csv')
indiv_tab=pd.read_csv('Analysis_files/SPARK_burden.csv')
indiv_tab.index=indiv_tab.Sample.to_list()

# Annotate families with whether they can be used in SNV analyses
fam_df['Family_SNV_analysis']=np.nan
fam_df.loc[(fam_df.Mother.isin(indiv_tab.Sample.to_list())) & (fam_df.Father.isin(indiv_tab.Sample.to_list())) & (fam_df.proband.isin(indiv_tab.Sample.to_list())), 'Family_SNV_analysis']='X'

fam_df['burden_correlation']=np.nan
corr_table=pd.read_csv('../5_MultiCohort/Analysis_files/2_all_spouse_table.csv')
used_pairs=corr_table.spousepair.value_counts()
used_pairs=used_pairs[used_pairs==2]
fam_df['testpair']=fam_df.Father+'.'+fam_df.Mother
fam_df.loc[((fam_df.spousepair.isin(used_pairs.index.to_list())) | (fam_df.testpair.isin(used_pairs.index.to_list()))), 'burden_correlation']='X'
fam_df.loc[(fam_df.Mother.isnull()) | (fam_df.Father.isnull()), 'burden_correlation']=np.nan

# Add proband ROH to table
roh = pd.read_csv('Analysis_files/ROH_filtered_summary.csv')
roh.index=roh.Sample.to_list()
fam_df['ROH_NSEG']=fam_df.proband.map(roh.NSEG.to_dict())

fam_df['ROH_analysis']=np.nan
fam_df.loc[(~fam_df.ROH_NSEG.isnull()) & (~fam_df.kinship_coeff.isnull()), 'ROH_analysis']='X'

fam_df['Used']=np.nan
fam_df.loc[fam_df.Parent_child_pheno=='X', 'Used']='X'
fam_df.loc[fam_df.Spouse_pheno=='X', 'Used']='X'
fam_df.loc[fam_df.Family_SNV_analysis=='X', 'Used']='X'
fam_df.loc[fam_df.ROH_analysis=='X', 'Used']='X'
fam_df.loc[fam_df.burden_correlation=='X', 'Used']='X'

# Also annotate with whether family is used as trio, parent-child pair, or spouse pair
fam_df['type']=''
fam_df.loc[(fam_df.Spouse_pheno=='X') | (fam_df.burden_correlation=='X'), 'type']='spousepair'
fam_df.loc[(fam_df.Parent_child_pheno=='X') & (fam_df.type!='spousepair'), 'type']='parent-child pair'
fam_df.loc[(fam_df.ROH_analysis=='X') | (fam_df.Family_SNV_analysis=='X'), 'type']='trio'
fam_df.loc[(fam_df.Parent_child_pheno=='X') & (fam_df.type=='spousepair'), 'type']='trio'
print(fam_df[fam_df.Used=='X'].type.value_counts())


# Save
fam_df.to_csv('SPARK_family_table.csv', index=False)