import pandas as pd

# Get burden of variants in AM families
df=pd.read_csv('tables/4_am_private_in_family.csv')

syn=pd.DataFrame(df[df.Mut_type=='synonymous']['Sample'].value_counts())
syn.columns=['synonymous_burden']

var=pd.DataFrame(df[df.Mut_type!='synonymous']['Sample'].value_counts())
var.columns=['variant_burden']

burden=pd.merge(syn, var, left_index=True, right_index=True)
print(burden)

# Save to file
burden.to_csv('tables/5_burden_table.csv')
