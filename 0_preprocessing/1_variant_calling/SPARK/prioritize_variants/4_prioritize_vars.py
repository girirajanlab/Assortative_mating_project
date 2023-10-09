import pandas as pd

# Variant IDs are from 1_SPARK/7_geneset_examples.py
int_vars=pd.read_csv('Analysis_files/7_variant_ids.csv')
vars=pd.read_csv('tables/2_slim_variants.csv')
vars['vid']=vars.Chr+'_'+vars.Pos.astype(str)+'_'+vars.Ref+'_'+vars.Alt

vars=vars[vars.vid.isin(int_vars.vid.to_list())]
print(vars)

vars.to_csv('tables/4_var_samples.csv', index=False)
