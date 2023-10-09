import pandas as pd
import os

# Confirm that all samples have a table with their variants

df=pd.read_csv('samples.csv')
df['var_tab']='tables/by_sample/'+df.Last_two_digits.astype(str).str.zfill(2)+'/'+df.Sample.astype(str)+'.tsv'

df['created_successfully']=df.var_tab.apply(lambda f: os.path.exists(f) and os.path.getsize(f)>0)
print(df[~df.created_successfully])

# All files exist and are not empty
