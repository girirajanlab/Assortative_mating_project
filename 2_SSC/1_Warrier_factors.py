import pandas as pd

# Calculate the factors from Warrier et al. Nature Genetics 2022: https://www.nature.com/articles/s41588-022-01072-5#Sec10
# Questions used in factors are found in Supp. Table 3
# Supp Table 3 is copied in Analysis_files/Warrier_tabS3.csv
qs=pd.read_csv('Analysis_files/Warrier_tabS3.csv')
qs['Factor']=qs.Factor.str.split('\\. ', expand=True)[1]
factors=list(qs.Factor.unique())

# For each "factor" calculate mean score of all questions
df=pd.DataFrame({'individual':['0']})
for f in factors:
    files=list(qs[qs.Factor==f]['SSC_filename'].unique())
    cols=list(qs[qs.Factor==f]['Items'].unique())
    samp_df=pd.DataFrame({'individual':['0']})
    for fi in files:
        fdf=pd.read_csv(fi)
        if 'scq' in fi:
            rel_cols=['q'+i[1:] for i in cols if 'q'+i[1:] in fdf.columns.to_list()]
        else:
            rel_cols=[i for i in cols if i in fdf.columns.to_list()]
        fdf=fdf[['individual']+rel_cols]
        samp_df=pd.merge(samp_df, fdf, on='individual', how='outer')
    samp_df.dropna(how='any', axis=0, inplace=True)
    sum_cols=[i for i in samp_df.columns.to_list() if i!='individual']
    samp_df[f]=samp_df[sum_cols].sum(axis=1)/len(sum_cols)
    df=pd.merge(df, samp_df[['individual', f]], on='individual', how='outer')
        
df=df[df.individual!='0']
# Save to file
df.to_csv('Analysis_files/1_ASD_factors.csv', index=False)