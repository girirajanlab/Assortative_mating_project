import pandas as pd

# Get burden counts for all participants in SPARK
df=pd.read_csv('../WES1/tables/11_gene_ids.csv')

# Get burden by sample
burden=pd.DataFrame(df.Sample.value_counts())
burden.columns=['burden']
burden['Sample']=burden.index.to_list()

# Save to file
burden.to_csv('analysis_files/12_sample_burden.csv', index=False)
