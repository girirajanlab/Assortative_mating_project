import pandas as pd

# Get the number and proportion of each phenotype in the assessed pairs in UK Biobank
df=pd.read_csv('../2_spouse_table/UKB_family_table.csv')
df=df[df.phenotype_analysis=='X']

phenos=['Anxiety', 'BPD_mania', 'Depression', 'Personality disorder']

count=[]
total=[]
freq=[]
for p in phenos:
    c=df[df[p+'_male']==1].shape[0]+df[df[p+'_female']==1].shape[0]
    t=df[~df[p+'_male'].isnull()].shape[0]+df[~df[p+'_female'].isnull()].shape[0]
    count.append(c)
    total.append(t)
    freq.append(c*100/t)

freq_df=pd.DataFrame({'phenotype':phenos, 'count':count, 'total_assessed':total, 'frequency':freq})
print(freq_df)

# Save to file
freq_df.to_csv('Result_tables/2_phenotype_frequency.csv', index=False)