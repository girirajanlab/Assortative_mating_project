import pandas as pd

# Calculate the frequency of the phenotypes used for spouse comparisons in SPARK
phenos=['mood_anx', 'mood_bipol', 'mood_dep', 'pers_dis']

df=pd.read_csv('Analysis_files/1_family_phenotypes.csv')
df=df[df.spouse=='X']

count=[]
total=[]
freq=[]
for p in phenos:
    t=df[~df[p+'_mother'].isnull()].shape[0]+df[~df[p+'_father'].isnull()].shape[0]
    c=df[df[p+'_mother']==1].shape[0]+df[df[p+'_father']==1].shape[0]
    count.append(c)
    total.append(t)
    freq.append(c*100/t)

freq_df=pd.DataFrame({'phenotype':phenos, 'count':count, 'total_assessed':total, 'frequency':freq})
print(freq_df)

# Save to file
freq_df.to_csv('Result_tables/5_phenotype_frequency.csv', index=False)