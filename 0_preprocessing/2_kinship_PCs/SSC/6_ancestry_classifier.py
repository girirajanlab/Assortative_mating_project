import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.metrics import classification_report, confusion_matrix

# Assign samples ancestry after training on the HapMap samples
df=pd.read_csv('list_files/5_SSC_hapmap_combined_output.txt', delim_whitespace=True, header=None,
		names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6'])

# Add HapMap3 ancestry
# Ancestry from: https://www.broadinstitute.org/files/shared/mpg/hapmap3/relationships_w_pops_051208.txt
hapmap=pd.read_csv('hapmap3/relationships_w_pops_051208.txt', sep='\t')
hapmap.index=hapmap.IID.to_list()
df['known_ancestry']=df.IID.map(hapmap.population.to_dict())

# Combine ancestries into larger ancestry groups
anc_dict={'ASW':'AFR', 'LWK':'AFR', 'MKK':'AFR', 'YRI':'AFR',
		'CEU':'EUR', 'TSI':'EUR',
		'CHB':'EAS', 'CHD':'EAS', 'JPT':'EAS',
		'GIH':'GIH', 'MEX':'MEX'}
df['ancestry_short']=df.known_ancestry.map(anc_dict)

# Split samples based on cohort
df['cohort']='HapMap3'
# Get SSC sample IDs
ssc_samps=pd.read_csv('SSC_microarray_samples.update.txt', sep='\t', header=None, names=['FID', 'SEQ_ID', 'FID2', 'IID'])
df.loc[df.IID.isin(ssc_samps.SEQ_ID.to_list()), 'cohort']='SSC'

hapmap=df[df.cohort=='HapMap3']
predict=df[df.cohort=='SSC']

# Predict using SVM classifier
pcs=[i for i in df.columns.to_list() if 'PC' in i]
X=hapmap[pcs]
y=hapmap.ancestry_short

# Split into train and test
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.10)

# Standardize features:
scaler = StandardScaler()
scaler.fit(X_train)

X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

# Fit data:
classifier = svm.SVC(decision_function_shape='ovr')
classifier.fit(X_train, y_train)

# Predict y data with classifier:
y_predict = classifier.predict(X_test)

print(confusion_matrix(y_test, y_predict))
print(classification_report(y_test, y_predict))

new_prediction=classifier.predict(scaler.transform(predict[pcs]))
predict['predicted_ancestry']=new_prediction

print(predict)
print(predict.predicted_ancestry.value_counts())
predict.index=predict.IID.to_list()

# Add predicted ancestry to FAM file
fam=pd.read_csv('plink_files/3_ld_het_filter/SSC_QC_final.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
fam['predicted_ancestry']=fam.IID.map(predict.predicted_ancestry.to_dict())

# Compare with self-reported ancestry
sr=pd.read_csv('analysis_files/ssc_background_hx.csv')
keep_cols=['individual']+[i for i in sr.columns.to_list() if 'race' in i]
sr=sr[keep_cols]
sr['family']=sr.individual.str.split('.', expand=True)[0]
sr['mother']=sr.family+'.mo'
sr['father']=sr.family+'.fa'

race_cols=[i for i in sr.columns.to_list() if 'race' in i and 'father' not in i and 'mother' not in i]
mdf=sr.copy()
fdf=sr.copy()

for rc in race_cols:
	mdf[rc]=mdf['mother_'+rc]
	fdf[rc]=fdf['father_'+rc]
mdf['individual']=mdf.mother
fdf['individual']=fdf.father

sr_indiv=pd.concat([sr[['individual']+race_cols], mdf[['individual']+race_cols], fdf[['individual']+race_cols]])

# Interpret self-reported ancestry
sr_indiv['ancestry']='unknown'
sr_indiv.loc[sr_indiv.race_asian==1, 'ancestry']='EAS'
sr_indiv.loc[sr_indiv.race_african_amer==1, 'ancestry']='AFR'
sr_indiv.loc[sr_indiv.race_native_amer==1, 'ancestry']='NAT'
sr_indiv.loc[sr_indiv.race_native_hawaiian==1, 'ancestry']='NAH'
sr_indiv.loc[sr_indiv.race_white==1, 'ancestry']='EUR'
sr_indiv.loc[sr_indiv.race_other==1, 'ancestry']='other'
sr_indiv.loc[sr_indiv.race_more_than_one==1, 'ancestry']='multiple'

# Add self-reported ancestry
fam['ID']=fam.IID.map(dict(zip(ssc_samps.SEQ_ID.to_list(), ssc_samps.IID.to_list())))
sr_indiv['ID']=sr_indiv.individual
sr_indiv['self_reported_ancestry']=sr_indiv.ancestry
fam=pd.merge(fam, sr_indiv[['ID', 'self_reported_ancestry']], on='ID', how='outer')
print(fam)

fam['check']=fam.predicted_ancestry+'.'+fam.self_reported_ancestry
print(fam.check.value_counts())

# Save to file
fam.to_csv('Results/6_ancestry.csv', index=False)

# European samples
fam=fam[(fam.check=='EUR.unknown') | (fam.self_reported_ancestry=='EUR')]
fam[['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']].to_csv('list_files/6_european_samples.fam', index=False, sep=' ')

# European spouse pairs
fam['ID2']=fam.IID.map(dict(zip(ssc_samps.SEQ_ID.to_list(), ssc_samps.IID.to_list())))
fam=fam[(fam.ID2.str.contains('mo')) | (fam.ID2.str.contains('fa'))]
spouses=fam.FID.value_counts()
spouses=spouses[spouses==2]
fam=fam[fam.FID.isin(spouses.index.to_list())]
fam[['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']].to_csv('list_files/6_european_spouses.fam', index=False, sep=' ')
