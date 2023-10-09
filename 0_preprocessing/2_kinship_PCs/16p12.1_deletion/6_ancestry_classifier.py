import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.metrics import classification_report, confusion_matrix

# Assign samples ancestry after training on the HapMap samples
df=pd.read_csv('list_files/5_16p12_hapmap_combined_output.txt', delim_whitespace=True, header=None,
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
df.loc[(df.IID.str.contains('SG')) & (df.FID.str.contains('PSU')), 'cohort']='16p12'

print(df)

hapmap=df[df.cohort=='HapMap3']
predict=df[df.cohort=='16p12']

# Predict using SVM classifier
pcs=[i for i in df.columns.to_list() if 'PC' in i]
X=hapmap[pcs]
y=hapmap.ancestry_short

print(hapmap)
print(hapmap[pcs])
print(hapmap.ancestry_short.value_counts())

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
fam=pd.read_csv('list_files/3.valid.sample', sep=' ', header=None, names=['FID', 'IID'])
fam['predicted_ancestry']=fam.IID.map(predict.predicted_ancestry.to_dict())

# Compare with self-reported ancestry
sr=pd.read_csv('analysis_files/16p12_known_ancestry.csv')
sr=sr[['IID', 'self_reported_ancestry']]
sr.index=sr.IID.to_list()
sr.fillna('.', inplace=True)

print(sr.self_reported_ancestry.unique())

# Interpret self-reported ancestry
sr['ancestry']='unknown'
sr.loc[sr.self_reported_ancestry!='.', 'ancestry']='other'
sr.loc[sr.self_reported_ancestry.str.contains('/'), 'ancestry']='multiple'
sr.loc[sr.self_reported_ancestry.isin(['Caucasian', 'White', 'European', 'British', 'Australian', 'English', 'Australian/Italian']), 'ancestry']='EUR'
sr.loc[sr.self_reported_ancestry.isin(['Chinese', 'Asian']), 'ancestry']='EAS'
sr.loc[sr.self_reported_ancestry=='African', 'ancestry']='AFR'
sr.loc[sr.self_reported_ancestry=='Indian', 'ancestry']='GIH'

# Add self-reported ancestry
sr['self_reported_ancestry']=sr.ancestry
fam=pd.merge(fam, sr[['IID','self_reported_ancestry']], on='IID', how='outer')
print(fam)

fam['check']=fam.predicted_ancestry+'.'+fam.self_reported_ancestry
print(fam.check.value_counts())

# Save
fam.to_csv('Results/6_ancestry.csv', index=False)

# Save only European samples to file
fam=fam[(fam.check=='EUR.unknown') | (fam.self_reported_ancestry=='EUR')]
fam[['FID', 'IID']].to_csv('list_files/6_european_samples.fam', sep=' ', index=None, header=None)

# Also filter out just European spouse pairs
am_samps=pd.read_csv('/data5/16p12_WGS/annotations/coding_annotations/synonymous/Analysis_files/Table_S1.csv')
am_samps=am_samps[['Sample', 'Spouse']]
am_samps=am_samps[~am_samps.Spouse.isnull()]
am_map=pd.read_csv('/data5/16p12_WGS/annotations/coding_annotations/synonymous/Analysis_files/sample_code_map.csv')
am_map.index=am_map.Sample.to_list()
am_samps['Sample_SG']=am_samps.Sample.map(am_map.Sample_SG.to_dict())
am_samps['Spouse_SG']=am_samps.Spouse.map(am_map.Sample_SG.to_dict())

am_samps=am_samps[(am_samps.Sample_SG.isin(fam.IID.to_list())) & (am_samps.Spouse_SG.isin(fam.IID.to_list()))]
spouses=fam[fam.IID.isin(am_samps.Sample_SG.to_list()+am_samps.Spouse_SG.to_list())]
print(spouses)
spouses[['FID', 'IID']].to_csv('list_files/6_european_spouses.fam', sep=' ', index=None, header=None)
