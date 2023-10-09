import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.metrics import classification_report, confusion_matrix

predicted=pd.DataFrame()
for i in range(1, 1001):
	print('Run', str(i))
	# Assign samples ancestry after training on the HapMap samples
	df=pd.read_csv('Results/faststructure_output/'+str(i)+'_output.txt', delim_whitespace=True, header=None,
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
	df.loc[df.FID.str.contains('SF'), 'cohort']='SPARK'

	hapmap=df[df.cohort=='HapMap3'].copy()
	predict=df[df.cohort=='SPARK'].copy()

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

	predicted=pd.concat([predicted, predict])

predicted.index=predicted.IID.to_list()

# Add predicted ancestry to FAM file
fam=pd.read_csv('../plink_files/3_ld_het_filter/SPARK_QC_final.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
fam['predicted_ancestry']=fam.IID.map(predicted.predicted_ancestry.to_dict())

# Compare with self-reported ancestry
sr=pd.read_csv('../analysis_files/individuals_registration.csv')
keep_cols=['subject_sp_id', 'family_sf_id']+[i for i in sr.columns.to_list() if 'race' in i]
sr=sr[keep_cols]
sr.index=sr.subject_sp_id.to_list()

# Interpret self-reported ancestry
sr['ancestry']='unknown'
sr.loc[sr.race_asian==1, 'ancestry']='EAS'
sr.loc[sr.race_african_amer==1, 'ancestry']='AFR'
sr.loc[sr.race_native_amer==1, 'ancestry']='NAT'
sr.loc[sr.race_native_hawaiian==1, 'ancestry']='NAH'
sr.loc[sr.race_white==1, 'ancestry']='EUR'
sr.loc[sr.race_other==1, 'ancestry']='other'
sr.loc[sr.race_more_than_one_calc==1, 'ancestry']='multiple'

# Add self-reported ancestry
sr['IID']=sr.subject_sp_id
sr['self_reported_ancestry']=sr.ancestry
fam=pd.merge(fam, sr[['IID', 'self_reported_ancestry']], on='IID', how='outer')
print(fam)

print(fam.self_reported_ancestry.value_counts())
print(fam.predicted_ancestry.value_counts())

fam['check']=fam.predicted_ancestry+'.'+fam.self_reported_ancestry
print(fam.check.value_counts())

# Save ancestry to file
fam.to_csv('Results/2_ancestry.csv', index=False)

fam=fam[(fam.check=='EUR.unknown') | (fam.self_reported_ancestry=='EUR')]
fam.to_csv('Results/2_european_samples.csv', index=False)

# Save only european spouse pairs to file
pairs=pd.read_csv('../analysis_files/individuals_registration.csv')
pairs=pairs[(~pairs.biomother_sp_id.isnull()) & (~pairs.biofather_sp_id.isnull())]
pairs=pairs[(pairs.biomother_sp_id.isin(fam.IID.to_list())) & (pairs.biofather_sp_id.isin(fam.IID.to_list()))]

fam=fam[fam.IID.isin(pairs.biomother_sp_id.to_list()+pairs.biofather_sp_id.to_list())]

fam[['FID', 'IID']].to_csv('Results/2_european_spouses.txt', sep=' ', index=False)
