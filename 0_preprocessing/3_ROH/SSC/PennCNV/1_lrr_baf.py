import pandas as pd
import subprocess

# Isolate the LRR and BAF values from the final reports for probands used in ROH analysis
df=pd.read_csv('../Analysis_files/1_roh_samples.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])

# Get sample filenames
files=pd.read_csv('Analysis_files/file_locations.list', header=None, names=['filename'])
files['IID']=files.filename.str.split('/', expand=True)[7].str.split('_', expand=True)[0]
files['seq_id']=files.filename.str.split('/', expand=True)[7].str.split('_', expand=True)[1]+'_'+files.filename.str.split('/', expand=True)[7].str.split('_', expand=True)[2]
files['batch']=files.filename.str.split('/', expand=True)[6]

# Add in sexes from fam files
fam_files=['1Mv1/SSC_1Mv1_binary.fam', '1Mv3/SSC_1Mv3_binary.fam', 'Omni2.5/SSC_Omni2.5_binary.fam']
ff_df=pd.DataFrame()
for ff in fam_files:
	fdf=pd.read_csv(ff, sep=' ', header=None, names=['FID', 'seq_id', 'Father', 'Mother', 'Sex', 'Phenotype'])
	ff_df=pd.concat([ff_df, fdf])
files=pd.merge(files, ff_df[['seq_id', 'Sex']], on='seq_id', how='left')
files['Sex']=files.Sex.map({1:'male', 2:'female', 0:'unknown'})

df=pd.merge(df[['IID']], files, on='IID', how='left')

print(df)

for idx, row in df.iterrows():
	filename=row['filename']
	samp=row['IID']

	file=pd.read_csv(filename, skiprows=10, sep='\t')[['SNP Name', 'Log R Ratio', 'B Allele Freq']]
	file.columns=['Name', samp+'.Log R Ratio', samp+'.B Allele Freq']

	file.to_csv("input_files/"+row['batch']+"/"+samp+".txt", index=False, sep='\t')

df['input_filename']="input_files/"+df.batch+"/"+df.IID+".txt"

# Make lists of input files and a sexfile for each batch
batches=list(df.batch.unique())
for b in batches:
	subdf=df[df.batch==b]
	subdf[['input_filename']].to_csv('Analysis_files/'+b+'_files.list', index=False, header=False)
	subdf[['input_filename', 'Sex']].to_csv('Analysis_files/'+b+'_sexfile.list', index=False, header=False, sep='\t')

# Also make a SNP location file for each batch
# Take one files from each batch and extract the SNP locations
for b in batches:
	filename=df[df.batch==b]['filename'].to_list()[0]
	file=pd.read_csv(filename, skiprows=10, sep='\t')[['SNP Name', 'Chr', 'Position']]
	file.columns=['Name', 'Chr', 'Pos']
	file.to_csv('Analysis_files/'+b+'_snp_pos.txt', index=False, sep='\t')

# Save entire dataframe
df.to_csv('Analysis_files/1_sample_info.csv', index=False)
