import pandas as pd
import csv
import subprocess
import random

# CNVs only need to be called on proband samples used for ROH analysis
# Take the sample list from the ROH analysis
samps=pd.read_csv('../Result_files/3.1_1kb_window.hom.indiv', sep='\s+')
fam=pd.read_csv('../Analysis_files/1_roh_samples.fam', sep=' ', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
df=pd.merge(samps[['FID', 'IID']], fam[['IID', 'Sex']], on='IID', how='left')
df.Sex=df.Sex.map({1:'male', 2:'female'})

# Get the LRR and BAF filename for each sample
files=pd.read_csv('Analysis_files/file_locations.list', header=None, names=['filename'])
files['IID']='SP0'+files.filename.str.split('SP0', expand=True)[1].str.split('_', expand=True)[0].str.split('.', expand=True)[0]

# If a single sample appears more than once, take only the most recent data
repeats=files.IID.value_counts()
repeats=repeats[repeats>1]

files['keep']=False
files.loc[~files.IID.isin(repeats.index.to_list()), 'keep']=True
files.loc[files.filename.str.contains('iWES'), 'keep']=True
files=files[files.keep]

# Get batch of data
files['batch']=''
files.loc[files.filename.str.contains('iWES'), 'batch']='iWES'
files.loc[files.filename.str.contains('Data1'), 'batch']='Data1'
files.loc[files.filename.str.contains('Pilot'), 'batch']='Pilot'

df=pd.merge(df, files[['IID', 'filename', 'batch']], on='IID', how='left')
print(df)
print(df.batch.value_counts())

# All samples are available in iWES except one - skip that sample and only assess others
print(df)
df=df[df.batch=='iWES']
df.reset_index(drop=True, inplace=True)
print(df)

# Add gunzip command to filenames
df['command']='`echo \"Name\\t'+df.IID+'.Log R Ratio\\t'+df.IID+'.B Allele Freq\" ; (gunzip -c '+df.filename+' | tail -n +2)`'
# Save commands to file
df[df.batch=='iWES'][['command']].to_csv('Analysis_files/iwes_command.list', index=False, header=False, doublequote=False, quoting=csv.QUOTE_NONE)

# To make PFB files, we need some unzipped files
# Select 2000 random files form iWES
random.seed(205)
idx=random.choices(df.index.to_list(), k=2000)
subdf=df.iloc[idx]
files=subdf.filename.to_list()
inputs=[]
for f in files:
	samp=f.split('/')[-1].split('_')[0].split('.')[0]
	outname='PFB_inputs/iWES/'+f.split('/')[-1].split('.gz')[0]
	with open(outname, "w") as outfile:
		cmd=subprocess.run(["gunzip", "-c", f], capture_output=True, text=True)
		new=cmd.stdout.replace("Name\tLog_R_Ratio\tB_Allele_Freq", "Name\t"+samp+".Log R Ratio\t"+samp+".B Allele Freq\n")
		outfile.write(new)
	inputs.append(outname)
with open("Analysis_files/iWES_pfb_input.list", "w") as lstfile:
	lstfile.write("\n".join(inputs))

# Make sex file
df[['command', 'Sex']].to_csv('Analysis_files/iwes_sexfile.list', index=False, header=False, sep='\t', doublequote=False, quoting=csv.QUOTE_NONE)

# Save full dataframe for later use
df.to_csv('Analysis_files/1_samples.csv', index=False)

# While here, also extract the SNP positions for each array
iwes=pd.read_csv('download/resources/GSA-24v2-0_A2.csv', skiprows=7)
iwes=iwes[~iwes.MapInfo.isnull()]
iwes['Pos']=iwes.MapInfo.astype(int)
iwes[['Name', 'Chr', 'Pos']].to_csv('Analysis_files/iWES_snp_pos.txt', sep='\t', index=False)
