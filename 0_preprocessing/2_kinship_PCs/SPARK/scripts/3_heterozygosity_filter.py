import pandas as pd
import sys
import numpy as np

if len(sys.argv)!=4:
	print('You must provide 3 file names to this function - QC file, valid sample file, and invalid sample file')
	exit()

qc=sys.argv[1]
valid=sys.argv[2]
invalid=sys.argv[3]

df=pd.read_csv(qc, delim_whitespace=True)

mean_F=np.mean(df.F.to_numpy())
sd_F=np.std(df.F.to_numpy())

upper=mean_F+(2*sd_F)
lower=mean_F-(2*sd_F)

df['keep']=False
df.loc[(df.F>=lower) & (df.F<=upper), 'keep']=True

print(df.keep.value_counts())

df[df.keep][['FID', 'IID']].to_csv(valid, sep=' ', index=False)
df[~df.keep][['FID', 'IID']].to_csv(invalid, sep=' ', index=False)
