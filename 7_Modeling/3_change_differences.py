import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Plot changes in disease liability in carriers as a function of changing parameters
df=pd.read_csv('Analysis_files/2_carrier_lineages.csv')

df=df[df.group.isin(['carrier', 'child_carrier', 'grandchild_carrier_carrier_lineage'])]

# Normalize the liability by the run
for r in list(df.run.unique()):
	df.loc[df.run==r, 'disease_liability']=df[df.run==r]['disease_liability']-df[df.run==r]['disease_liability'].mean()

# For each run, save the parameters and average and SD disease liability in each group to file
out=[]
for r in list(df.run.unique()):
	params=[r, df[df.run==r]['h2'].to_list()[0], df[df.run==r]['corr'].to_list()[0], df[df.run==r]['n_loci'].to_list()[0], df[df.run==r]['common_prop'].to_list()[0], df[df.run==r]['effect_ratio'].to_list()[0], df[df.run==r]['pop_size'].to_list()[0]]
	for g in 'carrier', 'child_carrier', 'grandchild_carrier_carrier_lineage':
		u=df[(df.run==r) & (df.group==g)]['disease_liability'].mean()
		sd=df[(df.run==r) & (df.group==g)]['disease_liability'].std()
		params.append(('%.6f (%.6f)' % (u, sd)))
	out.append(params)
outdf=pd.DataFrame(out, columns=['run', 'h2', 'corr', 'n_loci', 'common_prop', 'effect_ratio', 'pop_size', 'g1', 'g2', 'g3'])
outdf.to_csv('Result_tables/3_paramter_differences.csv', index=False)
