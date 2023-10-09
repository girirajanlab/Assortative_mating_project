import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

import seaborn as sns

df=pd.DataFrame()
for i,gen in enumerate(['founders', 'generation2', 'generation3']):
	g=pd.read_csv('Simulation_outputs/1/'+gen+'.csv', index_col=0)
	g['generation']=i+1
	df=pd.concat([df, g[['generation', 'disease_liability']]])
df.reset_index(inplace=True)
print(df)

ax=sns.kdeplot(data=df, x='disease_liability', hue='generation', palette='viridis_r', cut=0)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))

def zoom_in(ax_main, loc, xlim, ylim):
	axins=ax_main.inset_axes(loc)
	sns.kdeplot(data=df, x='disease_liability', hue='generation', cut=0, palette='viridis_r', legend=False, ax=axins)
	axins.set_xlim(xlim[0], xlim[1])
	axins.set_ylim(ylim[0], ylim[1])
	axins.set_xlabel('')
	axins.set_ylabel('')
	ax_main.indicate_inset_zoom(axins)

# Add zoom-ins
mid=np.mean(df.disease_liability.to_numpy())
sd=np.std(df.disease_liability.to_numpy())
mini=np.min(df.disease_liability.to_numpy())
maxi=np.max(df.disease_liability.to_numpy())

zoom_in(ax, [0.05, 0.5, 0.25, 0.25], [mini, mid-3*(sd)], [0, 0.003])
zoom_in(ax, [0.65, 0.5, 0.25, 0.25], [mid+3*(sd), (-1)*mini], [0, 0.003])

plt.savefig('Figures/5_pop_liability.pdf')
