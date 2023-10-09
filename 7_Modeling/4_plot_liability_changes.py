import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.patches import Rectangle, Circle
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42

import seaborn as sns

df=pd.read_csv('Analysis_files/2_carrier_lineages.csv', index_col=0)
df=df[df.run==1]

print(df)

# For each type of person, get the average genetic liability, variant burden, and proportion with high risk vars
groups=['carrier', 'partner',
	'child_carrier', 'child_carrier_partner', 'child_noncarrier', 'child_noncarrier_partner',
	'grandchild_carrier_carrier_lineage', 'grandchild_noncarrier_carrier_lineage', 'grandchild_noncarrier_noncarrier_lineage', 'grandchild_carrier_noncarrier_lineage']

lst=[]
for g in groups:
	subdf=df[df.group==g]
	if subdf.shape[0]==0:
		continue
	ave_lia=np.mean(subdf.disease_liability.to_numpy())
	prop_high_risk=subdf[subdf.he_rare_carrier==1].shape[0]/subdf.shape[0]
	n=subdf.shape[0]
	lst.append([g, ave_lia, prop_high_risk, n])

group_df=pd.DataFrame(lst, columns=['group', 'average_liability', 'proportion_carriers', 'num_samples'])

cmap=plt.cm.Blues
norm=matplotlib.colors.Normalize(vmin=min(group_df.average_liability.to_list()), vmax=max(group_df.average_liability.to_list()))

def plot_shape(ax, shape, liability, carrier, location):
	col=cmap(norm(liability))
	if shape=='square':
		ax.add_artist(Rectangle((location[0], location[1]), 0.1, 0.1, facecolor=col, edgecolor='k'))
	elif shape=='circle':
		ax.add_artist(Circle((location[0]+0.05, location[1]+0.05), 0.05, facecolor=col, edgecolor='k'))

	# Annotate a circle for the proportion of high risk variant carriers
	# The size of the circle  indicates the proportion
	ax.add_artist(Circle((location[0]+0.05, location[1]+0.05), 0.025*(np.sqrt(carrier)), facecolor='white', edgecolor='k'))

pdf=PdfPages('Figures/4_lineage_pedigree.pdf')
fig, axs = plt.subplots(ncols=2, figsize=(13,6))
ax=axs[0]
patch_dict={'carrier':['square', [0.6, 0.8]],
			'partner':['circle', [0.4, 0.8]],
			'child_carrier':['circle', [0.6, 0.6]],
			'child_carrier_partner':['square', [0.8, 0.6]],
			'child_noncarrier':['square', [0.4, 0.6]],
			'child_noncarrier_partner':['circle', [0.2, 0.6]],
			'grandchild_noncarrier_noncarrier_lineage':['circle', [0.2, 0.4]],
			'grandchild_noncarrier_carrier_lineage':['circle', [0.6, 0.4]],
			'grandchild_carrier_carrier_lineage':['square', [0.8, 0.4]]}
for key in patch_dict.keys():
	plot_shape(ax, patch_dict[key][0], group_df[group_df.group==key]['average_liability'].to_list()[0], group_df[group_df.group==key]['proportion_carriers'].to_list()[0], patch_dict[key][1])
plot_shape(ax, 'square', group_df[group_df.group=='grandchild_noncarrier_noncarrier_lineage']['average_liability'].to_list()[0],
				group_df[group_df.group=='grandchild_noncarrier_noncarrier_lineage']['proportion_carriers'].to_list()[0], [0.4, 0.4])

# Add lines connecting the pairs
ax.plot([0.45, 0.65], [0.85, 0.85], zorder=0, color='k', lw=2)
ax.plot([0.65, 0.85], [0.65, 0.65], zorder=0, color='k', lw=2)
ax.plot([0.25, 0.45], [0.65, 0.65], zorder=0, color='k', lw=2)
# Add lines connecting pairs to children
ax.plot([0.55, 0.55], [0.75, 0.85], zorder=0, color='k', lw=2)
ax.plot([0.75, 0.75], [0.55, 0.65], zorder=0, color='k', lw=2)
ax.plot([0.35, 0.35], [0.55, 0.65], zorder=0, color='k', lw=2)

# Weight the lines to specific children by the number of children in the category
def child_l(x1, x2, y1, y2, n, gen):
	width=2
	ax.plot([x2, x1], [y1, y1], zorder=0, color='k', lw=width)
	ax.plot([x2, x2], [y2, y1], zorder=0, color='k', lw=width)

child_l(0.55, 0.45, 0.75, 0.65, group_df[group_df.group=='child_noncarrier']['num_samples'].to_list()[0], 2)
child_l(0.55, 0.65, 0.75, 0.65, group_df[group_df.group=='child_carrier']['num_samples'].to_list()[0], 2)
child_l(0.35, 0.25, 0.55, 0.45, group_df[group_df.group=='grandchild_noncarrier_noncarrier_lineage']['num_samples'].to_list()[0]/2, 3)
child_l(0.35, 0.45, 0.55, 0.45, group_df[group_df.group=='grandchild_noncarrier_noncarrier_lineage']['num_samples'].to_list()[0]/2, 3)
child_l(0.75, 0.65, 0.55, 0.45, group_df[group_df.group=='grandchild_noncarrier_carrier_lineage']['num_samples'].to_list()[0], 3)
child_l(0.75, 0.85, 0.55, 0.45, group_df[group_df.group=='grandchild_carrier_carrier_lineage']['num_samples'].to_list()[0], 3)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()

fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axs[1], orientation='horizontal', label='disease liability')
pdf.savefig()
plt.close()

# Add a plot showing the range of liabilities for carriers in the carrier lineage
import seaborn as sns

import scipy.stats as stats

rel_groups=['carrier', 'child_carrier', 'grandchild_carrier_carrier_lineage']
sns.boxplot(data=df[df['group'].isin(rel_groups)], x='disease_liability', y='group', color='white', fliersize=0)
plt.xlim(-2.5, 4)
pdf.savefig()
plt.close()

pdf.close()

group_df.to_csv('Result_tables/4_lineage_proportions.csv', index=False)


# T-tests between carriers in each generation
stat_lst=[]
for g1 in rel_groups:
	for g2 in rel_groups:
		if rel_groups.index(g1)>=rel_groups.index(g2):
			continue
		res=stats.ttest_ind(df[df.group==g1].disease_liability.to_list(), df[df.group==g2].disease_liability.to_list(), alternative='less')
		stat_lst.append([g1, g2, df[df.group==g1].shape[0], df[df.group==g2].shape[0], res.statistic, res.pvalue])
statdf=pd.DataFrame(stat_lst, columns=['group1', 'group2', 'n_group1', 'n_group2', 'statistic', 'p'])
statdf.to_csv('Result_tables/4_ttests.csv', index=False)
