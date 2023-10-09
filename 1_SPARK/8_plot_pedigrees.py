import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Wedge
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42

from matplotlib.backends.backend_pdf import PdfPages

# Make pedigrees for families with variants

def plot_shape(ax, sex, location):
	if sex=='Male':
		ax.add_artist(Rectangle((location[0], location[1]), 0.1, 0.1, facecolor='None', edgecolor='k', zorder=2))
	elif sex=='Female':
		ax.add_artist(Circle((location[0]+0.05, location[1]+0.05), 0.05, facecolor='None', edgecolor='k', zorder=2))

pheno_dict={'asd':'#F8CBAD',
			'behav_adhd':'#B4C7E7',
			'mood_anx':'#C5E0B4',
			'mood_dep':'#FFE699',
			'sleep_dx':'#CCCCFF',
			'mood_ocd':'#FFBDDE',
			'mood_bipol':'#B7FFFD'}
def plot_phenotypes(ax, sex, location, pheno_lst):
	if sex=='Male':
		# For males, draw rectangles with width=.1/(number of phenotypes)
		w=.1/len(pheno_lst)
		locx=location[0]
		for p in pheno_lst:
			ax.add_artist(Rectangle((locx, location[1]), w, 0.1, facecolor=pheno_dict[p], edgecolor=None, zorder=1))
			locx+=w
	if sex=='Female':
		# For females, draw pie slices for phenotypes
		d_theta=360/len(pheno_lst)
		theta=90
		for p in pheno_lst:
			ax.add_artist(Wedge((location[0]+0.05, location[1]+0.05), 0.05, theta, theta+d_theta, facecolor=pheno_dict[p], edgecolor=None, zorder=1))
			theta+=d_theta

# Define locations for families
m_loc=[0.6, 0.8]
f_loc=[0.4, 0.8]
# Child locations depend on the number of children
def child_locs(n_kids, y=0.6):
	total_width=(0.1*n_kids)+(.1*(n_kids-1))
	start_loc=.55-(total_width/2)
	locs=[]
	for k in range(n_kids):
		locs.append([start_loc, y])
		start_loc+=.2
	return locs

def child_lines(n_kids, y1=0.75, xmid=0.55):
	total_width=(0.1*n_kids)+(.1*(n_kids-1))-.1
	
	ax.plot([xmid-(total_width/2), xmid+(total_width/2)], [y1, y1], zorder=0, color='k', lw=2)
	start=xmid-(total_width/2)
	for k in range(n_kids):
		ax.plot([start, start], [y1, y1-.1], zorder=0, color='k', lw=2)
		start+=.2

# Load in variant, phenotype, and family information
var_df=pd.read_csv('Result_tables/7_geneset_variants_restricted.csv')
pheno_df=pd.read_csv("SPARK data/SPARK_Collection_version5/basic_medical_screening.csv")
pheno_df=pheno_df[['subject_sp_id', 'family_sf_id', 'biomother_sp_id', 'biofather_sp_id', 'sex', 'age_at_eval_months', 'asd', 'behav_adhd', 'mood_anx', 'mood_dep', 'sleep_dx', 'mood_ocd', 'mood_bipol']]
pheno_df.sort_values(by=['family_sf_id', 'age_at_eval_months'], inplace=True)

sample_vars=pd.read_csv('Analysis_files/7_var_samples.csv')

var_samps=list(var_df.Sample.unique())
var_samps.sort()

pdf=PdfPages('Figures/8_family_pedigrees.pdf')
for v in var_samps:
	# Get family members
	fam=pheno_df[pheno_df.subject_sp_id==v]['family_sf_id'].to_list()[0]
	fam_df=pheno_df[pheno_df.family_sf_id==fam]
	if v not in fam_df.subject_sp_id.to_list():
		print(fam_df, v)
	mid=fam_df[fam_df.subject_sp_id==v]['biomother_sp_id'].to_list()[0]
	fid=fam_df[fam_df.subject_sp_id==v]['biofather_sp_id'].to_list()[0]
	# Check for 3 generation
	gpar=fam_df[fam_df.subject_sp_id.isin([mid, fid])]['biomother_sp_id'].to_list()+fam_df[fam_df.subject_sp_id.isin([mid, fid])]['biofather_sp_id'].to_list()
	gpar=[i for i in gpar if i==i and i in fam_df.subject_sp_id.to_list()]
	if len(gpar)>0:
		print(fam, gpar)
	kids=[v]+[i for i in fam_df.subject_sp_id.to_list() if i not in [mid, fid, v]]

	# For each member, plot shape and phenotypes
	fig, ax = plt.subplots(figsize=(6,6), layout='constrained')
	for fm in [mid, fid]+kids:
		if fm==mid:
			loc=m_loc
		elif fm==fid:
			loc=f_loc
		else:
			loc=child_locs(len(kids))[kids.index(fm)]
		
		sex=fam_df[fam_df.subject_sp_id==fm]['sex'].to_list()[0]
		phenos=fam_df[fam_df.subject_sp_id==fm][pheno_dict.keys()].sum(axis=0)>0
		phenos=phenos[phenos].index.to_list()
		
		plot_shape(ax, sex, loc)
		if len(phenos)>0:
			plot_phenotypes(ax, sex, loc, phenos)
		
		
	
	# Add lines connecting members
	# Parent line
	ax.plot([0.45, 0.65], [0.85, 0.85], zorder=0, color='k', lw=2)
	ax.plot([0.55, 0.55], [0.75, 0.85], zorder=0, color='k', lw=2)
	child_lines(len(kids))
	
	# Annotate variants (just gene name)
	mvars='\n'.join(var_df[(var_df.Sample==v) & (var_df.inheritance=='M')]['Gene_symbol'].to_list())
	fvars='\n'.join(var_df[(var_df.Sample==v) & (var_df.inheritance=='F')]['Gene_symbol'].to_list())
	
	vids=var_df[var_df.Sample==v]['vid'].to_list()
	
	ax.text(m_loc[0], m_loc[1]-0.02, mvars, va='top')
	ax.text(f_loc[0], f_loc[1]-0.02, fvars, va='top')
	pro_loc=child_locs(len(kids))[0]
	ax.text(pro_loc[0], pro_loc[1]-0.02, fvars+'\n'+mvars, va='top')
	
	# Also annotate the kids
	if len(kids)>1:
		for i, k in enumerate(kids[1:]):
			kvars='\n'.join(sample_vars[(sample_vars.Sample==k) & (sample_vars.vid.isin(vids))]['Gene_symbol'].to_list())
			loc=child_locs(len(kids))[i+1]
			ax.text(loc[0], loc[1]-0.02, kvars, va='top')
	
	# Annotate proband
	ax.annotate('P', (pro_loc[0]-0.01, pro_loc[1]-0.01), (pro_loc[0]-0.05, pro_loc[1]-0.05), arrowprops=dict(facecolor='black', width=2))
	
	# Annotate family ID
	plt.title(fam)
	
	# Clean up
	ax.set_xlim(0, 1)
	ax.set_ylim(0, 1)
	ax.set_axis_off()
	
	pdf.savefig()
	plt.close()
	
pdf.close()