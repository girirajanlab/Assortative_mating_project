import pandas as pd
import numpy as np
import random
import scipy.stats as stats
import sys

batch_id=sys.argv[1]

h2=float(sys.argv[2])
env_lia=np.sqrt(1-h2)
corr=float(sys.argv[3])
n_loci=int(sys.argv[4])
n_common=int(float(sys.argv[5])*n_loci)
n_rare=n_loci-n_common
effect_factor=float(sys.argv[6])
pop_n=int(sys.argv[7])

# Set random seed for reproducability
random.seed(205)

# Generate loci with MAF and effect size
common_loci=[str(i) for i in range(1, n_common+1)]
rare_loci=[str(i) for i in range(n_common+1, n_common+n_rare+1)]
loci=common_loci+rare_loci

# Vary the effects of rare and common variants so each class contributes equally to disease variance
eff_common=np.sqrt(h2/(n_common*(1+1/effect_factor)))
eff_rare=eff_common*np.sqrt(n_common/(effect_factor*n_rare))

maf=[random.uniform(0.05, 0.5) for i in range(n_common)]+[random.uniform(0, 0.001) for i in range(n_rare)]
effect=[random.normalvariate(0, eff_common) for i in range(n_common)]+[random.normalvariate(0, eff_rare) for i in range(n_rare)]

# Save variant information to file for later use
var_df=pd.DataFrame({'ID':loci, 'MAF':maf, 'effect':effect})
top_decile=var_df[var_df.MAF<=0.001]['effect'].quantile(q=0.9)
high_effect_rare=var_df[(var_df.MAF<=0.001) & (var_df.effect>=top_decile)]['ID'].to_list()
var_df.to_csv('Simulation_outputs/'+batch_id+'/variant_information.csv', index=False)
del var_df

# Calculate the liability in a person from the additive effects of alleles at all loci
def cal_liability(df_in):
	df=df_in[loci].copy()

	# Standardize counts at each locus: (allele_count-2(MAF))/sqrt(2MAF*(1-MAF))
	df=df.subtract(2*pd.Series(maf, index=loci), axis='columns')
	df=df.div(np.sqrt((2*pd.Series(maf, index=loci))*(1-pd.Series(maf, index=loci))), axis='columns')
	df=df.mul(effect, axis='columns')

	# Check the liability contributions from rare and common variants
	df['common_liability']=df[common_loci].sum(axis=1)
	print('common', df.common_liability.mean(), df.common_liability.std(), df.common_liability.var())
	df['rare_liability']=df[rare_loci].sum(axis=1)
	print('rare', df.rare_liability.mean(), df.rare_liability.std(), df.rare_liability.var())

	df['genetic_liability']=df[loci].sum(axis=1)
	df['environmental_liability']=[random.normalvariate(0, env_lia)for i in range(df.shape[0])]
	print('environment', df.environmental_liability.mean(), df.environmental_liability.std(), df.environmental_liability.var())
	df['disease_liability']=df.genetic_liability+df.environmental_liability
	print('total',  df.disease_liability.mean(), df.disease_liability.std(), df.disease_liability.var())

	return df[['common_liability', 'rare_liability', 'genetic_liability', 'environmental_liability', 'disease_liability']]

# Assign partners based on maintaining a specific correlation in liability
def assign_partners(df_in):
	df=df_in[['disease_liability']].copy()

	# Split into males and females (equal numbers of each)
	mdf=df.iloc[0:int(df.shape[0]/2)].copy()
	fdf=df.iloc[int(df.shape[0]/2):].copy()

	# Assign random liabilities corresponding to a correlation of corr with the male liabilities
	mdf['start_val']=fdf.disease_liability.to_numpy()
	mdf=mdf.subtract(mdf.mean(axis=0))
	id_mat=np.identity(mdf.shape[0])
	q, r=np.linalg.qr(mdf['disease_liability'].to_numpy()[:, None])
	p_lst=[]
	for idx in range(q.shape[0]):
		p_lst.append(q*q[idx])
	p=np.reshape(np.array(p_lst), (mdf.shape[0], mdf.shape[0]))
	mdf['orth_val']=np.matmul((np.subtract(id_mat, p)), mdf['start_val'].to_numpy())
	sq_df=mdf[['disease_liability', 'orth_val']].pow(2)
	id_mat2=np.diag(1/np.sqrt(sq_df.sum(axis=0).to_numpy()))
	scaled_col=pd.DataFrame(np.matmul(mdf[['disease_liability', 'orth_val']].to_numpy(), id_mat2), columns=['scaled_gl', 'scaled_rand'])
	mdf['final_vals']=scaled_col.scaled_rand+scaled_col.scaled_gl*(1/np.tan(np.arccos(corr)))

	# Order frame by these values
	mdf.sort_values(by='final_vals', inplace=True)
	# Order female values and merge
	fdf.sort_values(by='disease_liability', inplace=True)
	print('Liability correlation is', stats.pearsonr(mdf.disease_liability.to_numpy(), fdf.disease_liability.to_numpy())[0])
	pairs=dict(zip(mdf.index.to_list()+fdf.index.to_list(), fdf.index.to_list()+mdf.index.to_list()))
	df['partner']=df.index.map(pairs)
	return df.partner

# Create offspring from pairs
# Each parent has a probability of passing on a variant corresponding to their genotype (homozygous variant:100%, heterozygous:50%, homozygous reference:0%)
def create_kids(df_in):
	df=df_in[loci].copy()
	# Get a dictionary of partner pairs
	partners=dict(zip(df.index.to_list(), ['.'.join(sorted(list(i))) for i in list(zip(df.index.astype(str).to_list(), df_in.partner.astype(int).astype(str).to_list()))]))
	# Map alleles to probability
	df.replace({1:0.5, 2:1}, inplace=True)
	# Each pair will have 2 kids - decide on alleles for each child independently
	kid_df=pd.DataFrame()
	for i in range(2):
		# To decide if allele is passed on, compare probability to random numbers
		comp_df=pd.DataFrame(np.random.uniform(low=0, high=1, size=df.shape), columns=df.columns, index=df.index)
		booldf=df.gt(comp_df).astype(int)
		booldf['partner_id']=df_in.index.map(partners)

		child_df=booldf.groupby(['partner_id']).sum()
		kid_df=pd.concat([kid_df, child_df])
	kid_df['parents']=kid_df.index.to_list()
	kid_df.sort_values(by='parents', inplace=True)
	kid_df.reset_index(drop=True, inplace=True)
	return kid_df

def create_founders(n=100000):
	# Define parental alleles based on MAF
	allele_probs=pd.DataFrame(np.repeat(np.reshape(np.array(maf), (1, n_common+n_rare)), n, axis=0), columns=loci)
	a1=allele_probs.gt(pd.DataFrame(np.random.uniform(low=0, high=1, size=allele_probs.shape), columns=allele_probs.columns, index=allele_probs.index)).astype(int)
	a2=allele_probs.gt(pd.DataFrame(np.random.uniform(low=0, high=1, size=allele_probs.shape), columns=allele_probs.columns, index=allele_probs.index)).astype(int)

	df=a1.add(a2)

	# Define liability
	df=pd.merge(df, cal_liability(df), left_index=True, right_index=True)

	# Assign partners
	df['partner']=assign_partners(df)

	# Annotate samples for presence of a rare variant with large effect size
	df['he_rare_carrier']=(df[high_effect_rare].sum(axis=1)>0).astype(int)

	# Save
	df[['common_liability', 'rare_liability', 'genetic_liability', 'environmental_liability', 'disease_liability', 'partner', 'he_rare_carrier']].to_csv('Simulation_outputs/'+batch_id+'/founders.csv')

	# Get kids
	return create_kids(df)

def generation(df, gen_num, genmax=3):
	# Define liability
	df=pd.merge(df, cal_liability(df), left_index=True, right_index=True)

	# Assign partners
	if gen_num>=genmax:
		df['partner']=0
	else:
		df['partner']=assign_partners(df)

	# Annotate samples for presence of a rare variant with large effect size
	df['he_rare_carrier']=(df[high_effect_rare].sum(axis=1)>0).astype(int)

	# Save
	df[['common_liability', 'rare_liability', 'genetic_liability', 'environmental_liability', 'disease_liability', 'partner', 'he_rare_carrier', 'parents']].to_csv('Simulation_outputs/'+batch_id+'/generation'+str(gen_num)+'.csv')

	# Get kids
	if gen_num<genmax:
		return create_kids(df)

olddf=create_founders(n=pop_n)
for i in range(2):
	newdf=generation(olddf, i+2)
	olddf=newdf
