import pandas as pd

# Gather the disease liability for rare variant carriers across generations

def gen1(df_in):
	df=df_in.copy()
	carriers=df[df.he_rare_carrier==1].index.to_list()
	partners=df[df.he_rare_carrier==1].partner.to_list()

	df=df[df.index.isin(carriers+partners)]

	# If both carrier and partner have a high effect variant, assign the one with the lower index as the carrier
	overlap=list(set(carriers).intersection(set(partners)))
	df['group']='carrier'
	df.loc[df.index.isin(partners), 'group']='partner'
	df.loc[(df.index.isin(overlap)) & (df.index<df.partner), 'group']='carrier'

	df['generation']=1
	return df

def gen2(df_in, gen1):
	df=df_in.copy()
	df['par1']=df.parents.astype(str).str.split('.', expand=True)[0].astype(int)
	df['par2']=df.parents.astype(str).str.split('.', expand=True)[1].astype(int)

	kids=df[(df.par1.isin(gen1)) | (df.par2.isin(gen1))].index.to_list()
	partners=df[df.partner.isin(kids)].index.to_list()

	df=df[(df.index.isin(kids)) | (df.index.isin(partners))]

	df['group']='child'
	df.loc[~df.index.isin(kids), 'group']='partner'

	# If both partners are children of high-risk variant carrier, designate one who did not inherit the variant as the partner
	carriers=df[df.he_rare_carrier==1].index.to_list()
	df.loc[(df.partner.isin(carriers)) & (~df.index.isin(carriers)), 'group']='partner'

	# If both partners are still labeled as children, assign the one with the lower index as the child
	kids=df[df.group=='child'].index.to_list()
	df.loc[(df.index.isin(kids)) & (df.partner.isin(kids)) & (df.index<df.partner), 'group']='partner'

	# Refine labels to include carrier status
	df.loc[(df.group=='child') & (df.he_rare_carrier==1), 'group']='child_carrier'
	df.loc[(df.group=='child') & (df.he_rare_carrier==0), 'group']='child_noncarrier'

	carrier_kids=df[df.group=='child_carrier'].index.to_list()
	df.loc[(df.group=='partner') & (df.partner.isin(carrier_kids)), 'group']='child_carrier_partner'
	df.loc[(df.group=='partner') & (~df.partner.isin(carrier_kids)), 'group']='child_noncarrier_partner'

	df['generation']=2
	return df

def gen3(df_in, child_c, child_nc):
	df=df_in.copy()
	df['par1']=df.parents.astype(str).str.split('.', expand=True)[0].astype(int)
	df['par2']=df.parents.astype(str).str.split('.', expand=True)[1].astype(int)

	kids=df[(df.par1.isin(child_c+child_nc)) | (df.par2.isin(child_c+child_nc))].index.to_list()

	df=df[df.index.isin(kids)]

	# Refine labels
	df['group']='grandchild_noncarrier_noncarrier_lineage'
	df.loc[df.he_rare_carrier==1, 'group']='grandchild_carrier_carrier_lineage'
	df.loc[((df.par1.isin(child_c)) | (df.par2.isin(child_c))) & (df.he_rare_carrier==0), 'group']='grandchild_noncarrier_carrier_lineage'
	# Shouldn't exist, but can check
	df.loc[(df.group=='grandchild_noncarrier_noncarrier_lineage') & (df.he_rare_carrier==1), 'group']='grandchild_carrier_noncarrier_lineage'

	df['generation']=3
	return df

comp_df=pd.DataFrame()
for i in range(1, 29):
	print(i)
	g1=pd.read_csv('Simulation_outputs/'+str(i)+'/founders.csv', index_col=0)
	df=gen1(g1)
	del g1

	g2=pd.read_csv('Simulation_outputs/'+str(i)+'/generation2.csv', index_col=0)
	df=pd.concat([df, gen2(g2, df.index.to_list())])
	del g2

	g3=pd.read_csv('Simulation_outputs/'+str(i)+'/generation3.csv', index_col=0)
	cc=df[(df.generation==2) & (df.group=='child_carrier')].index.to_list()
	cnc=df[(df.generation==2) & (df.group=='child_noncarrier')].index.to_list()
	df=pd.concat([df, gen3(g3, cc, cnc)])

	df['run']=i

	comp_df=pd.concat([comp_df, df])

# Annotate runs with parameter information
para=pd.read_csv('0_model_parameters.txt', sep=' ', header=None, names=['h2', 'corr', 'n_loci', 'common_prop', 'effect_ratio', 'pop_size'])
para['run']=para.index.astype(int)+1

outdf=pd.merge(comp_df, para, on='run')

# Save
outdf.to_csv('Analysis_files/2_carrier_lineages.csv')
