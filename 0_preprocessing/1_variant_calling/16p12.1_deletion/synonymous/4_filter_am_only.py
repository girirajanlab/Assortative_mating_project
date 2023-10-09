import pandas as pd

# Combine synonymous and other variants
syn=pd.read_csv('tables/3_gencode.csv')
syn['Mut_type']='synonymous'
var=pd.read_csv('../nonsynonymous/Rare_Deleterious_Exonic_Variants.csv')

df=pd.concat([syn, var])
print(df)

# Remove mitochrondrial variants
df=df[df.Chrom!='chrM']

#Filter for variants private in related individuals
ped=pd.read_csv('Analysis_files/Feb_2023.fam', sep='\t')
def get_parents(samp):
	parents=ped[ped.IID==samp]['Mother'].to_list()+ped[ped.IID==samp]['Father'].to_list()
	return [i for i in parents if i !='0']
def get_kids(samp):
	kids=ped[(ped.Mother==samp) | (ped.Father==samp)]['IID'].to_list()
	return kids
def get_siblings(samp):
	pars=get_parents(samp)
	sibs=ped[((ped.Mother.isin(pars)) | (ped.Father.isin(pars))) & (ped.IID!=samp)]['IID'].to_list()
	return sibs
def related_lst(samp):
	relatives=get_parents(samp)
	kids=get_kids(samp)+get_siblings(samp)
	last_len=0
	while len(relatives)!=last_len:
		last_len=len(relatives)
		to_add=[]
		for r in relatives:
			to_add+=[i for i in get_parents(r) if i not in relatives]
			kids+=[i for i in get_kids(r) if i not in kids]
		relatives+=list(set(to_add))
	last_len=0
	kids=list(set(kids))
	while len(kids)!=last_len:
		last_len=len(kids)
		to_add=[]
		for k in kids:
			to_add+=[i for i in get_kids(k) if i not in kids]
		kids+=list(set(to_add))
	out=list(set(relatives+kids))
	out.sort()
	return out

samps=list(df.Sample.unique())
samps.sort()
rel_dict={}
for s in samps:
	rels=related_lst(s)
	rel_dict[s]=rels

df['FID']=df.Sample.map(dict(zip(ped.IID.to_list(), ped.FID.to_list())))
uniq_vars=list(df.variant_id.unique())
private_vars=[]
for v in uniq_vars:
	if len(list(df[df.variant_id==v]['FID'].unique()))>1:
		continue
	var_samps=df[df.variant_id==v]['Sample'].to_list()
	for s in var_samps:
		if len(list(set(var_samps) - set(rel_dict[s]+[s])))!=0:
			continue
	private_vars.append(v)

df=df[df.variant_id.isin(private_vars)]

print(df)
print(df.cohort_count.value_counts())

# Save to file
df.to_csv('tables/4_combined_private_in_family.csv', index=False)

# Filter for families used in AM paper
am=pd.read_csv('Analysis_files/Table_S1.csv')
map=pd.read_csv('Analysis_files/sample_code_map.csv')
map.index=map.Sample.to_list()
am['Sample_SG']=am.Sample.map(map.Sample_SG.to_dict())

df=df[df.Sample.isin(am.Sample_SG.to_list())]
print(df)

# Save to file
df.to_csv('tables/4_am_private_in_family.csv', index=False)
