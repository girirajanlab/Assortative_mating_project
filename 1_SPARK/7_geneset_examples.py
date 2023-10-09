import pandas as pd

# Get variants inherited from parents and associated with disease in probands whose parents have phenotypes
# Get families with relevant parent and child phenotypes
phenos=pd.read_csv('../Analysis_files/1_family_phenotypes.csv')
phenos=phenos[phenos.spouse=='X']

m_cols=[i for i in phenos.columns.to_list() if '_mother' in i]
phenos['m_sum']=phenos[m_cols].sum(axis=1)
f_cols=[i for i in phenos.columns.to_list() if '_father' in i]
phenos['f_sum']=phenos[f_cols].sum(axis=1)
p_cols=[i for i in phenos.columns.to_list() if '_proband' in i]
phenos['p_sum']=phenos[p_cols].sum(axis=1)+1

phenos['type']=''
phenos.loc[(phenos.m_sum>0) & (phenos.f_sum>0), 'type']='both'
phenos.loc[(phenos.m_sum==0) & (phenos.f_sum==0), 'type']='neither'

print(phenos.type.value_counts())
phenos=phenos[(phenos.type=='both') & (phenos.p_sum>1)]

# Get probands in families
df=pd.read_csv('../SPARK_family_table.csv')
df=df[df.FID.isin(phenos.FID.to_list())]

# Subset variants
vars=pd.read_csv('Analysis_files/SPARK_inherited_genset_annotated.csv')
vars=vars[vars.Sample.isin(df.proband.to_list())]

genesets=['ASD_risk_genes_TADA_FDR0.3', 'ASD_coexpression_networks_Willsey2013', 'PSD_Genes2Cognition', 'Developmental_delay_DDD',
            'CHD8_targets_Cotney2015_Sugathan2014', 'FMRP_targets_Darnell2011', 'Geisinger_DBD_Tier', 'DD_G2P', 'SFARI_gene_score', 'SZDB_schizophrenia', 'Wang_Epilepsy']

def keep_pros(vars):
        fpros=vars[vars.inheritance=='F']['Sample'].to_list()
        mpros=vars[vars.inheritance=='M']['Sample'].to_list()
        keep_pros=list(set(fpros).intersection(set(mpros)))
        return keep_pros

# Get probands who have inherited variants in genes in multiple disease gene sets
vars=vars[vars[genesets].sum(axis=1)>3]
vars=vars[vars.Sample.isin(keep_pros(vars))]

for gs in genesets:
    print(vars[gs].value_counts())
print(vars)

print(vars.Gene_symbol.value_counts())
print(vars.Sample.value_counts())

# Save to file
vars.to_csv('Result_tables/7_geneset_variants.csv', index=False)

# Restrict to probands who have no more than than 4 total variants in this list
pro_counts=vars.Sample.value_counts()
pro_counts=pro_counts[pro_counts<=4]
vars=vars[vars.Sample.isin(pro_counts.index.to_list())]

print(vars.Gene_symbol.value_counts())
print(vars.Sample.value_counts())

# Save to file
vars.to_csv('Result_tables/7_geneset_variants_restricted.csv', index=False)

# Save just the variant IDs to file
vars['vid']=vars.Chr+'_'+vars.Pos.astype(str)+'_'+vars.Ref+'_'+vars.Alt
vars.vid.to_csv('Analysis_files/7_variant_ids.csv', index=False)