import pandas as pd


# Get spouse pair IDs
howe=pd.read_csv('dataset.txt')

# Link the pseudoIDs to the UKB IDs
bridge=pd.read_csv('pseudoid_to_spouseid.txt', sep=' ', header=None, names=['ID', 'pseudoID'])
bridge.index=bridge.pseudoID.to_list()

howe['IID']=howe.ID_15825.map(bridge.ID.to_dict())
howe=howe[~howe.IID.isnull()]
howe.IID=howe.IID.astype(int)

# Filter for complete pairs
spouses=howe.Couple.value_counts()
comp_pairs=spouses[spouses==2].index.to_list()
howe=howe[howe.Couple.isin(comp_pairs)]

# Make fam file for bam filtering
howe[['IID', 'IID']].to_csv('list_files/spouses.fam', sep=' ', index=False, header=None)

# Save all spouses as file
howe.to_csv('list_files/1_spouses.csv', index=False)
