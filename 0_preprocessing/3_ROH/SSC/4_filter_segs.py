import pandas as pd

# Filter homozygous regions to remove those overlapping CNVs
hom=pd.read_csv('Result_files/3.1_1kb_window.hom', sep='\s+')
hom=hom[['IID', 'CHR', 'POS1', 'POS2']]
hom['CHR']='chr'+hom['CHR'].astype(str)

cnv=pd.read_csv('PennCNV/tables/7_size_filtered.txt', sep='\t')
cnv=cnv[cnv.Type=='DEL']

print(hom)

out=[]
for idx,row in hom.iterrows():
	# Check if sample has any CNVs in the same region as the ROH
	chr=row['CHR']
	start=int(row['POS1'])
	end=int(row['POS2'])

	samp=row['IID']

	overlap=cnv[(cnv.Sample==samp) & (cnv.Chr==chr) & (((cnv['Start']>=start) & (cnv['Start']<=end)) | ((cnv['End']>=start) & (cnv['End']<=end)) | ((cnv['End']>=end) & (cnv['Start']<=start)))].copy()
	if overlap.shape[0]==0:
		out.append(row.to_list())
		continue

	overlap['ROH_start']=start
	overlap['ROH_end']=end

	# Remove unneeded columns
	overlap=overlap[['Chr', 'Start', 'End', 'ROH_start', 'ROH_end']]

	# Ensure dels go in ascending order
	overlap.sort_values(by='Start', inplace=True)

	# If an ROH has multiple dels, check if the distance between the dels is < 1Mb
	for oidx, orow in overlap.iterrows():
		for oidx2, orow2 in overlap.iterrows():
			if oidx2 <= oidx:
				continue
			if orow2['Start']-orow['End']<=1000000:
				# Merge close-together calls
				overlap.at[oidx, 'End']=orow2['End']
				overlap.at[oidx2, 'Start']=orow['Start']

	overlap.drop_duplicates(inplace=True)

	# Check if deletion is in the middle of the ROH or on the outsides
	# For our purposes "middle" is defined as more than 1 Mb away from each end (we only need to keep the edges if they are more than 1Mb)
	overlap['edge']=''
	overlap.loc[(overlap.Start>=start) & (overlap.End>=(end-1000000)), 'edge']='right'
	overlap.loc[(overlap.Start<=(start+1000000)) & (overlap.End<=end), 'edge']='left'
	overlap.loc[(overlap.Start<=(start+1000000)) & (overlap.End>=(end-1000000)), 'edge']='entire'
	overlap.loc[(overlap.Start>=(start+1000000)) & (overlap.End<=(end-1000000)), 'edge']='middle'
	overlap.edge=pd.Categorical(overlap.edge, ['left', 'right', 'middle', 'entire'])

	# If any CNV completely overlaps the ROH, remove the ROH
	if 'entire' in overlap.edge.to_list():
		continue

	overlap.sort_values(by='edge', inplace=True)

	# Get the overlap
	starts=[0]*10
	ends=[0]*10
	starts[0]=start
	ends[0]=end
	loc=0
	for idx2, row2 in overlap.iterrows():
		edge=row2['edge']
		if edge=='middle':
			# First, ensure deletion is still in the middle with new start and end
			if row2['Start']-starts[loc]<=1000000:
				# Now on the left
				if row2['End']>starts[loc]:
					starts[loc]=row2['End']
			elif ends[loc]-row2['End']<=1000000:
				# Now on the right
				if row2['Start']<ends[loc]:
					ends[loc]=row2['Start']
			else:
				# If a deletion is in the middle of a segment, it will make 2 new segments
				starts[loc+1]=row2['End']
				ends[loc+1]=ends[loc]
				ends[loc]=row2['Start']
				loc+=1
		if edge=='left':
			# If deletion is on the left edge, the deletion end becomes the new start
			if row2['End']>starts[loc]:
				starts[loc]=row2['End']
		if edge=='right':
			# If deletion is on the right edge, the deletion start becomes the new end
			if row2['Start']<ends[loc]:
				ends[loc]=row2['Start']

	if overlap[overlap.edge=='middle'].shape[0]>1:
		print(overlap)
		print(starts, ends)

	for i in range(len(starts)):
		if ends[i]-starts[i]>=1000000:
			out.append([samp, chr, starts[i], ends[i]])

new_roh=pd.DataFrame(out, columns=['Sample', 'Chr', 'Start', 'End'])
new_roh['Length']=new_roh.End-new_roh.Start

# Save to file
new_roh.to_csv('Result_files/4_ROH_segs_filtered.csv', index=False)

# Also calculate summary values for each person
samps=list(hom.IID.unique())
samps.sort()
summ_lst=[]
for s in samps:
	nseg=new_roh[new_roh.Sample==s].shape[0]
	kb=sum(new_roh[new_roh.Sample==s]['Length'].to_list())/1000
	summ_lst.append([s, nseg, kb])
summ_df=pd.DataFrame(summ_lst, columns=['Sample', 'NSEG', 'KB'])
summ_df['KB_AVG']=summ_df.KB/summ_df.NSEG
print(summ_df)

# Save to file
summ_df.to_csv('Result_files/4_ROH_summary.csv', index=False)
