import pandas as pd
import subprocess
import datetime

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Variant calls
var_file = pd.read_csv('tables/1_combined_calls.csv')

# PED file
ped = pd.read_csv('tables/spark_fams.csv')

# gVCF locations
vcfs = pd.read_csv('tables/samples.csv')

def find_in_vcf(vcf, chr, pos, alt):
	chr_num=chr.split('r')[1]
	# Use bcftools to get relevant line in parent file
	command = "bcftools view -r %s:%s -H %s" % (chr_num, pos, vcf)
	rel_line = subprocess.run(command, capture_output = True, shell = True).stdout.decode()

	# If no matching line is found
	if len(rel_line) == 0:
		print('No matching line found! - ', command)
		return [-1, -1]

	line = rel_line.split('\t')
	format = line[8].split(':')
	info = line[9].split(':')
	# Look for allele depth information
	# NOTE: AD will only be reported in cases where parent has non-reference allele
	# Check if 'AD' is reported
	if 'AD' in format:
		ad_loc = format.index('AD')
		ad = info[ad_loc]
		ad = [int(i) for i in ad.split(',')]
		# Be careful of multi-allelic sites!
		# Check the length of AD
		# If it is more than 2, only return the read depth information for the reference allele and specific alternate allele
		# Alleles are always in the order Ref, Alt1, Alt2, Alt3, . . . , <NON_REF>
		if len(ad) > 2:
			alts = line[4].split(',')
			# Check if the alternate allele is in the list
			# If not, return the AD for <NON_REF> (should be 0)
			if alt in alts:
				alt_loc = alts.index(alt)+1
			else:
				print(alts)
				print(alt)
				alt_loc = alts.index('<NON_REF>')+1
			ad = [ad[0], ad[alt_loc]]
		return ad

	# If AD is not in parent file, they only have reference allele
	# Confirm this is true by checking genotype info
	geno_loc = format.index('GT')
	gt = info[geno_loc]

	if gt != '0|0':
		print('Parent has non-reference allele but is missing AD information!', command)
		return [-1, -1]
	# If they are homozygous for the reference, return an alt allele depth of 0 and the read depth as the depth of the reference
	dp_loc = format.index('DP')
	dp = info[dp_loc]
	return [int(dp), 0]

def interp_ad(ad):
	# To be called as de novo, read depth in parents needs to be at least 15
	# Additionally, the number of alt reads in the parents should be 0
	# To be called as inherited, at least 20% of the reads in the parent need to be the alternate allele
	alt_read = int(ad[1])
	ref_read = int(ad[0])

	parent=False
	unknown=False

	if alt_read==0 and ref_read==0:
		unknown = True
	elif float(alt_read)/(alt_read+ref_read) >= 0.2:
		parent = True
	elif alt_read > 0:
		unknown = True
	if sum([ref_read, alt_read]) < 15:
		unknown = True

	return (parent, unknown)


def add_inherit(row, father = False, mother = False, unknown = False):
	if int(row.name)%10000 == 0:
		print('Started line ' + str(row.name) + " at: " + str(datetime.datetime.now()))

	iid = row['Sample']

	# Get parents
	father_code = ped[ped.spid==iid]['father'].to_string(index = False, header = False).strip()
	mother_code = ped[ped.spid==iid]['mother'].to_string(index = False, header = False).strip()

	# Check if we have WGS data available for the parents
	if mother_code not in vcfs.Sample.to_list() and father_code not in vcfs.Sample.to_list():
		return '.'
	if mother_code == '0' and father_code == '0':
		return '.'


	# Get carrier status of parents
	if mother_code == '0':
		mother_return = '0'
	else:
		mother_return='M'

	if father_code == '0':
		father_return = '0'
	else:
		father_return='F'

	# Check if parents have the same variant
	vid = row['variant_id']
	other_samps = var_file[var_file.variant_id == vid]['Sample'].to_list()

	if mother_code in other_samps:
		mother = True
	if father_code in other_samps:
		father = True

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# If neither parent has variant, see if variant was present in original VCF
	chr = row['Chr']
	pos = row['Pos']
	alt = row['Alt']

	unknown1=False
	unknown2=False
	if mother_code != '0' and mother_code in vcfs.Sample.to_list():
		mother_vcf = vcfs[vcfs.Sample == mother_code]['VCF_filename'].to_string(index = False, header = False).strip()
		ad = find_in_vcf(mother_vcf, chr, pos, alt)
		if ad:
			mother, unknown1 = interp_ad(ad)

	if father_code != '0' and father_code in vcfs.Sample.to_list():
		father_vcf = vcfs[vcfs.Sample == father_code]['VCF_filename'].to_string(index = False, header = False).strip()
		ad = find_in_vcf(father_vcf, chr, pos, alt)
		if ad:
			father, unknown2 = interp_ad(ad)

	unknown=False
	if unknown1 or unknown2:
		unknown=True

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# To be called de novo, make sure we have WGS available for both parents
	if mother_code=='0' or father_code=='0':
		unknown = True

	if mother_code not in vcfs.Sample.to_list() or father_code not in vcfs.Sample.to_list():
		unknown = True

	if unknown:
		return '.'
	else:
		return 'de novo'

var_file['inheritance'] = var_file.apply(lambda row: add_inherit(row), axis = 1)

var_file.to_csv('SPARK_Variants_wInheritance.txt', sep = '\t', index = False)
