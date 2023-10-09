#!/bin/python3

# Remove a leading . from the INFO field in VCF
in_vcf='vcfs/all_annovar.hg38_multianno.vcf'
fin = open(in_vcf, 'r')
output=open('vcfs/7_remove_dot.vcf', 'w')

for line in fin:
	if line.startswith('#'):
		output.write(line)
		continue
	sline = line.split('\t')
	info = sline[7]
	if info.startswith('.;'):
		info = info[2:]
	sline[7] = info
	line = '\t'.join(sline)
	output.write(line)

fin.close()
output.close()
