
outfile=open('0_model_parameters.txt', 'w')

h2='0.5'
corr='0.29'
n_loci='10000'
common_prop='0.9'
effect_ratio='1'
pop_n='100000'

outfile.write(' '.join([h2, corr, n_loci, common_prop, effect_ratio, pop_n])+'\n')

for i in range(2, 9):
	if i!=5:
		outfile.write(' '.join([str(i/10), corr, n_loci, common_prop, effect_ratio, pop_n])+'\n')

for i in range(20, 50, 5):
	outfile.write(' '.join([h2, str(i/100), n_loci, common_prop, effect_ratio, pop_n])+'\n')

for i in [2,3]:
	outfile.write(' '.join([h2, corr, str(10**i), common_prop, effect_ratio, pop_n])+'\n')

for i in range(1, 9):
	outfile.write(' '.join([h2, corr, n_loci, str(i/10), effect_ratio, pop_n])+'\n')

for i in [0.25, 0.5, 2, 4]:
	outfile.write(' '.join([h2, corr, n_loci, common_prop, str(i), pop_n])+'\n')

outfile.write(' '.join([h2, corr, n_loci, common_prop, effect_ratio, str(10**4)])+'\n')

outfile.close()
