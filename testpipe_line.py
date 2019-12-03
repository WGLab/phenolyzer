import os, sys


output = sys.argv[1]

os.system('mkdir -p {}'.format(output))

for s in os.listdir('testing_data'):
	if(not os.path.isdir('{}/{}'.format('testing_data',s))):
		print(s)
		continue

	for c in os.listdir('{}/{}'.format('testing_data', s)):
		cmd = 'perl disease_annotation.pl {} -f -p -ph -logistic -out {}/{}/out'.format('{}/{}/{}'.format('testing_data', s,c), output, c)
		qsub = 'echo "{}" | qsub -cwd -V -l h_vmem=6G -e {}.e -o {}.o'.format(cmd, c,c)
		print(qsub)
		os.system(qsub)



