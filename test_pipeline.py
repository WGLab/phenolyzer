import os

test_set_folder = './test_sets/'

test_set_out = 'testsetoutput/'


for tset in os.listdir(test_set_folder):
	for case in os.listdir(test_set_folder+ tset + '/'):
		cmd = 'perl disease_annotation.pl {} -f -p -ph -logistic -out {}{}_{}/out'.format(test_set_folder+ tset + '/' + case, test_set_out, tset, case)
		qsub_cmd = 'echo "{}" | qsub -cwd -V -e {}.e -o {}.o'.format(cmd, tset+'_' + case, tset + '_'+ case)
		print(qsub_cmd)
		os.system(qsub_cmd)

