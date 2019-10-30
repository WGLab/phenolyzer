import os
import sys


test_folders = []

test_data=sys.argv[1]

for i in os.listdir(test_data):
        test_folders.append(test_data+'/'+i+ '/')
print(test_folders)
for i in test_folders:
        name = i.split('/')[1]
        for j in os.listdir(i):
                cmd = 'perl disease_annotation.pl -f {} -hpo_gene_weight 0.1 -prediction -phenotype -logistic -nproc 12 -out test/{}/out'.format(i+j, j)
                qsub_cmd = 'echo "{}" | qsub -cwd -V -l h_vmem=2G -pe smp 3 -e {}{} -o {}{}'.format(cmd, 'err/', j, 'log_ppln/', j)
                os.system(qsub_cmd)


