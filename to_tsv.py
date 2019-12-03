import os, subprocess, sys

outputfolder = sys.argv[1]

probe_gene = 'testing_data/probe_info'

AJHG_CSH_probe_gene = {}
CU_probe_gene = {}
DGD_probe_gene = {}

for line in open(probe_gene, 'r'):
    data = line.rstrip('\n').split('\t')
    if(data[0] == 'AJHG'):
        AJHG_CSH_probe_gene[data[1]] = data[2].strip(' ')
    elif(data[0] == 'CSH'):
        AJHG_CSH_probe_gene[data[1]] = data[2].strip(' ')
    elif(data[0] == 'Columbia_U'):
        CU_probe_gene[data[1]] = data[2].strip(' ')
    elif(data[0] == 'DGD'):
        DGD_probe_gene[data[1]] = data[2].strip(' ')


AJHG_total_num = 83
CSH_total_num = 72
CU_total_num = 27
DGD_total_num = 85
TAF1_total_num = 14

AJHG_old_phenolyzer = [0,0,3.6,7.2,7.2,7.2,7.2,7.2]
CSH_old_phenolyzer = [20.8,23.6,23.6,25,25,25,25]
AJHG_CSH_old_phenolyzer = [9.0, 11.0, 19.4, 23.2, 23.2, 26.5, 41.9]
CU_old_phenolyzer  = [7.4, 7.4, 18.5, 37.0, 37.0, 37.0, 63.0]
DGD_old_phenolyzer  =  [10.6, 10.6, 22.4, 28.2, 28.2, 28.2, 44.7]
TAF1_old_phenolyzer = [85.7,85.7,100,100,100,100,100]

AJHG_CSH_tops = [0,0,0,0,0,0,0]
CU_tops = [0,0,0,0,0,0,0]
DGD_tops = [0,0,0,0,0,0,0]
TAF1_tops = [0,0,0,0,0,0,0]

def tops(set_tops,rank):
    if(rank <= 10):
        set_tops[0] += 1
    if(rank <= 25):
        set_tops[1] += 1
    if(rank <= 50):
        set_tops[2] += 1
    if(rank <= 100):
        set_tops[3] += 1
    if(rank <= 250):
        set_tops[4] += 1
    if(rank <= 500):
        set_tops[5] += 1
    if(rank <= 1000):
        set_tops[6] += 1
    return set_tops


for c in os.listdir(outputfolder):
    if(c.startswith('PMID')):
        gene = AJHG_CSH_probe_gene[c]
        result = subprocess.run(['grep', gene, '-w','{}/{}/out.final_gene_list'.format(outputfolder, c)],stdout=subprocess.PIPE)
        value = result.stdout.decode().split('\t')[0]
        if(value != ''):
            AJHG_CSH_tops = tops(AJHG_CSH_tops, int(value))
        
    if(c.startswith('sample')):
        gene = DGD_probe_gene[c]
        result = subprocess.run(['grep', gene, '-w','{}/{}/out.final_gene_list'.format(outputfolder, c)],stdout=subprocess.PIPE)
        value = result.stdout.decode().split('\t')[0]
        if(value != ''):
            DGD_tops = tops(DGD_tops, int(value))

        

    if(c.startswith('case')):
        gene = 'TAF1'
        result = subprocess.run(['grep', gene, '-w','{}/{}/out.final_gene_list'.format(outputfolder, c)],stdout=subprocess.PIPE)
        value = result.stdout.decode().split('\t')[0]
        if(value != ''):
            TAF1_tops = tops(TAF1_tops, int(value))
    


    if(c.startswith('Columbia')):
        gene = CU_probe_gene[c]
        result = subprocess.run(['grep', gene, '-w','{}/{}/out.final_gene_list'.format(outputfolder, c)],stdout=subprocess.PIPE)
        value = result.stdout.decode().split('\t')[0]
        if(value != ''):
            CU_tops = tops(CU_tops, int(value))
for i in range(len(AJHG_CSH_tops)):
    AJHG_CSH_tops[i] = round( AJHG_CSH_tops[i]/(AJHG_total_num + CSH_total_num) * 100, 1)

for i in range(len(CU_tops)):
    CU_tops[i] = round( (CU_tops[i]/(CU_total_num) * 100), 1)

for i in range(len(TAF1_tops)):
    TAF1_tops[i] = round( TAF1_tops[i]/(TAF1_total_num) * 100, 1)


for i in range(len(DGD_tops)):
    DGD_tops[i] = round( DGD_tops[i]/(DGD_total_num) * 100, 1)

def output_tsv(new, old, f_name):
    os.makedirs('rankings',exist_ok=True)
    with open('rankings/'+f_name + '.tsv', 'w+') as fw:
        fw.write('\tnew\told Phenolyzer\n')
        fw.write('Top 10\t{}\t{}\n'.format(str(new[0]), str(old[0]) ))
        #fw.write('Top 25\t{}\t{}\n'.format(str(new[1]), str(pheno[1]) ))
        fw.write('Top 50\t{}\t{}\n'.format(str(new[2]), str(old[2]) ))
        fw.write('Top 100\t{}\t{}\n'.format(str(new[3]), str(old[3]) ))
        #fw.write('Top 250\t{}\t{}\n'.format(str(p2g[4]), str(pheno[4]) ))
        #fw.write('Top 500\t{}\t{}\n'.format(str(p2g[5]), str(pheno[5]) ))
        fw.write('Top 1000\t{}\t{}\n'.format(str(new[6]), str(old[6]) ))

output_tsv(AJHG_CSH_tops,AJHG_CSH_old_phenolyzer, 'AJHG_CSH')


output_tsv(DGD_tops,DGD_old_phenolyzer, 'DGD')
output_tsv(TAF1_tops,TAF1_old_phenolyzer, 'TAF1')
output_tsv(CU_tops,CU_old_phenolyzer, 'CU')

