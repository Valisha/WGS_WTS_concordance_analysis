#! /nfs/sw/python/python-3.6.1/bin/python3
import pandas as pd
import sys
import subprocess
#input_data = sys.argv[1]

#df = pd.read_csv('/gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv', sep='\t')
#print(df.columns)


mskcc_dragen_snv = open('/gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv', 'r')
mskcc_dragen = mskcc_dragen_snv.readlines()
counter = 0
sample_ids = []
gene_found = []
gene_list = []
line_pass = 0
line_total = 0
out1 = open('final1.tsv', 'w')
for line in mskcc_dragen:
    line_new = line.strip('\n').split('\t')
    if len(line_new) > 1:
        sample_id = line_new[0]
        chrs = line_new[1]
        if 'chr' not in line:
            #counter += 1
            #if sample_id not in sample_ids:
                #sample_ids.append(sample_id)
                #call_string = 'more /gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv | grep "chr" | grep "' + sample_id + '" | wc -l'
                #call_string1 = 'more /gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv | grep "chr" | grep "' + sample_id + '" | grep "PASS" | wc -l'
                #print(call_string1)
                #process_new1 = subprocess.Popen(call_string1, stdout=subprocess.PIPE, shell=True)
                #process2 = subprocess.Popen(call_string, stdout=subprocess.PIPE, shell=True)
                #subprocess.call(call_string, shell=True)
                #subprocess.call(call_string1, shell=True)
                #temp_out3, temp_err3 = process2.communicate()
                #temp_out4, temp_err4 = process_new1.communicate()
                #new4 = temp_out4.decode(sys.stdout.encoding).strip()
                #new3 = temp_out3.decode(sys.stdout.encoding).strip()
                #line_total += int(new3)
                #line_pass += int(new4)
                #out1.write('\t'.join([sample_id, new3, new4]) + '\n')

            gene = line_new[1]
            if gene not in gene_list and 'chr' not in line:
                if ' ' in gene:
                    gene = gene.split(' ')[0]
                gene_list.append(gene)
                for l in mskcc_dragen:
                    line_dragen = l.strip('\n').split('\t')
                    if len(line_dragen) > 1:
                        chrs_l = line_dragen[1]
                    else:
                        chrs_l = ''
                    if 'chr' in chrs_l:
                        #print(line_dragen)
                        if len(line_dragen) > 1:
                            gene_dragen = line_dragen[8]
                            if ',' in gene_dragen:
                                gene_new = gene_dragen.split(',')
                                for g in gene_new:
                                    #print(g, gene)
                                    if g != 'NA':
                                        if gene == g:
                                            #print(g)
                                            line_total += 1
                                            if g not in gene_found:
                                                line_total += 1
                                                gene_found.append(g)
                            else:
                                if gene == gene_dragen:
                                    #print(gene_dragen)
                                    if gene not in gene_found:
                                        gene_found.append(gene_dragen)
                                        line_total += 1
print(line_total)
#print(counter)

#for gene_new in gene_list:
#    if gene_new not in gene_found:
#        call_string = 'more /gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv | grep "chr" | grep "' + gene_new + '" | grep -v PASS'
#        process1 = subprocess.Popen(call_string, stdout=subprocess.PIPE, shell=True)
#        subprocess.call(call_string, shell=True)
#        temp_out1, temp_err1 = process1.communicate()
#        new2 = temp_out1.decode(sys.stdout.encoding).strip()
#        print(new2)
#out1.write(str(line_total) +"\t" + str(line_pass) + '\n')
#counter_gene = 0
#new = open('/gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-positive.tsv', 'r')
#counter_gene_pass = 0
out = open('final.tsv', 'w')
#for gene in gene_list:
    #call_string = 'more /gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-test.tsv | grep "chr" | grep "' + gene + '" | grep -v "PASS" | wc -l'
    #call_string_new = 'more /gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC-positive.tsv | grep "chr" | grep "' + gene + '" | wc -l'
    #process1 = subprocess.Popen(call_string, stdout=subprocess.PIPE, shell=True)
    #process_new = subprocess.Popen(call_string_new, stdout=subprocess.PIPE, shell=True)
    #subprocess.call(call_string, shell=True)
    #subprocess.call(call_string_new, shell=True)
    #temp_out1, temp_err1 = process1.communicate()
    #temp_out2, temp_err2 = process_new.communicate()
    #new2 = temp_out1.decode(sys.stdout.encoding).strip()
    #new = temp_out2.decode(sys.stdout.encoding).strip()
    #counter_gene += int(new2)
    #counter_gene_pass += int(new)
    #out.write('\t'.join([gene, new]) + '\n')
#out.write('\t'.join([str(counter_gene), str(counter_gene_pass)]) + '\n')
gene_dragen_list_new = []
with open('/gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/gene_dragen.txt', 'r') as gene_dragen_list:
    for gene in gene_dragen_list:
        gene = gene.strip('\n')
        #if ',' in gene:
        #    gene1 = gene.split(',')
        #    for g1 in gene1:
        #        if g1 != 'NA':
        #            gene_dragen_list_new.append(g1)
        #else:
        gene_dragen_list_new.append(gene)
gene_dragen_pass = []
with open('/gpfs/commons/groups/clinical/vshah/concordance_test_snv_RES/workflow-2/full_concordance_dragen_MSKCC_dragen-PASS.txt', 'r') as gene_dragen_PASS:
    for gene in gene_dragen_PASS:
        gene = gene.strip('\n')
        #if ',' in gene:
        #    gene1 = gene.split(',')
        #    for g1 in gene1:
        #        if g1 != 'NA':
        #            gene_dragen_pass.append(g1)
        #else:
        gene_dragen_pass.append(gene)
print(len(gene_dragen_list_new))
print(len(gene_dragen_pass))
from collections import Counter
gene_dragen_dict = Counter(gene_dragen_list_new)
gene_dragen_dict_PASS = Counter(gene_dragen_pass)
total = 0
for key, val in gene_dragen_dict.items():
    out.write('\t'.join([str(key), str(val), str(gene_dragen_dict_PASS[key])]) + '\n')
    #print('\t'.join([str(key), str(val), str(gene_dragen_dict_PASS[key])]) + '\n')
    total += val
total_pass = 0
for key, val in gene_dragen_dict_PASS.items():
    total_pass += val
print(total)
print(total_pass)
