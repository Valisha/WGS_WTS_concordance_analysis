#! /nfs/sw/python/python-3.6.1/bin/python3
import time
import os
import subprocess
import pandas as pd
from collections import defaultdict, Counter
import numpy as np

def parse_vcf_to_bed(input_dragen_vcf, sample_id):
    dragen_vcf = open(input_dragen_vcf, 'r')
    output_dragen_bed = open(str(sample_id+'.bed'), 'w')
    output_dragen_bed.write('\t'.join(['#CHR', 'START', 'END', 'INFO']) + '\n')
    lines = dragen_vcf.readlines()
    lines_new = []
    for line in lines:
        if '##' not in line and '#' not in line:
            line = line.strip('\n').split('\t')
            lines_new.append(line)
    dragen_table = pd.DataFrame(lines_new)
    dragen_table.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" ,"FORMAT", str(sample_id + 'T-D-W')]
    for i, v in enumerate(dragen_table['ID']):
        id_new = v.split(":")[3]
        start = id_new.split("-")[0]
        end = id_new.split("-")[1]
        if 'HET' in dragen_table.loc[i, 'INFO']:
            cn = int(np.round(float((dragen_table.loc[i, str(sample_id + 'T-D-W')]).split(":")[5])))
            if 'X' in dragen_table.loc[i, 'CHROM']:
                output_dragen_bed.write('\t'.join(['23', start, end, str(cn)]) + '\n')
            if 'X' not in dragen_table.loc[i, 'CHROM']:
                output_dragen_bed.write('\t'.join([dragen_table.loc[i, 'CHROM'][3:], start, end, str(cn)]) + '\n')
        elif 'HET' not in dragen_table.loc[i, 'INFO']:
            if 'X' in dragen_table.loc[i, 'CHROM']:
                cn = int((dragen_table.loc[i, str(sample_id + 'T-D-W')]).split(":")[1])
                output_dragen_bed.write('\t'.join(['23', start, end, str(cn)]) + '\n')
            if 'X' not in dragen_table.loc[i, 'CHROM']:
                cn = int((dragen_table.loc[i, str(sample_id + 'T-D-W')]).split(":")[1])
                output_dragen_bed.write('\t'.join([dragen_table.loc[i, 'CHROM'][3:], start, end, str(cn)]) + '\n')
    output_dragen_bed.close()


def run_bedtools(sample_id):
    call_string = str("bedtools intersect -b " + str(sample_id+'.bed') + " -a outcc.bed -wo > ../annotated_outputs_bedtools/" + str(sample_id + "_annotated_bedtools.tsv"))
    subprocess.call(call_string, shell=True)


def get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list):
    df_dragen = pd.read_csv(input_dragen, sep="\t")
    df_dragen.columns = ['CHR', 'START', 'END', 'GENE', 'CHR_dragen', 'START_dragen', 'END_dragen', 'CN_dragen', 'OVERLAP']
    df_dragen_new = df_dragen.loc[:,['CHR', 'START', 'END', 'GENE', 'CN_dragen']]
    gene_found = []
    dragen_cn_dict = {}
    for i, gene in enumerate(df_dragen_new['GENE']):
        if gene in genes_new:
            if gene not in gene_found and gene not in ignore_list:
                gene_found.append(gene)
                dragen_cn_dict.update({gene:[(df_dragen_new.loc[i, 'CN_dragen'])]})
            elif gene in gene_found and gene not in ignore_list:
                dragen_cn_dict[gene].append(df_dragen_new.loc[i, 'CN_dragen'])
    return dragen_cn_dict

def get_counter_per_gene_v5b_adjusted(input_ascat, genes_new, ignore_list):
    df_ascat= pd.read_csv(input_ascat, sep="\t")
    df_ascat_new = df_ascat.loc[:, ['Gene', 'CNV_copies']]
    gene_found = []
    ascat_cn_dict = {}
    for i, gene in enumerate(df_ascat_new['Gene']):
        if gene in genes_new:
            if ',' in df_ascat_new.loc[i, 'CNV_copies'] and gene not in gene_found and gene not in ignore_list:
                cn_list = df_ascat_new.loc[i, 'CNV_copies'].split(',')
                cn_list_round = [round(float(i)) for i in cn_list]
                gene_found.append(gene)
                ascat_cn_dict.update({gene:cn_list_round})
            if ',' not in df_ascat_new.loc[i, 'CNV_copies'] and gene not in gene_found and gene not in ignore_list:
                cn = df_ascat_new.loc[i, 'CNV_copies']
                cn_round = round(float(cn))
                gene_found.append(gene)
                ascat_cn_dict.update({gene:[cn_round]})
    return ascat_cn_dict


def compare(cn_dict, v5b_dict, sample_id):
    v5b_count = 0
    #print(v5b_dict)
    gene_found_in_dragen = []
    dragen_count = 0
    conc = open(str('concordance_report_v5b_adj_' + sample_id + '.txt'), 'w')
    conc.write("\t".join(["gene", "adjusted_CNV_v5b", "CN_dragen", "Concordant"]) + "\n")
    counter_conc = 0
    #print(cn_dict)
    found_keys = []
    for key, val in v5b_dict.items():
        v5b_count += len(val)
        if key in cn_dict.keys() and key not in found_keys:
            gene_found_in_dragen.append(key)
            if len(val) > 1:
                dragen_count += len(cn_dict[key])
                counter_val_v5b = Counter(v5b_dict[key])
                counter_val_dragen = Counter(cn_dict[key])
                counter_val_dragen = {int(k):int(v) for k,v in counter_val_dragen.items()}
                counter_val_v5b = {int(k):int(v) for k,v in counter_val_v5b.items()}
                if counter_val_dragen == counter_val_v5b:
                    counter_conc += len(cn_dict[key])
                    found_keys.append(key)
                    conc.write('\t'.join([key, str(counter_val_v5b), str(counter_val_dragen), 'Concordant']) + "\n")
                if v5b_dict[key] != cn_dict[key]:
                    counter_val_v5b = dict(Counter(v5b_dict[key]))
                    counter_val_dragen = dict(Counter(cn_dict[key]))
                    for key_v5b, val_v5b in counter_val_v5b.items():
                        if key_v5b in counter_val_dragen.keys() and key not in found_keys:
                            val_v5b = counter_val_v5b[key_v5b]
                            val_dragen = counter_val_dragen[key_v5b]
                            if int(counter_val_dragen[key_v5b]) == 0 and int(val_v5b) == 0:
                                counter_conc += val_v5b
                                found_keys.append(key)
                                conc.write('\t'.join([key, str(val_dragen), str(val_v5b), 'Concordant']) + "\n")
                            elif int(counter_val_dragen[key_v5b]) == 1 and int(val_v5b) == 1:
                                counter_conc += val_v5b
                                found_keys.append(key)
                                conc.write('\t'.join([key, str(val_dragen), str(val_v5b), 'Concordant']) + "\n")
                            elif int(counter_val_dragen[key_v5b]) == 2 and int(val_v5b) == 2:
                                counter_conc += val_v5b
                                found_keys.append(key)
                                conc.write('\t'.join([key, str(val_dragen), str(val_v5b), 'Concordant']) + "\n")
                            elif int(counter_val_dragen[key_v5b]) >= 3 and int(val_v5b) >= 3:
                                counter_conc += val_v5b
                                found_keys.append(key)
                                conc.write('\t'.join([key, str(val_dragen), str(val_v5b), 'Concordant']) + "\n")
                            elif int(counter_val_dragen[key_v5b]) != int(val_v5b):
                                if int(counter_val_dragen[key_v5b]) < int(val_v5b):
                                    counter_conc += val_dragen
                                    found_keys.append(key)
                                    conc.write('\t'.join([key, str(counter_val_v5b), str(counter_val_dragen), 'Partial-Concordant', str(val_dragen)]) + "\n")
                                elif int(val_dragen) > int(val_v5b):
                                    counter_conc += val_v5b
                                    found_keys.append(key)
                                    conc.write('\t'.join([key, str(counter_val_v5b), str(counter_val_dragen), 'Partial-Concordant', str(val_v5b)]) + "\n")
                        else:
                            conc.write('\t'.join([key, str(counter_val_v5b), str(counter_val_dragen), 'Discordant']) + "\n")
            if len(val) == 1:
                dragen_count += 1
                if int(val[0]) == 1 and int(cn_dict[key][0]) == 1:
                    counter_conc += 1
                    found_keys.append(key)
                    conc.write('\t'.join([key, str(val[0]), str(cn_dict[key][0]), 'Concordant']) + "\n")
                elif int(val[0]) == 0 and int(cn_dict[key][0]) == 0:
                    counter_conc += 1
                    found_keys.append(key)
                    conc.write('\t'.join([key, str(val[0]), str(cn_dict[key][0]), 'Concordant']) + "\n")
                elif int(val[0]) == 2 and int(cn_dict[key][0]) == 2:
                    counter_conc += 1
                    found_keys.append(key)
                    conc.write('\t'.join([key, str(val[0]), str(cn_dict[key][0]), 'Concordant']) + "\n")
                elif int(val[0]) >= 3 and int(cn_dict[key][0]) >= 3:
                    counter_conc += 1
                    found_keys.append(key)
                    conc.write('\t'.join([key, str(val[0]), str(cn_dict[key][0]), 'Concordant']) + "\n")
                else:
                    #print(key, val[0], cn_dict[key][0])
                    conc.write('\t'.join([key, str(val[0]), str(cn_dict[key][0]), 'Discordant']) + "\n")
    conc.write(str((counter_conc/v5b_count)*100) + "\n")
    print('percentage concordance ', str((counter_conc/(len(v5b_dict.keys())))*100))
    print('Number of dragen genes ', dragen_count)
    print('Number of v5b genes ',v5b_count )
    print(counter_conc)
    gene_not_found = v5b_dict.keys()-gene_found_in_dragen
    for g in gene_not_found:
        conc.write('\t'.join([g, 'not found in dragen']) + '\n')
    return str((counter_conc/(len(v5b_dict.keys())))*100), counter_conc, v5b_count
    
# write this into a file like gene, CN_dragen value if 1, CN_v5b value if 1 and then if there is a counter write it as a counter and all the concordance as the last line of the file 
if __name__ == "__main__":
    import time
    import os
    import subprocess
    import pandas as pd
    from collections import defaultdict, Counter
    start = time.time()
    census_genes = open("/gpfs/commons/groups/clinical/vshah/pipelines/vep_annotate/cosmic_actionable_old/cancer_gene_census.genes.txt", "r")
    census = census_genes.readlines()
    genes_new = []
    for line in census:
        gene = line.strip("\n")
        genes_new.append(gene)
    ignore_census_list = open("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/compbio/concordance/ignore-cancer-census.txt", "r")
    ignore_list = []
    ignore_lines = ignore_census_list.readlines()
    for gene in ignore_lines:
        gene = gene.strip("\n")
        ignore_list.append(gene)
    input_v5b = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/concordance/CA-0072T-D-W--CA-0072N-D-W.dna.rna.combined.adjLog2.v6.txt"
    input_dragen = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/annotated_outputs_bedtools/CA-0072_annotated_bedtools.tsv"
    input_dragen_vcf = "/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/CA-0072/CA-0072T-D-W--CA-0072N-D-W.cnv.vcf"
    df_v5b = pd.read_csv(input_v5b, sep="\t")
    parse_vcf_to_bed(input_dragen_vcf, 'CA-0072')
    run_bedtools('CA-0072')
    counter_dragen = get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list)
    counter_v5b = get_counter_per_gene_v5b_adjusted(input_v5b, genes_new, ignore_list)
    perc_concordance, dragen_count, v5b_count = compare(counter_dragen, counter_v5b, 'CA-0072')
