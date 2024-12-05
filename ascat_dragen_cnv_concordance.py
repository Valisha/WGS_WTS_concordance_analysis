#! /nfs/sw/python/python-3.6.1/bin/python3

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
        cn = (dragen_table.loc[i, str(sample_id + 'T-D-W')]).split(":")[1]
        output_dragen_bed.write('\t'.join([(dragen_table.loc[i, 'CHROM'])[3:], start, end, cn]) + '\n')
    output_dragen_bed.close()


def run_bedtools(sample_id):
    call_string = str("bedtools intersect -b " + str(sample_id+'.bed') + " -a outcc.bed -wo > ../annotated_outputs_bedtools/" + str(sample_id + "_annotated_bedtools.tsv"))
    subprocess.call(call_string, shell=True)

# write this into a file like gene, CN_dragen value if 1, CN_v5b value if 1 and then if there is a counter write it as a counter and all the concordance as the last line of the file 
if __name__ == "__main__":
    import time
    import os
    import subprocess
    import pandas as pd
    import conc_test_gene_ascat_2 as conc
    from collections import defaultdict, Counter
    start = time.time()
    census_genes = open("/gpfs/commons/groups/clinical/vshah/pipelines/vep_annotate/cosmic_actionable_old/cancer_gene_census.genes.txt", "r")
    census = census_genes.readlines()
    genes_new = []
    ignore_census_list = open("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/compbio/concordance/ignore-cancer-census.txt", "r")
    ignore_list = []
    ignore_lines = ignore_census_list.readlines()
    for gene in ignore_lines:
        gene = gene.strip("\n")
        ignore_list.append(gene)
    for line in census:
        gene = line.strip("\n")
        genes_new.append(gene)
    start = time.time()
    comp = open("compiled_cnv_concordance_ascat_ngs.tsv", "w")
    comp.write("\t".join(["sample_id", "dragen_count", "v5b_count", "percentage_concordane"]) + "\n")
    dragen = os.listdir("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/")
    ascat = open('/gpfs/commons/groups/clinical/vshah/somatic_cnv/concordance/ascat_ngs_list.txt', 'r')
    ascat_files = ascat.readlines()
    perc = []
    sample_ids = []
    for file_dragen in dragen:
        sample_ids.append(file_dragen)
    for file_ascat in ascat_files:
        if "Sample_" in file_ascat:
            file_ascat = file_ascat.strip('\n')
            counter = 0
            file1 = file_ascat.split(".")[0]
            file_new = file1.split("T")[0]
            file_new = file_new.split("Sample_")[1]
            if file_new in sample_ids:
                print(file_new)
                input_ascat = file_ascat
                input_dragen = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/annotated_outputs_bedtools/" + file_new + "_annotated_bedtools.tsv"
                input_dragen_vcf = "/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/" + file_new + "/" + file_new + "T-D-W--" + file_new + "N-D-W.cnv.vcf"
                df_v5b = pd.read_csv(input_ascat, sep="\t")
                parse_vcf_to_bed(input_dragen_vcf, file_new)
                run_bedtools(file_new)
                counter_dragen = conc.get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list)
                counter_v5b = conc.get_counter_per_gene_ascat(input_ascat, genes_new, ignore_list)
                perc_concordance, counter_conc, v5b_count = conc.compare(counter_dragen, counter_v5b, file_new)
                perc.append(perc_concordance)
                comp.write('\t'.join([file_new, str(counter_conc), str(v5b_count), str(perc_concordance)]) + '\n')
        if "Sample_" not in file_ascat:
            file_ascat = file_ascat.strip("\n")
            counter = 0
            file1 = file_ascat.split("/")[-1]
            file_new= file1.split("T")[0]
            if file_new in sample_ids:
                print(file_new)
                input_ascat = file_ascat
                input_dragen = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/annotated_outputs_bedtools/" + file_new + "_annotated_bedtools.tsv"
                input_dragen_vcf = "/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/" + file_new + "/" + file_new + "T-D--" + file_new + "N-D.cnv.vcf"
                try:
                    subprocess.call(str("gunzip " + input_dragen_vcf + ".gz"), shell=True)
                except:
                    pass
                df_v5b = pd.read_csv(input_ascat, sep="\t")
                parse_vcf_to_bed(input_dragen_vcf, file_new)
                run_bedtools(file_new)
                counter_dragen = conc.get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list)
                counter_v5b = conc.get_counter_per_gene_ascat(input_ascat, genes_new, ignore_list)
                perc_concordance, counter_conc, v5b_count = conc.compare(counter_dragen, counter_v5b, file_new)
                perc.append(perc_concordance)
                comp.write('\t'.join([file_new, str(counter_conc), str(v5b_count), str(perc_concordance)]) + '\n')
    stop = time.time()
    print(stop-start)
