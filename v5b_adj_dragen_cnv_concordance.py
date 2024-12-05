#! /nfs/sw/python/python-3.6.1/bin/python3

if __name__ == "__main__":
    import time
    import os
    import subprocess
    import pandas as pd
    import adjusted_log2_dragen_conc_test_HET as conc
    from collections import defaultdict, Counter
    start = time.time()
    counter_file = 0
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
    comp = open("compiled_cnv_concordance_v5b_adj_ngs.tsv", "w")
    comp.write("\t".join(["sample_id", "dragen_count", "v5b_count", "percentage_concordance"]) + "\n")
    dragen = os.listdir("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/")
    ascat = open('/gpfs/commons/groups/clinical/vshah/somatic_cnv/concordance/v5b_cn_adjusted.txt', 'r')
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
                conc.parse_vcf_to_bed(input_dragen_vcf, file_new)
                conc.run_bedtools(file_new)
                counter_dragen = conc.get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list)
                counter_v5b = conc.get_counter_per_gene_v5b_adjusted(input_ascat, genes_new, ignore_list)
                perc_concordance, counter_conc, v5b_count = conc.compare(counter_dragen, counter_v5b, file_new)
                perc.append(perc_concordance)
                comp.write('\t'.join([file_new, str(counter_conc), str(v5b_count), str((counter_conc/v5b_count)*100)]) + '\n')
                counter_file += 1
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
                conc.parse_vcf_to_bed(input_dragen_vcf, file_new)
                conc.run_bedtools(file_new)
                counter_dragen = conc.get_counter_per_gene_dragen(input_dragen, genes_new, ignore_list)
                counter_v5b = conc.get_counter_per_gene_v5b_adjusted(input_ascat, genes_new, ignore_list)
                perc_concordance, counter_conc, v5b_count = conc.compare(counter_dragen, counter_v5b, file_new)
                perc.append(perc_concordance)
                comp.write('\t'.join([file_new, str(counter_conc), str(v5b_count), str(perc_concordance)]) + '\n')
                counter_file += 1
    stop = time.time()
    print(counter_file)
    print(stop-start)
