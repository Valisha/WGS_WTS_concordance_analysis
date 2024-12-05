#! /nfs/sw/python/python-3.6.1/bin/python3
import pandas as pd
import sys
import os
import time 


def main_function(input_v5b, input_dragen, sample_id, genes_new):
    df_v5b = pd.read_csv(input_v5b, sep="\t", header=None)
    df_dragen = pd.read_csv(input_dragen, sep=",")

    df_v5b_new = df_v5b.iloc[:,[4, 1]]
    df_v5b_new.columns = ['transcript_id', 'log2_values']

    df_dragen_new = df_dragen.loc[:,['cytoband', 'copyNumber', 'FILTER', 'cancer_census_gene']]
    conc = open(str('concordance_report_' + sample_id + '.txt'), 'w')
    conc.write("\t".join(["cytoband_v5b", "log2values", "cytoband_dragen" , "CN_dragen", "FILTER", "Concordant"]) + "\n")
    cyto_found = []
    counter = 0
    gene_found = 0
    for i, v in enumerate(df_dragen_new['cytoband']):
        if "p" in v:
            v_new = v.split("p")
            cyto1 = v_new[1].split("-")[0]
            cytoband_new1 = "p" + cyto1
            for a, b  in enumerate(df_v5b_new['transcript_id']):
                if b == cytoband_new1 and b not in cyto_found:
                    cyto_found.append(b)
                    genes = df_dragen_new['cancer_census_gene'][i].split(',')
                    for gene in genes:
                        if gene in genes_new:
                            gene_found += 1
                            if float(df_v5b_new['log2_values'][a]) > 0.1 and float(df_dragen_new['copyNumber'][i]) == 3 and gene_found == 1:
                                #print(gene, str(df_v5b_new['log2_values'][a]))
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                            if float(df_v5b_new['log2_values'][a]) < 0.1 and float(df_dragen_new['copyNumber'][i]) < 2 and gene_found == 1:
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                            if 0.1 < float(df_v5b_new['log2_values'][a]) < 0.1 and float(df_dragen_new['copyNumber'][i]) == 2 and gene_found == 1:
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                    gene_found = 0
        if "q" in v:
            v_new = v.split("q")
            cyto1 = v_new[1].split("-")[0]
            cytoband_new1 = "q" + cyto1
            for a,b  in enumerate(df_v5b_new['transcript_id']):
                if b == cytoband_new1 and b not in cyto_found:
                    cyto_found.append(b)
                    genes = df_dragen_new['cancer_census_gene'][i].split(',')
                    for gene in genes:
                        if gene in genes_new:
                            gene_found += 1
                            if float(df_v5b_new['log2_values'][a]) > 0.1 and float(df_dragen_new['copyNumber'][i]) == 3 and gene_found == 1:
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                            if float(df_v5b_new['log2_values'][a]) < 0.1 and float(df_dragen_new['copyNumber'][i]) < 2 and gene_found == 1:
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                            if 0.1 < float(df_v5b_new['log2_values'][a]) < 0.1 and float(df_dragen_new['copyNumber'][i]) == 2 and gene_found == 1:
                                counter += 1
                                conc.write('\t'.join(df_v5b_new.iloc[a,:].map(str)) + '\t' + '\t'.join(df_dragen_new.iloc[i, :].map(str)) + '\t' + 'Concordant' + '\n')
                    gene_found = 0

    for i, v in enumerate(df_v5b_new['transcript_id']):
        if v not in cyto_found:
            conc.write('\t'.join(df_v5b.iloc[i,:].map(str)) + '\t' + 'Non-cordant' + '\n')
    percentage_concordance = (counter/(float(len(df_v5b_new['transcript_id']))))*100
    conc.write(str(percentage_concordance) + '\n')
    print(counter)
    print(float(len(df_v5b_new['transcript_id'])))
    return counter
if __name__ == "__main__":
    start = time.time()
    comp = open("compiled_cnv_concordance.tsv", "w")
    comp.write("\t".join(["sample_id", "dragen_count", "v5b_count", "percentage_concordane", "number_concordant"]) + "\n")
    dragen = os.listdir("/gpfs/commons/groups/clinical/vshah/somatic_cnv/output_csv_with_census")
    v5b = os.listdir("/gpfs/commons/groups/clinical/vshah/somatic_cnv/v5b_outputs")
    perc = []
    sample_ids = []
    census_genes = open("/gpfs/commons/groups/clinical/vshah/pipelines/vep_annotate/cosmic_actionable_old/cancer_gene_census.genes.txt", "r")
    census = census_genes.readlines()
    genes_new = []
    for genes in census:
        gene = genes.strip("\n")
        genes_new.append(gene)
    for file_dragen in dragen:
        file1 = file_dragen.split(".")[0]
        file_new = file1.split("T")[0]
        sample_ids.append(file_new)
    for file_v5b in v5b:
        counter = 0
        file1 = file_v5b.split(".")[0]
        file_new = file1.split("T")[0]
        if file_new in sample_ids:
            print(file_new)
            input_v5b = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/v5b_outputs/" + file_v5b
            input_dragen = "/gpfs/commons/groups/clinical/vshah/somatic_cnv/output_csv_with_census/" + file_new + ".csv"
            try:
                df_v5b = pd.read_csv(input_v5b, sep="\t")
                df_dragen = pd.read_csv(input_dragen, sep=",")
                counter = main_function(input_v5b, input_dragen, file_new, genes_new)
                print(counter)
                perc.append((counter/len(df_v5b))*100)
            except:
                pass
    
            comp.write("\t".join([file_new, str(len(df_dragen)), str(len(df_v5b)), str((counter/len(df_v5b))*100), str(counter)]) + '\n')
    stop = time.time()
    print(stop-start)
    #print(sum(perc)/len(perc))
