#! /nfs/sw/python/python-3.6.1/bin/python3

import pandas as pd
def make_truth_set(truth, sample_id):
    truth = pd.read_csv(truth, sep='\t')
    truth = truth.loc[:, ['SYMBOL', 'CHROM', 'POS', 'REF', 'ALT', 'FILTER','IMPACT', 'HGVSp', 'HGVSc', str(sample_id + '.AD'), str(sample_id + '.AF')]]
    parsed_truth_set = []
    truth_new = []
    genes_set = ['IDH1', 'CIC', 'MMP2', 'JAK1', 'EGFR', 'KIT', 'PAX7', 'CBLC', 'MYO5A']
    impacts = ['HIGH', 'MODERATE', 'LOW']
    for i, v in enumerate(truth['SYMBOL']):
        if type(v) != float and ',' in v:
            variants = v.split(',')
            for i_var, v_var in enumerate(variants):
                if type(v_var) != float and str(v_var) != 'NA' and str(v_var) != 'nan' and 'LOC' not in str(v_var) and str(v_var) in genes_set and str(truth.loc[i, 'IMPACT']).split(',')[i_var] in impacts:
                    truth.loc[i, 'SYMBOL'] = v_var
                    truth_new.append(truth.loc[i, :].values)
        if type(v) != float and ',' not in v:
            if type(v) != float and str(v) != 'NA' and str(v) != 'nan' and 'LOC' not in str(v) and str(v) in genes_set and str(truth.loc[i, 'IMPACT']) in impacts:
                truth_new.append(truth.loc[i, :].values)
    truth_df = pd.DataFrame(truth_new)
    truth_df.columns = truth.columns
    return truth_df

def compare(percent_df, truth_df, sample_id, sample_id_perc, average_sequenced_cov_normal, average_sequenced_cov_tumor):
    percent_df = pd.read_csv(percent_df, sep='\t')
    genes_found = []
    genes_set = ['IDH1', 'CIC', 'MMP2', 'JAK1', 'EGFR', 'KIT', 'PAX7', 'CBLC', 'MYO5A']
    percent_df = percent_df.loc[:, ['SYMBOL', 'CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'IMPACT', 'HGVSp', 'HGVSc', str(sample_id_perc + '.AD'), str(sample_id_perc + '.AF')]]
    for i, v in enumerate(truth_df['POS']):
        for i_perc, v_perc in enumerate(percent_df['POS']):
             if str(v_perc) == str(v):
                 genes_found.append(v_perc)
                 read = truth_df.loc[i, str(sample_id + '.AD')].split(',')
                 read_count_truth = float(read[0]) + float(read[1])
                 read_perc = percent_df.loc[i_perc, str(sample_id_perc + '.AD')].split(',')
                 read_count_perc = float(read_perc[0]) + float(read_perc[1])
                 if len(truth_df.loc[i, 'REF']) == 1 and len(truth_df.loc[i, 'ALT']) == 1:
                     out_file.write('\t'.join([str(sample_id), '\t'.join([str(f) for f in truth_df.loc[i, :].values]), str(read_count_truth)]) + '\n') 
                     out_file.write('\t'.join([str(sample_id_perc), '\t'.join([str(f) for f in percent_df.loc[i_perc, :]]), str(read_count_perc), 'SNP', 'CONCORDANT']) + '\n')
                 elif len(truth_df.loc[i, 'REF']) > 1 or len(truth_df.loc[i, 'ALT']) > 1:
                     out_file.write('\t'.join([str(sample_id), '\t'.join([str(f) for f in truth_df.loc[i, :].values]), str(read_count_truth)]) + '\n')
                     out_file.write('\t'.join([str(sample_id_perc), '\t'.join([str(f) for f in percent_df.loc[i_perc, :]]), str(read_count_perc), 'INDEL', 'CONCORDANT']) + '\n')
    for genes in genes_found:
        if genes not in truth_df['POS'].values:
            print(truth_df['SYMBOL'][truth_df['POS'] == genes])
            out_file.write('\t'.join([str(sample_id), str('\t'.join([str(f) for f in (truth_df[truth_df['POS'] == genes].values)])), 'Gene Not Found in Dilution Sample']) + '\n')

out_file = open('final_LOD_SNV_CA_G15.txt', 'w')
average_sequenced_cov_normal = {}
average_sequenced_cov_tumor = {}
sample_contam_normal = {}
sample_contam_tumor = {}
mapping = open("/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/mapping_metrics_csv.txt", "r")
mapping_list = mapping.readlines()
for line in mapping_list:
    line_new = line.strip('\n').split('/')[-1]
    sample = line_new.split('--')[0]
    if sample not in average_sequenced_cov_tumor.keys():
        mapping = pd.read_csv(str(line.strip('\n')), sep=',', header=None)
        mapping.columns = ['Section', 'Sample', 'Metric','Count/Ration/Time', 'Percentage/Seconds']
        df_summary_index_tumor = mapping.index[mapping['Section'] == "TUMOR MAPPING/ALIGNING SUMMARY"]
        df_summary_index_normal = mapping.index[mapping['Section'] == "NORMAL MAPPING/ALIGNING SUMMARY"]
        df_table_summary_tumor = mapping.iloc[df_summary_index_tumor, 2:]
        df_table_summary_normal = mapping.iloc[df_summary_index_normal, 2:]
        columns_of_interest = ['Estimated sample contamination', 'Average sequenced coverage over genome']
        sample_contam_normal.update({sample:list(df_table_summary_normal['Count/Ration/Time'][df_table_summary_normal['Metric'] == 'Estimated sample contamination'])[0]})
        sample_contam_tumor.update({sample:list(df_table_summary_tumor['Count/Ration/Time'][df_table_summary_tumor['Metric'] == 'Estimated sample contamination'])[0]})
        average_sequenced_cov_normal.update({sample:list(df_table_summary_normal['Count/Ration/Time'][df_table_summary_normal['Metric'] == 'Average sequenced coverage over genome'])[0]})
        average_sequenced_cov_tumor.update({sample:list(df_table_summary_tumor['Count/Ration/Time'][df_table_summary_tumor['Metric'] == 'Average sequenced coverage over genome'])[0]})
out_file.write('\t'.join(['Sample_ID', 'SYMBOL', 'CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'IMPACT', 'Consequence', 'HGVSp', 'HGVSc', 'AD', 'AF', 'ReadCount', 'SNP/INDEL']) + '\n')
g15_parsed = make_truth_set('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/truth_set/G15-35T-D--G15-35N-D.hard-filtered.annotated_filter_vep_aq-table.txt', 'G15-35T-D')
ca_parsed = make_truth_set('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/truth_set/CA-0061T-D-W--CA-0061N-D-W.hard-filtered.annotated_filter_vep_aq-table.txt', 'CA-0061T-D-W')
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/CA-0061T-D-10percent--CA-0061N-D-W.hard-filtered.annotated_filter_vep_aq-table.txt', ca_parsed, 'CA-0061T-D-W', 'CA-0061T-D-10percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/CA-0061T-D-20percent--CA-0061N-D-W.hard-filtered.annotated_filter_vep_aq-table.txt', ca_parsed, 'CA-0061T-D-W', 'CA-0061T-D-20percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/CA-0061T-D-30percent--CA-0061N-D-W.hard-filtered.annotated_filter_vep_aq-table.txt', ca_parsed, 'CA-0061T-D-W', 'CA-0061T-D-30percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/G15-35T-D-10percent--G15-35N-D.hard-filtered.annotated_filter_vep_aq-table.txt', g15_parsed, 'G15-35T-D', 'G15-35T-D-10percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/G15-35T-D-30percent--G15-35N-D.hard-filtered.annotated_filter_vep_aq-table.txt', g15_parsed, 'G15-35T-D', 'G15-35T-D-30percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
compare('/gpfs/commons/groups/clinical/vshah/limit_of_detection_study/LOD-SNV/percent_set/G15-35T-D-50percent--G15-35N-D.hard-filtered.annotated_filter_vep_aq-table.txt', g15_parsed, 'G15-35T-D', 'G15-35T-D-50percent', average_sequenced_cov_normal, average_sequenced_cov_tumor)
