#! /nfs/sw/python/python-3.6.1/bin/python3


def parse_vcf(vcf, sample_id):
    vcf_new = open(vcf ,'r')
    vcf_lines = vcf_new.readlines()
    index = 0
    vcf_list =[]
    for line in vcf_lines:
        if '#' not in line:
            line = line.strip('\n').split('\t')
            vcf_list.append(line)
    vcf_df = pd.DataFrame(vcf_list)
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_id]
    return vcf_df

def parse_format_column(vcf_df, sample_id):
    vcf_df[['GT', 'CN', 'MCN', 'CNQ', 'MCNQ', 'CNF', 'MCNF', 'SD', 'MAF', 'BC', 'AS', 'PE']] = vcf_df[sample_id].str.split(pat=":", expand=True)
    vcf_df.drop(['QUAL', 'INFO', 'FORMAT', sample_id], inplace=True, axis=1)
    return vcf_df

def get_gene_info(cnv_annotated, vcf_df, sample_id, dragen_purity, dragen_ploidy, cov_tumor, cov_normal, contam_tumor, contam_normal):
    out_tsv = open(str('merged_tsvs/' + sample_id + '.merged_annotated.tsv'), 'w')
    variants = ['FLT4', 'MYC', 'CDKN2A', 'PTEN','FGFR2', 'CCND2', 'RB1', 'BRCA2','IDH2', 'ERBB2', 'CCNE1', 'EGFR', 'JAK1', 'CIC', 'IDH1']
    cnv_annotated_table = pd.read_csv(cnv_annotated, sep='\t')
    cnv_annotated_table.columns = ['chr', 'start', 'end', 'Gene', 'chr_1', 'start_1', 'end_1', 'CN', 'overlap']
    cnv_table = cnv_annotated_table.loc[:,  ['Gene', 'chr_1', 'start_1', 'end_1']]
    out_tsv.write('\t'.join(['Sample_ID', 'Gene', str('\t'.join(vcf_df.columns)), 'Purity', 'Ploidy', 'Coverage_tumor', 'Coverage_normal', 'Contam_tumor', 'Contam_normal']) + '\n')
    for variant in variants:
        for i, v in enumerate(vcf_df['ID']):
            start_end = v.split(":")[3]
            start = start_end.split("-")[0]
            end = start_end.split("-")[1]
            for i_a, v_a in enumerate(cnv_table['start_1']):
                if int(start) == int(v_a) and str(cnv_table.loc[i_a, 'Gene']) == str(variant):
                    purity = dragen_purity[sample_id]
                    ploidy = dragen_ploidy[sample_id]
                    coverage_tumor = cov_tumor[sample_id]
                    coverage_normal = cov_normal[sample_id]
                    sample_contam_tumor = contam_tumor[sample_id]
                    sample_contam_normal = contam_normal[sample_id]
                    out_tsv.write(sample_id + '\t' + variant + '\t' + str('\t'.join([str(f) for f in vcf_df.loc[i, :].values])) + '\t' + str(purity) + '\t' + str(ploidy) + '\t' + str(coverage_tumor) + '\t' + str(coverage_normal) + '\t' + str(sample_contam_tumor) + '\t' + str(sample_contam_normal) + '\n')
    out_tsv.close()

def exists(var):
    return var in globals()

def check_concordance():
    df = pd.read_csv('project_10800_LOD_purity_ploidy.tsv', sep='\t')
    df['Concordance'] = pd.Series()
    samples = list(df.loc[:, 'Sample_ID'].values)
    #print(samples)
    sample_ids = []
    variants = ['FLT4', 'MYC', 'CDKN2A', 'PTEN','FGFR2', 'CCND2', 'RB1', 'BRCA2','IDH2', 'ERBB2', 'CCNE1', 'EGFR', 'JAK1', 'CIC', 'IDH1']
    for value in samples:
        value_new = value.split('T')[0]
        if value_new not in sample_ids:
            sample_ids.append(value_new)
    for value in sample_ids:
        for variant in variants:
            var_df = new_df[new_df['Gene'] == variant]
            truth_set = var_df[~var_df['Sample_ID'].str.endswith('percent')]
            try:
                truth_CN = int(truth_set['CN'].values[0])
            except IndexError:
                pass
            percent_set = var_df[var_df['Sample_ID'].str.endswith('percent')]
            for index, v in enumerate(percent_set.loc[:, 'Sample_ID']):
                i = df.index[(df['Sample_ID'] == v) & (df['Gene'] == variant)].values[0]
                #print(i)
                cns = int(percent_set['CN'][percent_set['Sample_ID'] == v])
                if int(truth_CN) == 0 and int(cns) == 0:
                    print(truth_CN, cns, percent_set['Sample_ID'][percent_set['Sample_ID'] == v]) 
                    df.loc[i, 'Concordance'] = 'Concordant'
                    #print(v, '\t'.join([str(f) for f in df.loc[i, :].values]))
                    #print(df.loc[i, :].values)
                elif int(truth_CN) == 1 and int(cns) == 1:
                    df.loc[i, 'Concordance'] = 'Concordant'
                    #print(v + '\t' + '\t'.join([str(f) for f in df.loc[i, :].values]) + '\n')
                    #print(df.loc[i, :].values)
                elif int(truth_CN) == 2 and int(cns) == 2:
                    df.loc[i, 'Concordance'] = 'Concordant'
                    #print(v, '\t'.join([str(f) for f in df.loc[i, :].values]))
                    #print(df.loc[i, :].values)
                elif int(truth_CN) >= 3 and int(cns) >= 3:
                    df.loc[i, 'Concordance'] = 'Concordant'
                    #print(v, '\t'.join([str(f) for f in df.loc[i, :].values]))
                    #print(df.loc[i, :].values)
                else:
                    df.loc[i, 'Concordance'] = 'Discordant'
                    #print(v, '\t'.join([str(f) for f in df.loc[i, :].values]))
                    #print(df.loc[i, :].values)
    #print(df)
    return df


def run_bedtools(sample_id):
    call_string = str("bedtools intersect -b " + str(sample_id+'.bed') + " -a outcc.bed -wo > annotated_bedtools_tsvs/" + str(sample_id + "_annotated_bedtools.tsv"))
    subprocess.call(call_string, shell=True)

if __name__ == "__main__":
    import os
    import sys
    import subprocess
    import pandas as pd
    import adjusted_log2_dragen_conc_test as adj
    cnv_vcf_list = open(sys.argv[1], 'r')
    mapping = open(sys.argv[2], 'r')
    cnv_vcf_lines = cnv_vcf_list.readlines()
    dragen_purity_ploidy = open(sys.argv[3], "r")
    dragen_purity_ploidy_lines = dragen_purity_ploidy.readlines()
    dragen_purity = {}
    dragen_ploidy = {}
    for line in dragen_purity_ploidy_lines:
        line_new = line.strip('\n').split('/')[-1]
        sample = line_new.split('--')[0]
        print(sample)
        if sample not in dragen_purity.keys():
            cnv = pd.read_csv(str(line.strip('\n')), sep=',', header=None, error_bad_lines=False)
            df = cnv.iloc[:, 2:]
            df.columns = ['Metrics', 'Counts', 'Percentage']
            dragen_ploidy.update({sample:list(df['Counts'][df['Metrics'] == 'Overall ploidy'])[0]})
            dragen_purity.update({sample:list(df['Counts'][df['Metrics'] == 'Estimated tumor purity'])[0]})
    average_sequenced_cov_normal = {}
    average_sequenced_cov_tumor = {}
    sample_contam_normal = {}
    sample_contam_tumor = {}
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
    for cnv_vcf in cnv_vcf_lines:
        vcf = cnv_vcf.strip('\n')
        sample_id = str((vcf.split('/')[-1]).split('--')[0])
        vcf_df = parse_vcf(vcf, sample_id)
        parsed_vcf_df = parse_format_column(vcf_df, sample_id)
        adj.parse_vcf_to_bed(vcf, sample_id)
        run_bedtools(sample_id)
        get_gene_info(str('annotated_bedtools_tsvs/' + sample_id + '_annotated_bedtools.tsv'), parsed_vcf_df, sample_id, dragen_purity, dragen_ploidy, average_sequenced_cov_tumor, average_sequenced_cov_normal, sample_contam_tumor, sample_contam_tumor)
    proj = os.listdir('merged_tsvs/')
    merged_df = pd.DataFrame()
    for lines in proj:
        if 'merged' in lines:
            df = pd.read_csv(str('merged_tsvs/' + lines), sep='\t')
            merged_df = merged_df.append(df)
    #print(merged_df)
    merged_df.to_csv('ONC21-9_LOD_purity_ploidy.tsv', sep='\t')
    final_df = check_concordance()
    final_df.to_csv('ONC21-9_LOD-final.tsv', sep='\t')
