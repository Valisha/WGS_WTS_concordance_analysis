#! /nfs/sw/python/python-3.6.1/bin/python3
import sys
import pandas as pd
import argparse
import json

def get_df_mapping_metrics(df_new, tumor_id, normal_id):
    df = pd.read_csv(df_new, sep=",", header=None)
    df.columns = ['Section', 'Sample', 'Metric','Count/Ration/Time', 'Percentage/Seconds']
    df_summary_index_tumor = df.index[df['Section'] == "TUMOR MAPPING/ALIGNING SUMMARY"]
    df_summary_index_normal = df.index[df['Section'] == "NORMAL MAPPING/ALIGNING SUMMARY"]
    df_table_summary_tumor = df.iloc[df_summary_index_tumor, 2:]
    df_table_summary_normal = df.iloc[df_summary_index_normal, 2:]
    columns_of_interest = ['Total input reads', 'Number of duplicate marked reads', 'Q30 bases', 'Estimated read length', 'Insert length: mean', 'Insert length: median', 'Provided sex chromosome ploidy', 'Estimated sample contamination', 'Average sequenced coverage over genome', 'Reads with MAPQ [40:inf)', 'Reads with MAPQ [30:40)']
    df_normal_interest_count = {}
    df_normal_interest_perc = {}
    for values in columns_of_interest:
        if values in ['Reads with MAPQ [40:inf)', 'Reads with MAPQ [30:40)', 'Number of duplicate marked reads', 'Q30 bases']:
            line_new_count = list(df_table_summary_normal['Percentage/Seconds'][df_table_summary_normal['Metric'] == values])[0]
            line_new_perc = list(df_table_summary_normal['Percentage/Seconds'][df_table_summary_normal['Metric'] == values])[0]
            df_normal_interest_count[values] = line_new_count
            str_perc = str("perc" + values)
            df_normal_interest_perc[str_perc] = line_new_perc
        else:
            line_new_count = list(df_table_summary_normal['Count/Ration/Time'][df_table_summary_normal['Metric'] == values])[0]
            line_new_perc = list(df_table_summary_normal['Percentage/Seconds'][df_table_summary_normal['Metric'] == values])[0]
            df_normal_interest_count[values] = line_new_count
            str_perc = str("perc" + values)
            df_normal_interest_perc[str_perc] = line_new_perc
    df_tumor_interest_count = {}
    df_tumor_interest_perc = {}
    for values in columns_of_interest:
        if values in ['Reads with MAPQ [40:inf)', 'Reads with MAPQ [30:40)', 'Number of duplicate marked reads', 'Q30 bases']:
            line_new_count = list(df_table_summary_tumor['Percentage/Seconds'][df_table_summary_tumor['Metric'] == values])[0]
            line_new_perc = list(df_table_summary_tumor['Percentage/Seconds'][df_table_summary_tumor['Metric'] == values])[0]
            str_perc = str("perc" + values)
            df_tumor_interest_perc[str_perc] = line_new_perc
            df_tumor_interest_count[values] = line_new_count
        else:
            line_new_count = list(df_table_summary_tumor['Count/Ration/Time'][df_table_summary_tumor['Metric'] == values])[0]
            line_new_perc = list(df_table_summary_tumor['Percentage/Seconds'][df_table_summary_tumor['Metric'] == values])[0]
            str_perc = str("perc" + values)
            df_tumor_interest_perc[str_perc] = line_new_perc
            df_tumor_interest_count[values] = line_new_count
    normal_id = str(normal_id)
    tumor_id = str(tumor_id)
    df_normal_count = pd.DataFrame(df_normal_interest_count, index=[normal_id])
    df_normal_perc = pd.DataFrame(df_normal_interest_perc, index=[normal_id])
    df_tumor_perc = pd.DataFrame(df_tumor_interest_perc, index=[tumor_id])
    df_tumor_count = pd.DataFrame(df_tumor_interest_count, index=[tumor_id])
    df_tumor_final_perc = df_tumor_perc.append(df_tumor_perc)
    df_tumor_final_count = df_tumor_count.append(df_normal_count)
    #df_final = pd.concat([df_tumor_final_perc, df_tumor_final_count], axis=1)
    #df_final.dropna()
    return df_tumor_final_count

def get_wgs_cov_df(df1, df2, tumor_id, normal_id):
    df_tumor = pd.read_csv(df1, sep=",", header=None, error_bad_lines=False)
    df_normal = pd.read_csv(df2, sep=",", header=None, error_bad_lines=False)
    df_tumor = df_tumor.iloc[:,2:]
    df_tumor.columns = ['Metrics', 'Counts']
    df_tumor_dict = {}
    for i, v in enumerate(df_tumor['Metrics']):
        if v == 'Average alignment coverage over genome' or v.startswith('PCT'):
            df_tumor_dict[v] = list(df_tumor['Counts'][df_tumor['Metrics'] == v])[0]
    df_tumor = pd.DataFrame(df_tumor_dict, index=[tumor_id])
    tumor_id = str(tumor_id)
    normal_id = str(normal_id)
    df_normal = df_normal.iloc[:,2:]
    df_normal.columns = ['Metrics', 'Counts']
    df_normal_dict = {}
    for i, v in enumerate(df_normal['Metrics']):
        if v == 'Average alignment coverage over genome' or v.startswith('PCT') or v == 'Reads with MAPQ [40:inf)' or v == 'Reads with MAPQ [30:40)':
            df_normal_dict[v] = list(df_normal['Counts'][df_normal['Metrics'] == v])[0]
    df_normal = pd.DataFrame(df_normal_dict, index=[normal_id])
    df_wgs_cov = df_tumor.append(df_normal)
    return df_wgs_cov

def get_cnv_metrics(df_cnv, tumor_id, normal_id):
    cnv = pd.read_csv(df_cnv, sep=',', header=None, error_bad_lines=False)
    df = cnv.iloc[:, 2:]
    df.columns = ['Metrics', 'Counts', 'Percentage']
    df_cnv_dict = {}
    df_cnv_dict['Overall ploidy'] = list(df['Counts'][df['Metrics'] == 'Overall ploidy'])[0]
    df_cnv_dict['Estimated tumor purity'] = list(df['Counts'][df['Metrics'] == 'Estimated tumor purity'])[0]
    df_cnv_dict['Diploid coverage'] = list(df['Counts'][df['Metrics'] == 'Diploid coverage'])[0]
    ids = str(tumor_id + "--" + normal_id)
    df_cnv_new = pd.DataFrame(df_cnv_dict, index=[ids])
    df_cnv_new = df_cnv_new.mask(df_cnv_new.eq('NaN')).dropna()
    print(df_cnv_new)
    return df_cnv_new

def get_at_gc_dropout(gc_normal, gc_tumor, tumor_id, normal_id):
    gcn = pd.read_csv(gc_normal, sep=',')
    gct = pd.read_csv(gc_tumor, sep=',')
    dft = gct.iloc[:, 2:]
    dfn = gcn.iloc[:, 2:]
    dft.columns = ['Metrics', 'Counts', 'Perc']
    dfn.columns = ['Metrics', 'Counts', 'Perc']
    df_gc_tumor = {}
    df_gc_normal = {}
    df_gc_tumor['AT_Dropout'] = list(dft['Counts'][dft['Metrics'] == 'AT Dropout'])[0]
    df_gc_tumor['GC_Dropout'] = list(dft['Counts'][dft['Metrics'] == 'GC Dropout'])[0]
    df_gc_normal['AT_Dropout'] = list(dfn['Counts'][dfn['Metrics'] == 'AT Dropout'])[0]
    df_gc_normal['GC_Dropout'] = list(dfn['Counts'][dfn['Metrics'] == 'GC Dropout'])[0]
    normal_id = str(normal_id)
    tumor_id = str(tumor_id)
    gc_normal_df = pd.DataFrame(df_gc_normal, index=[normal_id])
    gc_tumor_df = pd.DataFrame(df_gc_tumor, index=[tumor_id])
    gc_df = gc_tumor_df.append(gc_normal_df)
    return gc_df

def time_metrics(df_tumor_time, tumor_id, normal_id):
    time_tumor = pd.read_csv(df_tumor_time, sep=',')
    time = time_tumor.iloc[(len(time_tumor)-1),3]
    time_dict = {}
    time_dict['time_metrics_tumor_normal'] = time
    print(time)
    ids = str(tumor_id) + '--' + str(normal_id)
    time_tumor_df = pd.DataFrame(time_dict, index=[ids]) 
    return time_tumor_df
     
def merge_dataframes(df_tumor_final, df_wgs_cov, df_cnv_new, gc_df, time_tumor_df):
    df1 = pd.concat([df_tumor_final, df_wgs_cov, df_cnv_new, gc_df, time_tumor_df], axis=1)
    return df1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process QC metrics from dragen')
    parser.add_argument('--mapping_metrics', '-m',type=str, dest='df_mapping_metrics', help='mapping_metrics')
    parser.add_argument('--wgs_cov_tumor', '-wt', type=str, dest='df_wgs_cov_tumor', help='WGS tumor')
    parser.add_argument('--wgs_cov_normal', '-wn', type=str, dest='df_wgs_cov_normal', help='WGS normal')
    parser.add_argument('--output_tsv', '-ot', type=str, dest='output_tsv', help='output_tsv')
    parser.add_argument('--cnv_metrics', '-cnv', type=str, dest='df_cnv', help='cnv_metrics')
    parser.add_argument('--gc_normal', '-gn', type=str, dest='gc_normal', help='atgc_normal')
    parser.add_argument('--gc_tumor', '-gt', type=str, dest='gc_tumor', help='atgc_tumor')
    parser.add_argument('--time_metrics', '-tm', type=str, dest='time_tumor', help='time_metrics_file')
    args = parser.parse_args()
    path = str(args.df_mapping_metrics).split("/")[-1]
    sample_id = path.split("T")[0]
    tumor_id = str(sample_id + "T-1-D1")
    normal_id = str(sample_id + "N-1-D1")
    df_summary = get_df_mapping_metrics(args.df_mapping_metrics, tumor_id, normal_id)
    df_wgs_cov = get_wgs_cov_df(args.df_wgs_cov_tumor, args.df_wgs_cov_normal, tumor_id, normal_id)
    df_cnv_new = get_cnv_metrics(args.df_cnv, tumor_id, normal_id)
    df_gc = get_at_gc_dropout(args.gc_normal, args.gc_tumor, tumor_id, normal_id)
    time_tumor = time_metrics(args.time_tumor, tumor_id, normal_id)
    df_final = merge_dataframes(df_summary, df_wgs_cov, df_cnv_new, df_gc, time_tumor)
    df_final.to_csv(args.output_tsv, sep='\t')

