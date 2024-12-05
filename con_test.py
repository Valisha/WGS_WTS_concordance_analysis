#! /nfs/sw/python/python-3.6.1/bin/python3
#SBATCH --job-name=concordance_test
#SBATCH --export=ALL
#SBATCH --mem=50G
#SBATCH -n10
#SBATCH --output=changeme.%j.out
#SBATCH --account=clinical

import pandas as pd 
import re 
import sys
from itertools import groupby

def parse_sample(sample_df):
    # splits alteration column into two HGVSp and HGVSc
    new = sample_df['Alteration'].str.split(" ", expand=True)
    sample_df['HGVSc'] = new[1]
    sample_df['HGVSp'] = new[0]
    for i, v in enumerate(sample_df['HGVSp']):
        if 'p.' in v:
            sample_df.loc[i, 'HGVSp'] = str(v[2:])
    indexes = []
    if 'Additional Information' in sample_df.columns:
         for i, v in enumerate(sample_df['Additional Information']):
             if 'MAF' in str(v):
                 indexes.append(i)
         sample_df = sample_df.iloc[indexes, :]
    sample_df.drop('Alteration', inplace=True, axis=1)
    return sample_df

def parse_res(df_res):
    #['CHR', 'POS', 'REF', 'ALT', 'FILTER', 'Type','IMPACT', 'Gene', 'HGVSc', 'HGVSp', 'VAF', 'SQ']
    genes = []
    hgvsp = []
    hgvsc = []
    chrs = []
    pos = []
    ref = []
    alt = []
    filters = []
    impacts = []
    vaf = []
    sq = []
    types = []
    df_res_parsed = pd.DataFrame()
    for i, var in enumerate(df_res['HGVSp']):
        if str(var) != 'nan':
           p1 = var.split("p.") 
           if '=' not in str(p1[1]):
               hgvsp.append(str(p1[1]))
               genes.append(str(df_res.loc[i, 'Gene']))
               hgvsc.append(str(df_res.loc[i, 'HGVSc']))
               chrs.append(str(df_res.loc[i, 'CHR']))
               pos.append(str(df_res.loc[i, 'POS']))
               ref.append(str(df_res.loc[i, 'REF']))
               alt.append(str(df_res.loc[i, 'ALT']))
               filters.append(str(df_res.loc[i, 'FILTER']))
               impacts.append(str(df_res.loc[i, 'IMPACT']))
               vaf.append(str(df_res.loc[i, 'VAF']))
               sq.append(str(df_res.loc[i, 'SQ']))
               types.append(str(df_res.loc[i, 'Type']))
    df_res_parsed['CHR'] = chrs
    df_res_parsed['POS'] = pos
    df_res_parsed['REF'] = ref
    df_res_parsed['ALT'] = alt
    df_res_parsed['FILTER'] = filters
    df_res_parsed['Type'] = types
    df_res_parsed['IMPACT'] = impacts
    df_res_parsed['Gene'] = genes
    df_res_parsed['HGVSc'] = hgvsc
    df_res_parsed['HGVSp'] = hgvsp
    df_res_parsed['VAF'] = vaf
    df_res_parsed['SQ'] = sq
    return df_res_parsed
def analyse_variants_not_found(df_mskcc, num, df_res, input_vep_vcf):
    df_aa = {'Ala': 'A','Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
    val_list = list(df_aa.values())
    key_list = list(df_aa.keys())
    out_cordant = open(str('out_cordant_'+num+'.tsv'), 'w')
    out_noncordant = open(str('out_non_cordant_'+num+'.tsv'), 'w')
    new_variants = []
    res_variants = list(df_res.loc[:, 'HGVSp'])
    i = ''
    line_number_non_cordant = 0
    found_variants = []
    list_of_results = []
    missing = []
    line_number = 0
    not_found = 0
    old_variants = list(df_mskcc['HGVSp'])
    for i_mskcc, var in enumerate(df_mskcc['HGVSp']):
        res =re.findall(r'[A-Za-z]+|\d+|[@_!#$%^&*()<>?/\|}{~:]', var)
        aa_string = ''
        for i, v in enumerate(res):
            try: 
                 if bool(re.search(r'\d', v)) == False:
                    if len(v) == 1:
                        position = val_list.index(v)
                        var_new = key_list[position]
                        aa_string += var_new
                    new_val = ''
                    if len(v) > 1:
                        for value in v:
                            if value in val_list:
                                position = val_list.index(value)
                                var_new = key_list[position]
                                aa_string += var_new
                            if value not in val_list:
                                new_val += value
                    aa_string += new_val
                 if bool(re.search(r'\d', v)) == True:
                    aa_string += v
            except ValueError:
                print(var, "Not a variant")
        
        new_variants.append(aa_string)
        
        if aa_string in res_variants:
            line_number += 1
            print(aa_string)
            row_new = list(df_res[df_res['HGVSp'] == aa_string])
            row_no_nan = [i for i in row_new if i]
            out_cordant.write('\t'.join([str(line_number), str('\t'.join(row_no_nan)), str('\t'.join(list(df_mskcc.loc[i_mskcc, :])))]) + "\n")
            found_variants.append(aa_string)
    for nv in list(new_variants):
        if nv not in found_variants:
            missing.append(nv)
            line_number_non_cordant += 1
            index_nv = new_variants.index(nv)
            try:
                row_missing_mskcc = list(df_mskcc.loc[index_nv, :])
                row_missing_no_nan = [i for i in row_missing_mskcc if i]
            except KeyError:
                print('index not found')
            try:
                out_noncordant.write(str(line_number_non_cordant) + "\t" + str('\t'.join(row_missing_no_nan)) + "\n")
            except TypeError:
                new_list = []
                for value in row_missing_no_nan:
                    if str(value) != 'nan':
                        new_list.append(value)
                out_noncordant.write(str(line_number_non_cordant) + "\t" + str("\t".join(new_list)))
    with open(input_vep_vcf, 'r') as vcf:
        for missing_vars in missing:
            index = (df_mskcc.index[df_mskcc['HGVSp'] == missing_vars].tolist())
            res =re.findall(r'[A-Za-z]+|\d+|[@_!#$%^&*()<>?/\|}{~:]', missing_vars)
            if len(res) > 1:
                aa = ''.join([res[0], res[1]])
                for i_new in index:
                    for line in vcf:
                        if aa in  line and str(df_mskcc.loc[i_new, 'Gene']) in line:
                            print('matched aa: ', aa)

    print('missing:', missing) 
    print('length missing: ', len(missing))
    print('length found variants: ', len(found_variants))
    print('length mskcc: ', len(list(df_mskcc['HGVSp'])))
    return missing, len(missing), len(found_variants), len(new_variants)


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings('ignore')
    import sys
    import pandas as pd
    import re
    import subprocess
    import os
    from collections import Counter
    import time
    import numpy as np
    start = time.time()
    numbers_res = []
    filenames_res_outputs = os.listdir('/gpfs/commons/groups/clinical/vshah/concordance_test/no_multi')
    filenames_mskcc = os.listdir('/gpfs/commons/groups/clinical/vshah/concordance_test/final_outputs/samples_mskcc/')
    new_df = open(str('new_df.tsv'), 'w')
    new_df.write("Sample" + "\t" + "MSKCC Count" + "\t" + "Dragen Count" + "\t" + "Perc concordance" + "\t" + "Non cordant" + "\t" + "IMPACT"+ "\n")
    for file1 in filenames_res_outputs:
        split1 = file1.split("-")
        numbers_res.append(split1[1].split("T")[0])
    for file2 in filenames_mskcc:
        split2 = file2.split("_")
        num = split2[1]
        if num in numbers_res:
            sample_id = 'RES20-' + num + 'T-1-D1'
            print(sample_id)
            input_mskcc = '/gpfs/commons/groups/clinical/vshah/concordance_test/final_outputs/samples_mskcc/' + str('sample_'+num+'_mskcc.tsv')
            input_vep_vcf = '/gpfs/commons/groups/clinical/vshah/pipelines/vep_annotate/res20_annotation/no_multi/' + 'RES20-' + num + 'T-1-D1--RES20-' + num + 'N-1-D1.hard-filtered.annotated_filter_vep_table.txt'
            vcf = '/gpfs/commons/groups/clinical/vshah/pipelines/vep_annotate/res20_annotation/no_multi/' + 'RES20-' + num + 'T-1-D1--RES20-' + num + 'N-1-D1.hard-filtered.annotated_vep.vcf'
            try:
                df_mskcc = pd.read_csv(input_mskcc, sep="\t")
                df_res = pd.read_csv(input_vep_vcf, sep="\t")
                df_res = df_res.loc[:, ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'Consequence', 'IMPACT', 'SYMBOL', 'HGVSc', 'HGVSp', str(sample_id + '.AF'), str(sample_id + '.SQ')]]
                df_res.columns = ['CHR', 'POS', 'REF', 'ALT', 'FILTER', 'Type','IMPACT', 'Gene', 'HGVSc', 'HGVSp', 'VAF', 'SQ']
                df_mskcc_parsed = parse_sample(df_mskcc)
                df_res_parsed =  parse_res(df_res)
                missing, len_missing, len_found_variants, len_new_variants = analyse_variants_not_found(df_mskcc_parsed, num, df_res_parsed, vcf)
                new_df.write(str(sample_id) + "\t" + str(len_new_variants) + "\t" + str(len_found_variants) + "\t" + str(len_found_variants/len_new_variants) + "\t" + str(len_missing) + "\n")
            except:
                new_df.write(sample_id + "\t" + "Data Frame empty" + "\n")
    new_df.close()
    end = time.time()
    print('time taken to run', end-start)
