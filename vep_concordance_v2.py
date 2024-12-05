#! /nfs/sw/python/python-3.6.1/bin/python3
#SBATCH --job-name=concordance_test
#SBATCH --export=ALL
#SBATCH --mem=50G
#SBATCH -n10
#SBATCH --output=changeme.%j.out
#SBATCH --account=clinical
# sed -n '/#/!p' res_output/RES20-03T-1-D1--RES20-03N-1-D1.hard-filtered.annotated_filter_vep_table.txt | cut -f 7,9,16,17,66 | less > test_RES20-03.txt

def convert(v2):
    """
    convert three letter aa to one letter
    """
    df_aa = {'Ala': 'A','Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
    if v2 in df_aa.keys():
        return df_aa[v2]
def parse_res_output(df_res):
    """
    parses RES output to make it look similar to MSKCC
    """
    chr_list = []
    pos_list = []
    hgvsp_list = []
    hgvsc_list = []
    vaf_list = []
    alt_list = []
    ref_list = []
    impact_list = []
    gene_list = []
    filter_list = []
    type_list = []
    df_aa = {'Ala': 'A','Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': ''}
    for i in range(len(df_res['HGVSc'])):
        v = str(df_res.loc[i, 'HGVSc'])
        p = str(df_res.loc[i, 'HGVSp'])
        if "p." in str(p) and df_res.loc[i, 'Type'] != 'synonymous_variant' and "%" not in p:
            # checks is hgvsc or hgvsp is nan, removes synonymous variant, removes hgvsp value with %3D and takes only first value from HGVSp (if there are two)
            value = str( "c." + str(v.split("c.")[1]))
            # HGVSc corrected
            value_p1 = convert(str(p.split("p.")[1])[0:3])
            v2 = p.split("p.")[1]
            # first three amino acids in HGVSp converted to 1
            length = len(v2)
            # checks for the last three characters 
            v_new = (v2[length - 3: length])
            # changes stop gained in RES nonsense mutation
            if df_res.loc[i, 'Type'] == 'stop gained':
                df_res.loc[i, 'Type'] == 'Nonsense Mutation'
            if v_new in df_aa.keys():
                # it will be matched if there is no fs or ext in the HGVSp
                value_p2 = convert(v_new)
                p1 = str(value_p1 + v2[3: length - 3] + value_p2)
            if v_new not in df_aa.keys():
                # changes HGVSp according fs and ext in between
                res = re.findall(r'(\w+?)(\d+)', v2)
                out = [item for t in res for item in t]
                if str(out[0]) == "Met" and len(out) == 2:
                    value_p1 = convert(out[0])
                    p1 = str(value_p1 + out[1] + "?")
                if len(out) > 2:
                    value_p1 = convert(out[0])
                    ref_number = out[1]
                    value_p2 = convert(str(out[2])[0:3])
                    if value_p2 != None:
                        if len(str(out[2])[3:]) > 3:
                             length_out = len(out[2])
                             last_val = str(out[2])[length_out - 3:]
                             last_val = convert(last_val)
                             fs = str(out[2])[3:length_out - 3]
                             p1 = str(value_p1 + ref_number + value_p2 + fs + last_val + out[3])
                        else:
                             p1 = str(value_p1 + ref_number + value_p2)
            hgvsc_list.append(value)
            hgvsp_list.append(p1)
            chr_list.append(str(df_res.loc[i, 'CHR']))
            filter_list.append(str(df_res.loc[i, 'FILTER']))
            gene_list.append(str(df_res.loc[i, 'Gene']))
            type_list.append(str(df_res.loc[i, 'Type']))
            vaf_list.append(str(df_res.loc[i, 'VAF']))
            pos_list.append(str(df_res.loc[i, 'POS']))
            alt_list.append(str(df_res.loc[i, 'ALT']))
            ref_list.append(str(df_res.loc[i, 'REF']))
            impact_list.append(str(df_res.loc[i, 'IMPACT']))
            p1 = '' ; fs = '' ; last_val = '' ; length_out = '' ; value_p2 = '' ; ref_number = '' ; value_p1 = '' ; p1 = '' ; out = [] ; v_new = '' ; v2 = ''
    dict_res_new = {'Gene': gene_list, 'IMPACT': impact_list, 'Type': type_list, 'VAF': vaf_list, 'HGVSc': hgvsc_list, 'HGVSp': hgvsp_list, 'FILTER': filter_list, 'REF': ref_list, 'ALT': alt_list, 'POS': pos_list, 'CHR': chr_list}
    df_res_new = pd.DataFrame(dict_res_new)
    df_res_new.to_csv("df_res.tsv", sep="\t")
    return df_res_new


def parse_sample(sample_df):
    # splits alteration column into two HGVSp and HGVSc
    new = sample_df['Alteration'].str.split(" ", expand=True)
    sample_df['HGVSc'] = new[1]
    sample_df['HGVSp'] = new[0]
    indexes = []
    for i, v in enumerate(sample_df['Additional Information']):
         if 'MAF' in str(v):
             indexes.append(i)
    sample_df = sample_df.iloc[indexes, :]
    sample_df.drop('Alteration', inplace=True, axis=1)
    return sample_df

def getWord(word):
    return lyst.get(word, word)

def merge_dfs(sample_df, df_res, sample_id):
    print(list(sample_df['HGVSp']))
    print(list(df_res['HGVSp']))
    index_found = []
    na_counter = 0
    countcord = 0
    count_noncordant = 0
    index_a = []
    impact_noncordant_list = []
    part_string = str('partial_concordance' + sample_id + '.tsv')
    part_file = open(part_string, 'w')
    on_string = str('new_df_non_cordant' + sample_id + '.tsv')
    out_file_noncordant = open(on_string, 'w')
    out_string = "new_df_conc" + sample_id + ".tsv"
    out_file = open(out_string, 'w')
    #print("Count" + "\t" + "Gene_RES" + "\t" + "HGVSc_RES" + "\t" + "HGVSp_RES" + "\t" + "Type_RES"+  "\t" + "VAF_RES" +"\t" + "Gene_MSKCC" +"\t" + "Type_MSKCC"+"\t" + "Add_MSKCC" "\t" + "HGVSc_MSKCC" + "\t" + "HGVSp_MSKCC" + "\t" + "Concordance" + "\n")
    out_file.write('Count' + '\t' + '\t'.join(list(df_res.columns)))
    for i, v in enumerate(sample_df['HGVSp']):
        if sample_df.loc[i, 'Gene'] in list(df_res.loc[:, 'Gene']) and sample_df.loc[i, 'HGVSp'] in list(df_res.loc[:, 'HGVSp']) and i not in index_found:
            gene_list = list(np.where(df_res['Gene'] == sample_df.loc[i, 'Gene']))[0]
            df_gene = df_res.iloc[gene_list, :]
            a = list(np.where(df_res['HGVSp'] == sample_df.loc[i, 'HGVSp']))[0]
            countcord += 1
            if len(a) > 1:
                for a_val in a:
                    if df_res.loc[a_val, 'IMPACT'] != 'MODIFIER' and str(df_res.loc[a_val, :].any()) != "NaN" and str(df_res.loc[a_val, :].any()) != None and a_val not in index_a:
                        b = df_res.loc[a_val, 'HGVSp']
                        if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and v == b:
                            countcord += 1
                            index_a.append(a_val)
                            index_found.append(i)
                            # checks if both the variants are the same
                            out_file.write(str(len(index_found)) + '\t' +'\t'.join(sample_df.loc[i,:]) + "\t" + '\t'.join(df_res.loc[a_val,:]) + "\t" + 'Concordant' + "\n")
                            new_line = str(len(index_found)) + '\t' +'\t'.join(sample_df.loc[i,:]) + "\t" + '\t'.join(df_res.loc[a_val,:]) + "\t" + 'Concordant' + "\n"
                        else:
                            if sample_df.loc[i, 'HGVSp'][:-1] == df_res.loc[a_val, 'HGVSp'][:-1] and sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene']:
                                #print('Gene is the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a, :]))
                                part_file.write('Gene is the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                            if sample_df.loc[i, 'Gene'] != df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][:-1] == df_res.loc[a_val, 'HGVSp'][:-1]:
                                #print('Gene is not the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a, :]))
                                part_file.write('Gene is not the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                            if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][:-1] != df_res.loc[a_val, 'HGVSp'][:-1]:
                                #print('Gene is the same and ref is not the same' + '\t' + sample_df.loc[i, :] + '\t' + df_res.loc[a, :])
                                part_file.write('Gene is the same and ref is not the same' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                            if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][1:] == df_res.loc[a_val, 'HGVSp'][1:]:
                               #print('Gene is the same and ref is different, gene change is the same' + '\t' + sample_df.loc[i, :] + '\t' + df_res.loc[a, :])
                               part_file.write('Gene is the same and ref is different, gene change is the same' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                    if df_res.loc[a_val, :].any() == "NaN" or df_res.loc[a_val, :].any() == "NA" or df_res.loc[a_val, :].any() == "nan" or df_res.loc[a_val, :].any() == "na":
                        na_counter += 1

            if len(a) == 1:
                a_val = a[0]
                if df_res.loc[a_val, 'IMPACT'] != 'MODIFIER' and str(df_res.loc[a_val, :].any()) != "NaN" and str(df_res.loc[a_val, :].any()) != None and a_val not in index_a:
                    b = df_res.loc[a_val, 'HGVSp']
                    if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and v == b:
                        index_a.append(a_val)
                        countcord += 1
                        index_found.append(i)
                        # checks if both the variants are the same
                        out_file.write(str(len(index_found)) + '\t' +'\t'.join(sample_df.loc[i,:]) + "\t" + '\t'.join(df_res.loc[a_val,:]) + "\t" + 'Concordant' + "\n")
                        new_line = str(len(index_found)) + '\t' +'\t'.join(sample_df.loc[i,:]) + "\t" + '\t'.join(df_res.loc[a_val,:]) + "\t" + 'Concordant' + "\n"
                    else:
                        if sample_df.loc[i, 'HGVSp'][:-1] == df_res.loc[a_val, 'HGVSp'][:-1] and sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene']:
                            #print('Gene is the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a, :]))
                            part_file.write('Gene is the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                        if sample_df.loc[i, 'Gene'] != df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][:-1] == df_res.loc[a_val, 'HGVSp'][:-1]:
                            #print('Gene is not the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a, :]))
                            part_file.write('Gene is not the same and ref is the same: ' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                        if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][:-1] != df_res.loc[a_val, 'HGVSp'][:-1]:
                            #print('Gene is the same and ref is not the same' + '\t' + sample_df.loc[i, :] + '\t' + df_res.loc[a, :])
                            part_file.write('Gene is the same and ref is not the same' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
                        if sample_df.loc[i, 'Gene'] == df_res.loc[a_val, 'Gene'] and sample_df.loc[i, 'HGVSp'][1:] == df_res.loc[a_val, 'HGVSp'][1:]:
                            #print('Gene is the same and ref is different, gene change is the same' + '\t' + sample_df.loc[i, :] + '\t' + df_res.loc[a, :])
                            part_file.write('Gene is the same and ref is different, gene change is the same' + '\t' + '\t'.join(sample_df.loc[i, :]) + '\t' + '\t'.join(df_res.loc[a_val, :]) + '\n')
    df_new = df_res[~df_res.index.isin(index_a)]
    count = len(index_found)
    bad_df = df_res.index.isin(index_a)
    df_new = df_res[~bad_df]
    df_new_variants = list(df_new.loc[:, 'HGVSp'])
    # use these variants and genes to look up the missing values in the MSKCC table and print out the impact for it. and also the whole row in txt file
    mskcc_not_found_variants = sample_df[~sample_df.index.isin(index_found)]
    print('length of the variants that were not found in MSKCC: ', len(list(mskcc_not_found_variants.loc[:, 'HGVSp'])))
    print('Length of RES sample: ', len(df_res))
    print('Length of MSKCC sample: ', len(df_mskcc))
    print('Number of MSKCC indexes found in Dragen (based on index values): ', count)
    print('Number of Dragen indexes found in MSKCC (based on index values): ', len(index_a))
    print('Number of Dragen indexes not found in MSKCC (based on indexes not found): ', len(df_new))
    print('Number of variants not found in dragen: ', len(df_new.loc[:, 'HGVSp']))
    impact_counter = Counter(mskcc_not_found_variants.loc[:, 'Type'])
    print('IMPACT of the variants not found: ', dict(impact_counter))
    mskcc_not_found_variants = mskcc_not_found_variants.mask(mskcc_not_found_variants.eq('None')).dropna()
    for values in mskcc_not_found_variants.index:
        if type(mskcc_not_found_variants.loc[values, :].all()) == str:
            #mskcc_not_found_variants.loc[values, 'HGVSc'] = str( "c." + str(mskcc_not_found_variants.loc[values, 'HGVSc'].split("c.")[1]))
            out_file_noncordant.write(str('\t'.join(list(mskcc_not_found_variants.loc[values, :]))))
            #out_file_noncordant.write('\n')
    return mskcc_not_found_variants, index_found, index_a, df_new, impact_counter

if __name__ == "__main__":
    import sys
    import pandas as pd
    import re
    import subprocess
    import os
    from collections import Counter
    import time
    import numpy as np
    start = time.time()
    filenames_res_outputs = os.listdir('/gpfs/commons/groups/clinical/vshah/concordance_test/no_multi')
    filenames_mskcc = os.listdir('/gpfs/commons/groups/clinical/vshah/concordance_test/samples_mskcc')
    numbers_res = []
    new_df = open("new_df_conc_full_result.tsv", 'w')
    new_df.write("Sample" + "\t" + "MSKCC Count" + "\t" + "Dragen Count" + "\t" + "Perc concordance" + "\t" + "Non cordant" + "\t" + "IMPACT"+ "\n")
    for file1 in filenames_res_outputs:
        split1 = file1.split("-")
        numbers_res.append(split1[1].split("T")[0])
    for file2 in filenames_mskcc:
        split2 = file2.split("_")
        num = split2[1]
    #    numbers_res = []
    #    numbers_res = ['03', '06']
        if num in numbers_res:
            print(num)
            sample_id = 'RES20-' + num + 'T-1-D1'
            input_res_vep = "/gpfs/commons/groups/clinical/vshah/concordance_test/no_multi/RES20-" + num + "T-1-D1--RES20-" + num + "N-1-D1.hard-filtered.annotated_filter_vep_table.txt"
            df_res = pd.read_csv(input_res_vep, sep="\t")
            print('csv read')
            df_res = df_res.loc[:, ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'Consequence', 'IMPACT', 'SYMBOL', 'HGVSc', 'HGVSp', str(sample_id + '.AF'), str(sample_id + '.SQ')]]
            input_mskcc = '/gpfs/commons/groups/clinical/vshah/concordance_test/samples_mskcc/' + file2
            df_res.columns = ['CHR', 'POS', 'REF', 'ALT', 'FILTER', 'Type','IMPACT', 'Gene', 'HGVSc', 'HGVSp', 'VAF', 'SQ']
            df_mskcc = pd.read_csv(input_mskcc, sep="\t", low_memory=False)
            df_res_parsed = parse_res_output(df_res)
            df_mskcc_parsed = parse_sample(df_mskcc)
            mskcc_not_found_variants, index_found_mskcc, index_found_res, df_new, impact_counter = merge_dfs(df_mskcc_parsed, df_res_parsed, sample_id)
            new_df.write(str(sample_id) + "\t" + str(len(df_mskcc_parsed)) + "\t" + str(len(index_found_mskcc)) + "\t" + str(len(index_found_mskcc)/len(df_mskcc_parsed)) + "\t" + str(len(mskcc_not_found_variants)) + "\t" + str(impact_counter) + "\n")
    new_df.close()
    end = time.time()
    print('time taken to run', end-start)
