#! /nfs/sw/python/python-3.6.1/bin/python3

import pandas as pd
import os

sample_ids = []
overall_ploidy = []
estimated_tumor_purity =[]

dragen = os.listdir("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/")
for file_dragen in dragen:
    files = os.listdir(("/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/" + file_dragen))
    for file1 in files:
        if '.cnv_metrics.csv' in file1:
            sample_ids.append(file_dragen)
            dataframe = "/gpfs/commons/groups/clinical/Project_CLIV_13843_B01_SOM_WGS_AWS_Dragen/" + file_dragen + "/" + file1
            data = pd.read_csv(dataframe, sep=",")
            data.columns = ['summary type', 'blank', 'metrics', 'value', 'perc']
            data_new = data.loc[:, ['metrics', 'value']]
            for i, v in enumerate(data_new['metrics']):
                if v == 'Overall ploidy':
                    overall_ploidy.append(data_new.loc[i, 'value'])
                if v == 'Estimated tumor purity':
                    estimated_tumor_purity.append(data_new.loc[i, 'value'])

compiled_cnv_conc = pd.DataFrame()
compiled_cnv_conc['sample'] = sample_ids
compiled_cnv_conc['overall_ploidy'] = overall_ploidy
compiled_cnv_conc['estimated_tumor_purity'] = estimated_tumor_purity
#print(len(sample_ids))
#print(len(overall_ploidy))
#print(len(estimated_tumor_purity))
compiled_cnv_conc.to_csv("compiled_cnv_conc.tsv", sep="\t")
