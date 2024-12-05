#!./python3

import pandas as pd
import os

df_mskcc = pd.read_csv('/gpfs/commons/groups/clinical/vshah/tmb_corr/nygc_mskcc_res_report.csv', sep=',')
df_mskcc = df_mskcc.iloc[0:49, 0:7]
df_mskcc = df_mskcc.loc[:, ['Sample-MSKCC', 'MSI Status', 'TMB (mt/Mb)', 'NYGC-ID', 'MSI Sensor Score']]
df_mskcc.to_csv('df_mskcc.tsv', sep='\t')
