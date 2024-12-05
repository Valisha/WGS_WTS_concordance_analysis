#!./python3

import pandas as pd
import os
tmb_list = os.listdir()
TMB = []
TMB_nonsyn = []
sample = []
for tmb in tmb_list:
    if "tmb.metrics.csv" in tmb:
        df = pd.read_csv(tmb, sep= ",")
        s2 = tmb.split("T")[0]
        sample.append(s2)
        df.columns = ['SECTION', '', 'TMB', 'value']
        tmb_index = df.index[df['TMB'] == 'TMB']
        TMB.append(df.loc[tmb_index[0], 'value'])
        nonsyn_index = df.index[df['TMB'] == 'Nonsyn TMB']
        TMB_nonsyn.append(df.loc[nonsyn_index[0], 'value'])
df_tmb = pd.DataFrame()
df_tmb['Sample'] = sample
df_tmb['TMB'] = TMB
df_tmb['TMB_nonsyn'] = TMB_nonsyn
#df_tmb.to_csv("df_tmb_res20.tsv", sep="\t")

df_mskcc = pd.read_csv("df_mskcc.tsv", sep="\t")

df_msi = pd.read_csv("df_msi.tsv", sep="\t")

df1 = df_mskcc.merge(df_msi, left_on="NYGC-ID", right_on="Sample")
df2 = df1.merge(df_tmb, left_on='Sample', right_on='Sample')
df2.to_csv("msi_tmb_conc_RES_MSKCC.tsv", sep="\t")

