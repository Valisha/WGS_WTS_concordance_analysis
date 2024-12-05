#!./python3

import pandas as pd
import os
import json
msi_list = os.listdir()

#"TotalMicrosatelliteSitesAssessed": "45819",
#"TotalMicrosatelliteSitesUnstable": "596",
#"PercentageUnstableSites"
sites_assessed = []
sites_unstable = []
perc_unstable = []
sample = []
for msi in msi_list:
    if "microsat_output.json" in msi: 
        s2 = msi.split("T")[0]
        sample.append(s2)
        f = open(msi, 'r')
        data = json.load(f)
        sites_assessed.append(data["TotalMicrosatelliteSitesAssessed"])
        sites_unstable.append(data["TotalMicrosatelliteSitesUnstable"])
        perc_unstable.append(data["PercentageUnstableSites"])

new_df = pd.DataFrame()
new_df["Sample"] = sample
new_df["Microsites_assessed"] = sites_assessed
new_df["Microsites_unstable"] = sites_unstable
new_df["Perc_of_unstable_sites"] = perc_unstable
new_df["Perc_of_unstable_sites"] = pd.to_numeric(new_df["Perc_of_unstable_sites"], downcast = "float")
new_df["Perc_of_unstable_sites"] = new_df["Perc_of_unstable_sites"].round(decimals=2)
new_df.to_csv("df_msi.tsv", sep="\t")

