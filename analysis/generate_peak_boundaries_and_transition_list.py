import pandas as pd
import os

dir = "../data/IDX/"

infoFilename = dir + "IDX_rts_M3T_pos.csv"

filenameDir = dir + "IDX_RF_data/"

charge = 1

info = pd.read_csv(infoFilename)

filenames = [x for x in os.listdir(filenameDir) if ".mzML" in x]

out_dict = {}
transition_list_dict = {}

i = 0
y = 0
for index,row in info.iterrows():
    transition_list_dict[y] = {"Molecule List Name":infoFilename.replace(dir,""),
                               "Precursor Name":row["Name"],
                               "Precursor Mz":row["mz"],
                               "Precursor Charge":charge,
                               "Explicit Retention Time":(row["rt_start"] + row["rt_end"])/2}

    y += 1
    for fn in filenames:
       out_dict[i] = {"File Name":fn,"Peptide Modified Sequence":row["Name"],"Min Start Time":row["rt_start"],"Max End Time":row["rt_end"]}
       i += 1

pd.DataFrame.from_dict(transition_list_dict,orient="index").to_csv(infoFilename.replace(".csv","_transitionListSL.csv"),index=False)
pd.DataFrame.from_dict(out_dict,orient="index").to_csv(infoFilename.replace(".csv","_peakBoundaries.csv"),index=False)

