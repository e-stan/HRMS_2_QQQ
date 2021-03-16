from src.srm_helper import *
import pandas as pd
import os
import pickle as pkl

import multiprocessing.pool

class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NestablePool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(NestablePool, self).__init__(*args, **kwargs)

datadir = "data/IDX/IDX_MS2_data/"

polarity = "neg"

peakInfo = "data/IDX/IDX_rts_M3T_neg.csv"

api_key = "pnlPwLgBGycAFZW7Dany2LwkDy265u"

libFile = "D:/github_clean/DecoID/databases/HMDB_experimental.db"

rfFile = "data/IDX/areas_rf_" + polarity + ".csv"

peakInfo = pd.read_csv(peakInfo)

rfData = pd.read_csv(rfFile,index_col=0)

numCores = 15

rt_start = .3
rt_end = 1.0

growingSRM = pd.DataFrame()
growingBreakdown = {}


if __name__ == "__main__":
    args = []
    #multiprocessing.set_start_method("spawn")
    p = NestablePool(numCores)
    for index,row in peakInfo.iterrows():
        filename = row["Name"] + "_ddMS2_" + polarity + ".mzML"
        if os.path.isfile(datadir + filename):
            tmpPeak = peakInfo[peakInfo["Name"] == row["Name"]]
            tmpPeak.at[index,"rt_start"] = rt_start
            tmpPeak.at[index,"rt_end"] = rt_end
            args.append([datadir+filename,tmpPeak,10,libFile,1,api_key,0,0,13.9,5])


    results = p.starmap(createSRMsCE,args)
    p.close()
    p.join()

    for srmTable,breakdownCurves in results:
        if len(srmTable) > 0:
            growingSRM = pd.concat((growingSRM,srmTable),ignore_index=True,axis=0)
            growingBreakdown.update(breakdownCurves)

    growingSRM = addRFtoSRM(growingSRM,rfData)

    #replace RT
    for index,row in peakInfo.iterrows():
        tmp = growingSRM[growingSRM["Name"] == row["Name"]]
        for index2,row2 in tmp.iterrows():
            growingSRM.at[index2,"rt_start"] = row["rt_start"]
            growingSRM.at[index2,"rt_end"] = row["rt_end"]



    growingSRM.to_csv("data/IDX/SRM_IDX_M3T_" + polarity + ".csv")

    pkl.dump(growingBreakdown,open("data/IDX/IDX_M3T_" + polarity + '.pkl',"wb"))

    plotBreakdownCurves(growingBreakdown,"data/IDX/breakdown_IDX_M3T_" + polarity + ".pdf")

#growingBreakdown = pkl.load(open("data/IDX/IDX_M3T_" + polarity + '.pkl', "rb"))

#plotBreakdownCurves(growingBreakdown, "data/IDX/breakdown_IDX_M3T_" + polarity + ".pdf")