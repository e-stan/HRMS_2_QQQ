import DecoID.DecoID as DecoID
import numpy as np
import os

#parameters
polarity = "Negative"
datadirs = ["../data/IDX/AX_stab_pool/neg/","../data/IDX/AX_stab_pool/pos/"]

libFiles = ["D:/github_clean/DecoID/databases/mzCloud_reference.db","D:/github_clean/DecoID/databases/IROA_updated_NCE40.db"]

if __name__ == "__main__":
    for datadir in datadirs:
        filenamesToCombine = []
        for libFile,lab in zip(libFiles,["mzCloud","in-house"]):
            decoID = DecoID.DecoID(libFile, "reference", numCores=20, label=lab, scoringFunc=DecoID.dotProductSpectra)
            filenamesToCombine += [datadir + x.replace(".mzML","") + lab for x in os.listdir(datadir) if ".mzML" in x]
            for filename in [x for x in os.listdir(datadir) if ".mzML" in x]:
                print(datadir + filename)
                peakfile = datadir + "peaks.csv"
                decoID.readData(datadir + filename,2,True,True,10,peakDefinitions=peakfile,frag_cutoff=0)
                decoID.identifyUnknowns(resPenalty=5.0,iso=True,rtTol=np.inf,dpThresh=80)
                decoID.searchSpectra("y",resPenalty=5.0,iso=True,rtTol=np.inf)
        decoID.combineResultsAutoAnnotate(filenamesToCombine,datadir + "combined_resultsTOP3_above50.csv",numHits=3,min_score=50)