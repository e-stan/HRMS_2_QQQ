from srm_helper import *
import pandas as pd
import numpy as np

### format csv files from whole transition list
totalTransitionInfoFn = "../data/IDX/M3T_transitions_ALTIS_optimized_allCpds.csv"

totalTransitions = pd.read_csv(totalTransitionInfoFn)

switcher = {"Positive":1,"Negative":-1}

totalTransitions["rt_start"] = np.array([.3 for _ in range(len(totalTransitions))])
totalTransitions["rt_end"] = np.array([1.0 for _ in range(len(totalTransitions))])

filt = totalTransitions[totalTransitions["Training"] == 1]
filt.to_csv("../data/IDX/target_transitions_to_learn_conv.csv",index=False)

targets = totalTransitions[totalTransitions["Training"] == 0]
cpds = []
toDrop = []
for index,row in targets.iterrows():
    if (row['Name'],row["Charge"]) in cpds:
        toDrop.append(index)
    else:
        cpds.append((row['Name'],row["Charge"]))

targets = targets.drop(toDrop)
goodCols = ["Name","rt_start","rt_end","mz","Charge"]
targets = targets[goodCols]
targets.to_csv("../data/IDX/targets_for_evaluation.csv",index=False)

if __name__ == "__main__":
    #params
    tol = .5
    ppmTol = 10
    numCores = 20

    #create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol,numCores=numCores)

    #set datafiles for learning conversion
    trainingData = pd.read_csv("../data/IDX/target_transitions_to_learn_conv.csv")
    msFilenames = ["../data/IDX/IDX_MS2_data/training_cpds/" + x for x in os.listdir(
        "../data/IDX/IDX_MS2_data/training_cpds/") if ".mzML" in x and "merged" not in x]
    print(msFilenames)

    #build conversion
    merged = srm_maker.buildConversion(msFilenames,trainingData,tic_cutoff=0,frag_cutoff=0,frag_ppm_tolerance=2 * 1e6 * .5/200)
    merged.to_csv("../data/IDX/target_transitions_to_learn_conv_optimizedValues.csv")

    #output conversion
    print(srm_maker.getConversionEquationString())

    #set datafiles to build srms
    targets = pd.read_csv("../data/IDX/targets_for_evaluation.csv")

    #data structures
    breakdownCurves = {}
    srm_table = pd.DataFrame()


    #construct SRM tables for individual files
    datadir = "../data/IDX/IDX_MS2_data/"
    for index,row in targets.iterrows():
        if row["Charge"] > 0:
            polarity = "pos"
        else:
            polarity = "neg"
        msFilename = datadir + row["Name"] + "_ddMS2_" + polarity + ".mzML"
        t = targets.loc[index:index+1,:]
        tmp_srm,tmp_breakdownCurves = srm_maker.createSRMsCE(msFilename,t)
        breakdownCurves.update(tmp_breakdownCurves)
        srm_table = pd.concat((srm_table,tmp_srm),axis=0,ignore_index=True)

    #output SRM file
    srm_table.to_csv(datadir + "generated_SRM_table.csv")

    #plot breakdown curves
    plotBreakdownCurves(breakdownCurves,datadir + "breakdown_curves.pdf")