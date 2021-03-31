from srm_helper import *
import pandas as pd

import numpy as np

### format csv files from whole transition list
totalTransitionInfoFn = "../data/IDX/6495_m3t_transitions.csv"

totalTransitions = pd.read_csv(totalTransitionInfoFn)

#totalTransitions["rt_start"] = totalTransitions["rt_start"] - .5
#totalTransitions["rt_end"] = totalTransitions["rt_end"] + .5

switcher = {"Positive":1,"Negative":-1}

filt = totalTransitions[totalTransitions["Training"] == 1]
filt.to_csv("../data/IDX/target_transitions_to_learn_conv_6495.csv",index=False)

targets = totalTransitions
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
targets.to_csv("../data/IDX/targets_for_evaluation_6495.csv",index=False)

if __name__ == "__main__":
    # params
    tol = .5  # MS2 fragment tolerance for QqQ optimized transitions
    ppmTol = 10  # m/z tolerance for HRMS data in ppm
    numCores = 2  # number of CPU cores to use

    # create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol, numCores=numCores)

    # set datafiles for learning conversion
    trainingData = pd.read_csv("../data/IDX/target_transitions_to_learn_conv_6460.csv")

    msFilenames = ["../data/IDX/IDX_MS2_data/M3T_10uM_pos_DDA_10NCEs_25-35_50ms_5e4_DE5s_updatedRT.mzML",
                   "../data/IDX/IDX_MS2_data/M3T_10uM_neg_DDA_10NCEs_25-35_50ms_5e4_DE5s_updatedRT.mzML",
                   "../data/IDX/IDX_MS2_data/M3T_10uM_pos_DDA_10NCEs_25-35_80ms_1e4_DE5s_updatedRT_missing.mzML"]

    #build conversion
    merged = srm_maker.buildConversion(msFilenames, trainingData, tic_cutoff=0, frag_cutoff=0,
                                       frag_ppm_tolerance=2 * 1e6 * .5 / 200)
    merged.to_csv("../data/IDX/conversion_results_6495.csv")

    # output conversion
    print(srm_maker.getConversionEquationString())

    # set datafiles to build srms
    targets = pd.read_csv("../data/IDX/targets_for_evaluation_6495.csv")

    srm_table = pd.DataFrame()
    breakdownCurves = {}

    for msFilename in msFilenames:

        # create SRM table
        srm_table1, breakdownCurves1 = srm_maker.createSRMsCE(msFilename, targets)

        srm_table = pd.concat((srm_table,srm_table1),axis=0,ignore_index=True)
        breakdownCurves.update(breakdownCurves1)

    # output SRM file
    srm_table.to_csv("../data/IDX/generated_SRM_table_6495.csv")

    # plot breakdown curves
    plotBreakdownCurves(breakdownCurves,"../data/IDX/breakdown_curves_6495.pdf")