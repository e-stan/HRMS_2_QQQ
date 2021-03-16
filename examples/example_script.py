from src.srm_helper import *
import pandas as pd
import numpy as np

### format csv files from whole transition list
totalTransitionInfoFn = "data/IDX/M3T_transitions_ALTIS_optimized.csv"

totalTransitions = pd.read_csv(totalTransitionInfoFn)

switcher = {"Positive":1,"Negative":-1}
for polarity in ["Positive","Negative"]:
    filt = totalTransitions[totalTransitions["Charge"] == switcher[polarity]]
    filt["rt_start"] = np.array([.1 for _ in range(len(filt))])
    filt["rt_end"] = np.array([1.0 for _ in range(len(filt))])

    training = filt[filt["Training"] == 1]
    training.to_csv("data/IDX/target_transitions_to_learn_conv_"+polarity[:3].lower()+".csv",index=False)
    targets = filt[filt["Training"] == 0]
    cpds = []
    toDrop = []
    for index,row in targets.iterrows():
        if (row['Name'],row["Charge"]) in cpds:
            toDrop.append(index)
        else:
            cpds.append((row['Name'],row["Charge"]))
    targets = targets.drop(toDrop)
    goodCols = ["Name","Formula","rt_start","rt_end","mz"]
    targets = targets[goodCols]
    targets.to_csv("data/IDX/targets_for_evaluation_"+polarity[:3].lower()+".csv",index=False)

if __name__ == "__main__":
    #params
    polarity = "neg"
    tol = .5
    ppmTol = 10
    numCores = 20

    #create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol,numCores=numCores)

    #set datafiles for learning conversion
    trainingData = pd.read_csv("data/IDX/target_transitions_to_learn_conv_"+polarity+".csv")
    msFilename = "data/IDX/IDX_MS2_data/training_cpds/merged_" + polarity + ".mzML"

    #build conversion
    merged = srm_maker.buildConversion(msFilename,trainingData,tic_cutoff=0,frag_cutoff=0,frag_ppm_tolerance=2 * 1e6 * .5/200)

    #output conversion
    print(srm_maker.getConversionEquationString())

    #set datafiles to build srms
    targets = pd.read_csv("data/IDX/targets_for_evaluation_"+polarity+".csv")
    msFilename = "data/IDX/IDX_MS2_data/merged_" + polarity + ".mzML"

    #build and output srms
    srm_table,breakdownCurves = srm_maker.createSRMsCE(msFilename,targets)
    srm_table.to_csv(msFilename.replace(".mzML","_generated_SRM_table.csv"))

    #plot breakdown curves
    plotBreakdownCurves(breakdownCurves,msFilename.replace(".mzML","_breakdown_curves.pdf"))