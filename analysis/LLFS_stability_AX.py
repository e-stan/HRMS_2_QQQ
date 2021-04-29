from srm_helper import *
import pandas as pd

if __name__ == "__main__":
    # params
    tol = .5  # MS2 fragment tolerance for QqQ optimized transitions
    ppmTol = 10  # m/z tolerance for HRMS data in ppm
    numCores = 5  # number of CPU cores to use

    # create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol, numCores=numCores)

    # set datafiles for learning conversion
    trainingData = pd.read_csv("../data/IDX/M3T_transitions_ALTIS_optimized_allCpds.csv")

    msFilenames = ["../data/IDX/IDX_MS2_data/M3T_10uM_pos_DDA_10NCEs_25-35_50ms_5e4_DE5s_updatedRT.mzML",
                   "../data/IDX/IDX_MS2_data/M3T_10uM_neg_DDA_10NCEs_25-35_50ms_5e4_DE5s_updatedRT.mzML",
                   "../data/IDX/IDX_MS2_data/M3T_10uM_pos_DDA_10NCEs_25-35_80ms_1e4_DE5s_updatedRT_missing.mzML"]

    # build conversion
    merged = srm_maker.buildConversion(msFilenames, trainingData, tic_cutoff=0, frag_cutoff=0,
                                       frag_ppm_tolerance=2 * 1e6 * .5 / 200)

    merged.to_csv("../data/IDX/AX_stab_pool/conversion_results_AX.csv")

    # output conversion
    print(srm_maker.getConversionEquationString())

    # set datafiles to build srms
    targets = pd.read_csv("../data/IDX/AX_stab_pool/combined_peak_list.csv")

    # filename for HRMS MS/MS of targets
    msFilenames = ["../data/IDX/AX_stab_pool/pos/"+x for x in os.listdir("../data/IDX/AX_stab_pool/pos/") if ".mzML" in x] + ["../data/IDX/AX_stab_pool/neg/"+x for x in os.listdir("../data/IDX/AX_stab_pool/neg/") if ".mzML" in x]

    print(msFilenames)

    srm_table = pd.DataFrame()
    breakdownCurves = {}

    for msFilename in msFilenames:

        # create SRM table
        srm_table1, breakdownCurves1 = srm_maker.createSRMsCE(msFilename, targets)

        srm_table = pd.concat((srm_table,srm_table1),axis=0,ignore_index=True)
        breakdownCurves.update(breakdownCurves1)

    # output SRM file
    srm_table.to_csv("../data/IDX/AX_stab_pool/generated_SRM_table.csv")

    # plot breakdown curves
    plotBreakdownCurves(breakdownCurves,"../data/IDX/AX_stab_pool/breakdown_curves.pdf")