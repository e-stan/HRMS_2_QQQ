from srm_helper import *
import pandas as pd

if __name__ == "__main__":
    # params
    tol = .5  # MS2 fragment tolerance for QqQ optimized transitions
    ppmTol = 10  # m/z tolerance for HRMS data in ppm
    numCores = 2  # number of CPU cores to use

    # create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol, numCores=numCores)

    # set datafiles for learning conversion
    trainingData = pd.read_csv("target_transitions_to_learn_conv.csv")

    msFilenames = ["training.mzML"]
    print(msFilenames)

    # build conversion
    merged = srm_maker.buildConversion(msFilenames, trainingData, tic_cutoff=0, frag_cutoff=0,
                                       frag_ppm_tolerance=2 * 1e6 * .5 / 200)
    merged.to_csv("conversion_results.csv")

    # output conversion
    print(srm_maker.getConversionEquationString())

    # set datafiles to build srms
    targets = pd.read_csv("targets.csv")

    # filename for HRMS MS/MS of targets
    msFilename = "targets.mzML"

    # create SRM table
    srm_table, breakdownCurves = srm_maker.createSRMsCE(msFilename, targets)

    # output SRM file
    srm_table.to_csv("generated_SRM_table.csv")

    # plot breakdown curves
    plotBreakdownCurves(breakdownCurves,"breakdown_curves.pdf")