# HRMS_2_QQQ

Metabolomics software for the automated developement of SRM tables for HRMS MS/MS Data

## System Requirements

Package has been tested with Python 3.7 on Windows 10

Package has the following dependencies:

numpy (v1.18.1)

sklearn (v0.22.1)

pandas (v1.0.1)

DecoID (v0.2.8)

matplotlob (v3.1.3)

In order to process vendor formatted data without manual conversion, MS-Convert (http://proteowizard.sourceforge.net/tools.shtml) needs to be installed and added to PATH. 

## Installation

### Installation with ```pip```:

```
pip install srm_helper
```
PyPI:

https://pypi.org/project/srm_helper/

### Manual installation from source:

```
git clone https://github.com/e-stan/HRMS_2_QQQ.git
pip install DecoID/src/.
```

## Demo

Demo data available under HRMS_2_QQQ/examples/

Example usage to build a conversion between two instruments and apply to list of targets:

```
from srm_helper import *
import pandas as pd

if __name__ == "__main__":
    #params
    tol = .5 #MS2 fragment tolerance for QqQ optimized transitions
    ppmTol = 10 #m/z tolerance for HRMS data in ppm
    numCores = 2 #number of CPU cores to use

    #create srm_maker object
    srm_maker = SRM_maker(ppm=ppmTol,numCores=numCores)

    #set datafiles for learning conversion
    trainingData = pd.read_csv("target_transitions_to_learn_conv.csv")
    
    msFilenames = ["training.mzML"]
    print(msFilenames)

    #build conversion
    merged = srm_maker.buildConversion(msFilenames,trainingData,tic_cutoff=0,frag_cutoff=0,frag_ppm_tolerance=2 * 1e6 * .5/200)
    merged.to_csv("conversion_results.csv")

    #output conversion
    print(srm_maker.getConversionEquationString())

    #set datafiles to build srms
    targets = pd.read_csv("targets.csv")

    #filename for HRMS MS/MS of targets
    msFilename = "targets.mzML"

    #create SRM table
    srm_table,breakdownCurves = srm_maker.createSRMsCE(msFilename,targets)

    #output SRM file
    srm_table.to_csv("generated_SRM_table.csv")

    #plot breakdown curves
    plotBreakdownCurves(breakdownCurves,"breakdown_curves.pdf")

```

expected output files are included in the examples directory




