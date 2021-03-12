import os
import numpy as np
from difflib import SequenceMatcher
from DecoID.DecoID import DecoID
import pandas as pd
from copy import deepcopy
import uuid
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.linear_model import LinearRegression


class SRM_maker():

    def __init__(self,ppm=10,numCores=2,CE_converter=lambda x: x,ms2_resolution = 2):
        self.ppm = ppm
        self.numCores = numCores
        self.CE_converter = CE_converter
        self.coeffs = [0,1,0]
        self.ms2_resolution = ms2_resolution


    def readMSMSData(self,msFile,targets,tic_cutoff,frag_cutoff):
        # make DecoID object
        uid = str(uuid.uuid1())
        tmpLib = uid + ".msp"
        open(tmpLib, "w")

        decID = DecoID(tmpLib, "reference", self.numCores, self.ms2_resolution)

        # write temporary peak file
        targets.to_csv(uid + ".csv", index=False)

        # read in file and save all spectra
        decID.readData(msFile,self.ms2_resolution, True, True, self.ppm, peakDefinitions=uid + ".csv", tic_cutoff=tic_cutoff,
                       frag_cutoff=frag_cutoff)
        samplesAll = deepcopy(decID.samples)

        # get charge
        polarity = samplesAll[0]["mode"]
        switcher = {"Positive": 1, "Negative": -1}
        polarity = switcher[polarity]

        # structure to hold spectra
        output_dict = {}

        # process spectra for each CE
        if len(decID.samples) > 0:

            # get unique CEs
            ces = list(set([x["CE"] for x in samplesAll]))
            ces.sort()

            # iterate over CEs
            for ce in ces:

                # parse relevatn samples
                decID.samples = [x for x in samplesAll if x["CE"] == ce]
                decID.label = str(ce)

                # deconvolve/process spectra
                decID.searchSpectra("n", np.inf, iso=False, rtTol=np.inf)

                # read in results
                results = pd.read_csv(msFile.replace(".mzML", decID.label + "_scanInfo.csv"))
                ind = 0

                # iterate over targets
                for _, row in targets.iterrows():
                    # parse results for single compound
                    tmp = results[results["#featureID"] == ind]

                    # get processed spectra
                    tmp = tmp[tmp["componentID"] == "original"]
                    if len(tmp) > 0:
                        spectrum = tmp.at[tmp.index.values[0], "spectrum"]
                        spectrum = {float(x.split(":")[0]): float(x.split(":")[1]) for x in spectrum.split()}
                        rts = [x["rt"] for x in decID.samples if x["group"] == ind]

                        # compute interpolated MS1 intensity
                        func = extractChromatogram(row["mz"], decID.ms1, [row["rt_start"], row["rt_end"]], self.ppm)
                        normalizer = np.sum([func(x) for x in rts])

                        # normalize spectra and save result
                        spectrum = {key: val / normalizer for key, val in spectrum.items()}
                        if row["Name"] not in output_dict:
                            output_dict[row["Name"]] = {}
                        output_dict[row["Name"]][ce] = spectrum

                        ind += 1

                # delete results
                os.remove(msFile.replace(".mzML", decID.label + "_scanInfo.csv"))
                os.remove(msFile.replace(".mzML", decID.label + "_decoID.csv"))
                os.remove(msFile.replace(".mzML", decID.label + ".DecoID"))


        os.remove(uid + ".csv")
        os.remove(uid + ".msp")

        return output_dict,polarity

    def getConversionEquationString(self):
        return "CE = " + str(self.coeffs[0]) + "*mz + " + str(self.coeff[1]) + "*ce + " + str(self.coeff[2])

    def buildConversion(self,msFile,targetTransitions,tic_cutoff=0,frag_cutoff=0,frag_ppm_tolerance=10):

        spectra,_ = self.readMSMSData(msFile,targetTransitions,tic_cutoff,frag_cutoff)
        optimals = {}
        for index,row in targetTransitions.iterrows():
            optimals[row["Name"]] = {"HRMS_CE":-1}
            if row["Name"] in spectra:
                tmp = {}
                for ce,spec in spectra[row["Name"]]:
                    for mz2,i in spec:
                        if 1e6 * np.abs(mz2-row["Product mz"])/row["Product mz"] < frag_ppm_tolerance:
                            if ce not in tmp:
                                tmp[ce] = 0
                            tmp[ce] += i
                keys = list(tmp.keys())
                keys.sort(key=lambda x: tmp[x],reverse=True)
                optimals["HRMS_CE"] = keys[0]

        merged = pd.concat((targetTransitions,pd.DataFrame.from_dict(optimals,orient="index")),axis=1)
        merged = merged[merged["HRMS_CE"] >= 0]
        X = np.array([merged["mz"].values,merged["HRME_CE"],[1 for _ in range(len(merged))]]).transpose()
        y = merged["CE"].values
        linreg = LinearRegression(fit_intercept=False)
        linreg.fit(X,y)
        self.CE_converter = lambda x: linreg.predict([[x[0],x[1],1]])[0]
        self.coeffs = linreg.coef_
        return merged


    def createSRMsCE(self,msFile,targets,tic_cutoff=0,frag_cutoff=0,lowestProduct = 13.9,numProduct=5,outfile="none"):

        output_dict,polarity = self.readMSMSData(msFile,targets,tic_cutoff,frag_cutoff)

        if len(output_dict) > 0:
            #find best fragments
            srm = {}
            breakdown_curves = {}
            ind = 0

            for cpd in output_dict:
                frag_int_dict = {}
                #get precursor mz
                precursorMz = targets[targets["Name"] == cpd]
                precursorMz = precursorMz.at[precursorMz.index.values[0],"mz"]

                #iterate over spectra at different ces
                for ce in output_dict[cpd]:
                    #collect fragment intensities and ces
                    for mz,i in output_dict[cpd][ce].items():
                        if precursorMz-mz > lowestProduct:
                            if mz not in frag_int_dict:
                                frag_int_dict[mz] = {}
                            frag_int_dict[mz][ce] = i

                #sort fragments by highest observed intensity
                topFrags = list(frag_int_dict.keys())
                topFrags.sort(key=lambda x:np.max(list(frag_int_dict[x].values())),reverse=True)

                #take top fragments
                topFrags = topFrags[:numProduct]
                breakdown_curves[cpd] = {mz:breakdown for mz,breakdown in frag_int_dict.items() if mz in topFrags}

                #get best ces
                cesOpt = []
                for frag in topFrags:
                    tmp = list(breakdown_curves[cpd][frag].keys())
                    tmp.sort(key = lambda x: breakdown_curves[cpd][frag][x],reverse=True)
                    cesOpt.append(tmp[0])

                #addToSRM
                tmp = targets[targets["Name"] == cpd]
                cols = tmp.columns.values
                i = tmp.index.values[0]

                for frag,ce in zip(topFrags,cesOpt):
                    srm[ind] = {c:tmp.at[i,c] for c in cols}
                    srm[ind]["Product mz"] = frag
                    srm[ind]["CE"] = self.CE_converter(ce)
                    srm[ind]["Normalized Intensity"] = breakdown_curves[cpd][frag][ce]
                    srm[ind]["Charge"] = polarity
                    ind += 1
            srm = pd.DataFrame.from_dict(srm,orient="index")

        else:
            srm = pd.DataFrame()
            breakdown_curves = {}

        if outfile != "none":
            srm.to_csv(outfile)

        return srm,breakdown_curves


def plotBreakdownCurves(breakdown_curves,filename):
    pp = PdfPages(filename)

    for cpd in breakdown_curves:
        pllt = plt.figure()
        plt.title(cpd)
        offset = 0
        index = 1
        offsets = []
        xtick = []
        labels = []
        fragments = list(breakdown_curves[cpd].keys())
        fragments.sort(key=lambda x: list(np.max(breakdown_curves[cpd][x].values())),reverse=True)
        for f in fragments:
            label = "Fragment " + str(index) + ": "
            CEs = list(breakdown_curves[cpd][f].keys())
            CEs.sort()
            index += 1
            plt.bar([x + offset for x in range(len(CEs))], [breakdown_curves[cpd][f][c] for c in CEs],
                    label=label + str(np.round(f, 2)))
            offsets.append(offset)
            xtick += [x + offset for x in range(len(CEs))]
            labels += [str(ce) for ce in CEs]
            offset += len(CEs) + 1


        plt.xticks(xtick, labels, rotation=45, size=7)
        plt.xlabel("NCE")
        plt.ylabel("Normalized Intensity")
        plt.legend()
        plt.tight_layout()

        pp.savefig(pllt)
    pp.close()


def addRFtoSRM(srmtable,RF_peak_areas):
    rfs = RF_peak_areas.columns.values

    for index,row in RF_peak_areas.iterrows():
        #get best RF
        tmp = list(rfs)
        tmp.sort(key=lambda x:RF_peak_areas.at[index,x],reverse=True)
        bestRF = tmp[0]

        #get relevant srm tables indices
        tmp = srmtable[srmtable["Name"] == index]

        for i in tmp.index.values:
            srmtable.at[i,"RF"] = bestRF

    return srmtable


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def extractChromatogram(mz, ms1, rtRange,mError):
    allRt = [x for x in list(ms1.keys()) if x >= rtRange[0] and x <= rtRange[1]]
    allRt.sort()
    chromatogram = {}
    allMzList = list()
    for s in ms1:
        allMzList += [x for x in ms1[s] if ms1[s][x] > 0]
    allMzList = np.array(list(set(allMzList)))


    mzOI = np.extract(np.abs((10 ** 6) * (allMzList - mz)) / mz < mError, allMzList)

    getInt = lambda rt: np.sum([ms1[rt][x] for x in mzOI if x in ms1[rt]])

    tempMs1 = []
    for rt in allRt:
        tempMs1.append([rt, getInt(rt)])

    tempMs1 = np.array(tempMs1)
    func = lambda x: np.interp([x],tempMs1[:,0],tempMs1[:,1])
    return func


