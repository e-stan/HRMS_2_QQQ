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
from multiprocessing import Pool
import pickle as pkl



class SRM_maker():

    def __init__(self,ppm=10,numCores=2,CE_converter=lambda x: x[1],ms2_resolution = 2):
        self.ppm = ppm
        self.numCores = numCores
        self.CE_converter = CE_converter
        self.coeffs = [0,1,0]
        self.ms2_resolution = ms2_resolution
        self.uid = str(uuid.uuid1())
        tmpLib = self.uid + ".db"
        pkl.dump({"Positive":{},"Negative":{}},open(tmpLib, "wb"))
        #with nostdout():
        self.decID = DecoID(tmpLib, "reference", self.numCores, self.ms2_resolution)
        os.remove(tmpLib)

    def readMSMSData(self,msFile,targets,tic_cutoff,frag_cutoff):
        # make DecoID object


        # write temporary peak file
        targets.to_csv(self.uid + ".csv", index=False)

        # read in file and save all spectra
        self.decID.readData(msFile,self.ms2_resolution, True, True, self.ppm, peakDefinitions=self.uid + ".csv", tic_cutoff=tic_cutoff,
                       frag_cutoff=frag_cutoff)

        # structure to hold spectra
        output_dict = {}
        polarity = 0

        # process spectra for each CE
        if len(self.decID.samples) > 0:

            # get charge
            polarity = self.decID.samples[0]["mode"]
            switcher = {"Positive": 1, "Negative": -1}
            polarity = switcher[polarity]

            self.decID.samples = [x for x in self.decID.samples if targets.at[targets.index.values[x["group"]],"Charge"] == polarity]
            samplesAll = deepcopy(self.decID.samples)


            # get unique CEs
            ces = list(set([x["CE"] for x in samplesAll]))
            ces.sort()

            # iterate over CEs

            ceList = []
            args = []
            gs = []
            for ce in ces:

                # parse relevant samples
                self.decID.samples = [x for x in samplesAll if x["CE"] == ce]
                self.decID.label = str(ce)

                groups = list(set([x["group"] for x in self.decID.samples]))

                for g in groups:
                    specs = [x["spectra"] for x in self.decID.samples if x["group"] == g]
                    args.append([specs,self.ms2_resolution])
                    ceList.append(ce)
                    gs.append(g)

            if len(args) > 0:
                p = Pool(min([self.numCores,len(args)]))
                results = p.starmap(sumSpectra,args,chunksize=int(len(args)/min([self.numCores,len(args)])))
                p.close()
                p.join()

                spectra = {(g,ce):spec for g,spec,ce in zip(gs,results,ceList)}


                names = []
                ceList = []
                args = []
                for (ind,ce),spectrum in spectra.items():
                    row = targets.iloc[ind,:]
                    rts = [x["rt"] for x in self.decID.samples if x["group"] == ind]
                    args.append([spectrum,rts,self.decID.ms1,row["mz"],row["rt_start"],row["rt_end"],self.ppm])
                    names.append(row["Name"])
                    ceList.append(ce)

                p = Pool(min([self.numCores,len(args)]))
                results = p.starmap(normalizeSpectrum,args,chunksize=int(len(args)/min([self.numCores,len(args)])))
                p.close()
                p.join()

                for name,ce,spectrum in zip(names,ceList,results):
                    if name not in output_dict:
                        output_dict[name] = {}
                    output_dict[name][ce] = spectrum


        os.remove(self.uid + ".csv")

        return output_dict,polarity

    def getConversionEquationString(self):
        return "CE = " + str(self.coeffs[0]) + "*mz + " + str(self.coeffs[1]) + "*ce + " + str(self.coeffs[2])

    def buildConversion(self,msFiles,targetTransitions,tic_cutoff=0,frag_cutoff=0,frag_ppm_tolerance=10):
        #get unique cpd targets
        cpds = []
        toDrop = []
        for index, row in targetTransitions.iterrows():
            if (row['Name'],row["Charge"]) in cpds:
                toDrop.append(index)
            else:
                cpds.append((row['Name'],row["Charge"]))
        uniquePrecursors = targetTransitions.drop(toDrop)

        #read files and store spectra
        spectra = {}
        for msFile in msFiles:
            print("reading data...")
            #get in and process spectra
            specs,polarity = self.readMSMSData(msFile,uniquePrecursors,tic_cutoff,frag_cutoff)
            if polarity not in spectra:
                spectra[polarity] = {}
            for cpd in specs:
                if cpd not in spectra[polarity]: spectra[polarity][cpd] = {}
                for ce in specs[cpd]:
                    if cpd in spectra[polarity] and ce in spectra[polarity][cpd]:
                        spectra[polarity][cpd][ce] = sumSpectra([spectra[polarity][cpd][ce],specs[cpd][ce]],self.ms2_resolution)
                    else:
                        spectra[polarity][cpd][ce] = specs[cpd][ce]
            print("read data")

        print("Finding target transitions")
        #iterate through results and look for fragments
        optimals = {}
        for index,row in targetTransitions.iterrows():
            optimals[index] = {"HRMS_CE":-1}
            #if MS/MS collected
            if row["Charge"] in spectra and row["Name"] in spectra[row["Charge"]]:
                #gather CE intensity results for transition
                tmp = {}
                for ce,spec in spectra[row["Charge"]][row["Name"]].items():
                    for mz2,i in spec.items():
                        if 1e6 * np.abs(mz2-row["Product mz"])/row["Product mz"] < frag_ppm_tolerance:
                            if ce not in tmp:
                                tmp[ce] = 0
                            tmp[ce] += i
                #get best CE
                keys = list(tmp.keys())
                if len(keys) > 0:
                    keys.sort(key=lambda x: tmp[x],reverse=True)
                    optimals[index]["HRMS_CE"] = keys[0]

        print("gathering optimal CEs")
        #merge CE information
        merged = pd.concat((targetTransitions,pd.DataFrame.from_dict(optimals,orient="index")),axis=1)
        merged = merged[merged["HRMS_CE"] >= 0]

        #get CE pairs
        X = np.array([merged["mz"].values,merged["HRMS_CE"],[1 for _ in range(len(merged))]]).transpose()
        y = merged["CE"].values

        #find equation
        linreg = LinearRegression(fit_intercept=False)
        linreg.fit(X,y)
        self.CE_converter = lambda x: linreg.predict([[x[0],x[1],1]])[0]
        self.coeffs = linreg.coef_

        converted_info = {}
        for index,row in merged.iterrows():
            converted_info[index] = {"Converted_CE":self.CE_converter([row["mz"],row["HRMS_CE"]])}

        merged = pd.concat((merged,pd.DataFrame.from_dict(converted_info,orient="index")),axis=1,ignore_index=False)

        print("built conversion")
        return merged


    def createSRMsCE(self,msFile,targets,tic_cutoff=0,frag_cutoff=0,lowestProduct = 13.9,numProduct=5,outfile="none"):
        print("reading data...")
        output_dict,polarity = self.readMSMSData(msFile,targets,tic_cutoff,frag_cutoff)
        print("read data")

        cpds = list(output_dict.keys())

        if len(output_dict) > 0:
            #find best fragments
            srm = {}
            p = Pool(min([self.numCores,len(cpds)]))
            args = []
            for cpd in cpds:
                #get precursor mz
                precursorMz = targets[targets["Name"] == cpd]
                precursorMz = precursorMz[precursorMz["Charge"] == polarity]
                precursorMz = precursorMz.at[precursorMz.index.values[0],"mz"]
                args.append([output_dict[cpd],precursorMz,lowestProduct,numProduct])

            print("starting to find transitions")
            result = p.starmap(findTransitions,args)
            p.close()
            p.join()
            breakdown_curves = {}
            ind = 0
            print("writing results")
            for cpd,(breakdown_curve,cesOpt,topFrags) in zip(cpds,result):

                breakdown_curves[(cpd,polarity)] = breakdown_curve

                #addToSRM
                tmp = targets[(targets["Name"] == cpd) & (targets["Charge"] == polarity)]
                cols = tmp.columns.values
                i = tmp.index.values[0]

                for frag,ce in zip(topFrags,cesOpt):
                    srm[ind] = {c:tmp.at[i,c] for c in cols}
                    srm[ind]["Product mz"] = frag
                    srm[ind]["CE"] = self.CE_converter([srm[ind]["mz"],ce])
                    srm[ind]["Normalized Intensity"] = breakdown_curves[(cpd,polarity)][frag][ce]
                    srm[ind]["Charge"] = polarity
                    ind += 1

            srm = pd.DataFrame.from_dict(srm,orient="index")

        else:
            srm = pd.DataFrame()
            breakdown_curves = {}

        if outfile != "none":
            srm.to_csv(outfile)

        return srm,breakdown_curves


def normalizeSpectrum(spectrum,rts,ms1,mz,rt_start,rt_end,ppm):

    # compute interpolated MS1 intensity
    #func = extractChromatogram(row["mz"], self.decID.ms1, [row["rt_start"], row["rt_end"]], self.ppm)
    func = extractChromatogram(mz, ms1, [rt_start, rt_end], ppm)

    allRts = np.linspace(rt_start,rt_end,1000)
    totalPeakArea = np.trapz([func(x)[0] for x in allRts],allRts)

    normalizer = np.sum([func(x) for x in rts])

    # normalize spectra and save result
    return {key: totalPeakArea * val / normalizer for key, val in spectrum.items()}


def findTransitions(spectra,precursorMz,lowestProduct,numProduct):
    frag_int_dict = {}
    # iterate over spectra at different ces
    for ce in spectra:
        # collect fragment intensities and ces
        for mz, i in spectra[ce].items():
            if precursorMz - mz > lowestProduct:
                if mz not in frag_int_dict:
                    frag_int_dict[mz] = {}
                frag_int_dict[mz][ce] = i

    # sort fragments by highest observed intensity
    topFrags = list(frag_int_dict.keys())
    topFrags.sort(key=lambda x: np.max(list(frag_int_dict[x].values())), reverse=True)

    # take top fragments
    topFrags = topFrags[:numProduct]
    breakdown_curve = {mz: breakdown for mz, breakdown in frag_int_dict.items() if mz in topFrags}

    # get best ces
    cesOpt = []
    for frag in topFrags:
        tmp = list(breakdown_curve[frag].keys())
        tmp.sort(key=lambda x: breakdown_curve[frag][x], reverse=True)
        cesOpt.append(tmp[0])

    return breakdown_curve,cesOpt,topFrags


def sumSpectra(specs,tol):
    compositeSpectrum = {}
    for s in specs:
        for mz, i in s.items():
            found = False
            for mz2 in compositeSpectrum:
                if np.abs(mz - mz2) < (10 ** (-1 * tol)) / 2:
                    compositeSpectrum[mz2] += i
                    found = True
                    break
            if not found:
                compositeSpectrum[mz] = i
    return compositeSpectrum


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


def integrateTargets(ms1,targets,ppm):
    intensities = {}
    for index,row in targets:
        func = extractChromatogram(row["mz"],ms1,[row["rt_start"],row["rt_end"]],ppm)
        allRts = np.linspace(row["rt_start"],row["rt_end"],1000)
        totalPeakArea = np.trapz([func(x) for x in allRts],allRts)
        intensities[index] = {"Peak Area": totalPeakArea}

    intensities = pd.DataFrame.from_dict(intensities)
    return pd.concat((targets,intensities),axis=1,ignore_index=False)

