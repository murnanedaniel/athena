# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

# this script can be used to check the output ntuples of NSWPRDValAlg

from __future__ import print_function
import os, sys, ROOT, argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='checkNSWValTree', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inputFile', help='choose input ROOT file', default='NSWPRDValAlg.digi.ntuple.root', type=str)
    parser.add_argument('--checkPRD', help='check also the RPDs (for reco validation only)', action='store_true', default=False)
    parser.add_argument('--checkHits', help='check also the Hits (Hits are skimmed out after Digitisation step!)', action='store_true', default=False)
    parser.add_argument('--checkDigits', help='check also the Digits (Hits are skimmed out after Digitisation step!)', action='store_true', default=False)
    parser.add_argument('--checkSDO', help='check also the SDO (Hits are skimmed out after Digitisation step!)', action='store_true', default=False)
    Options = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    if not os.path.exists(Options.inputFile):
        print ('ERROR: File %s does not exist'%Options.inputFile)
        sys.exit(1)

    inputFile = ROOT.TFile(Options.inputFile, "READ")
    if not inputFile:
        print ('ERROR: Failed to open file %s'%Options.inputFile)
        sys.exit(1)
    inputTree = inputFile.Get("NSWValTree")
    if not inputTree:
        print ('ERROR: NSWValTree does not exist in file %s'%Options.inputFile)
        sys.exit(1)

    nEntries = inputTree.GetEntries()
    if nEntries==0:
        print ('ERROR: NSWValTree of file %s has 0 entries'%Options.inputFile)
        sys.exit(1)

    nHitsMM = 0
    nHitsSTGC = 0
    nDigitsMM = 0
    nDigitsSTGC = 0
    nSDOMM = 0
    nSDOSTGC = 0
    nRDOMM = 0
    nRDOSTGC = 0
    nPRDMM = 0
    nPRDSTGC = 0
    for i in range(nEntries):
        inputTree.GetEntry(i)
        nRDOMM       += inputTree.N_RDO_MM
        nRDOSTGC     += inputTree.N_RDO_sTGC

        if Options.checkPRD:
            nPRDMM += inputTree.PRD_MM
            nPRDSTGC += inputTree.PRD_sTGC

        if Options.checkHits:
            nHitsMM      += inputTree.Hits_MM_nHits
            nHitsSTGC    += inputTree.Hits_sTGC_nHits

        if Options.checkHits:
            nDigitsMM    += inputTree.N_Digits_MM
            nDigitsSTGC  += inputTree.N_Digits_sTGC

        if Options.checkSDO:
            nSDOMM       += inputTree.nSDO_MM
            nSDOSTGC     += inputTree.nSDO_sTGC

    if Options.checkHits and nHitsMM==0:
        print ('ERROR: NSWValTree of file %s has 0 MM Hits'%Options.inputFile)
        sys.exit(1)
    elif Options.checkHits and  nHitsSTGC==0:
        print ('ERROR: NSWValTree of file %s has 0 STGC Hits'%Options.inputFile)
        sys.exit(1)
    elif Options.checkDigits and nDigitsMM==0:
        print ('ERROR: NSWValTree of file %s has 0 MM Digits'%Options.inputFile)
        sys.exit(1)
    elif Options.checkDigits and nDigitsSTGC==0:
        print ('ERROR: NSWValTree of file %s has 0 STGC Digits'%Options.inputFile)
        sys.exit(1)
    elif Options.checkSDO and  nSDOMM==0:
        print ('ERROR: NSWValTree of file %s has 0 MM SDOs'%Options.inputFile)
        sys.exit(1)
    elif Options.checkSDO and  nSDOSTGC==0:
        print ('ERROR: NSWValTree of file %s has 0 STGC SDOs'%Options.inputFile)
        sys.exit(1)
    elif nRDOMM==0:
        print ('ERROR: NSWValTree of file %s has 0 MM RDOs'%Options.inputFile)
        sys.exit(1)
    elif nRDOSTGC==0:
        print ('ERROR: NSWValTree of file %s has 0 STGC RDOs'%Options.inputFile)
        sys.exit(1)
    elif Options.checkPRD and nPRDMM==0:
        print ('ERROR: NSWValTree of file %s has 0 MM PRDs'%Options.inputFile)
        sys.exit(1)
    elif Options.checkPRD and nPRDSTGC==0:
        print ('ERROR: NSWValTree of file %s has 0 STGC PRDs'%Options.inputFile)
        sys.exit(1)

    print ('INFO: All fine with file %s'%Options.inputFile)

    if Options.checkHits: 
        print ('INFO: Number of found MM hits:\t\t%i'%nHitsMM)
        print ('INFO: Number of found STGC hits:\t%i'%nHitsSTGC)

    if Options.checkDigits: 
        print ('INFO: Number of found MM digits:\t%i'%nDigitsMM)
        print ('INFO: Number of found STGC digits:\t%i'%nDigitsSTGC)

    if Options.checkSDO: 
        print ('INFO: Number of found MM SDOs:\t\t%i'%nSDOMM)
        print ('INFO: Number of found STGC SDOs:\t%i'%nSDOSTGC)

    print ('INFO: Number of found MM RDOs:\t\t%i'%nRDOMM)
    print ('INFO: Number of found STGC RDOs:\t%i'%nRDOSTGC)

    if Options.checkPRD:
        print ('INFO: Number of found MM PRDs:\t\t%i'%nPRDMM)
        print ('INFO: Number of found STGC PRDs:\t%i'%nPRDSTGC)

