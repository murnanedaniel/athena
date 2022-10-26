#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''
@file SetupEgammaMonitoring.py
@author B. Laforge
4 May 2020
@brief Example python configuration for the Run III AthenaMonitoring package
'''
from AthenaConfiguration.ComponentFactory import CompFactory
import TrigEgammaMonitoring.TrigEgammaMonitCategoryMT as egammaConf

def BookHistogramsPerRegions(thegroupe,theparttype,thename,title,path,xbins,xmin,xmax,thetype="TH1F",thecut="is_pt_gt_4gev"):
    '''
    Function to configure the egamma histograms per region (BARREL, CRACK, ENDCAP, FORWARD)
    The list is different for electrons, forward electrons and photons
    '''
    Regions = ['BARREL', 'CRACK', 'ENDCAP']
    for n in range(len(Regions)):
        hname = thename + "in" + Regions[n] #+ ";" + thename + "in" + Regions[n] #+ " " + thepartype + thereconame + thewithTrigger
        htitle = title + " " + theparttype + " " + Regions[n]
        thefinalcut = thecut+Regions[n] 
        thegroupe.defineHistogram(hname,title=htitle, path=path, xbins=xbins,xmin=xmin,xmax=xmax,type = thetype, cutmask = thefinalcut)

def BookHistogramsPerForwardRegions(thegroupe,theparttype,thename,title,path,xbins,xmin,xmax,thetype="TH1F",thecut="is_pt_gt_10gev"):
    '''
    Function to configure the egamma histograms per region (BARREL, CRACK, ENDCAP, FORWARD)
    The list is different for electrons, forward electrons and photons
    '''
    Regions = ['ENDCAP','FORWARD']
    for n in range(len(Regions)):
        hname = thename + "in" + Regions[n] # + ";" + thename + "in" + Regions[n] #+ " " + thepartype + thereconame + thewithTrigger
        htitle = title + " " + theparttype + " " + Regions[n]
        thefinalcut = thecut+Regions[n] 
        thegroupe.defineHistogram(hname,title=htitle, path=path, xbins=xbins,xmin=xmin,xmax=xmax,type = thetype, cutmask = thefinalcut)

def BookHistograms(groupe,reconame,particletype,withTrigger=""):
    '''
    Function to configure the egamma histograms
    The list is different for electrons and photons
    usage example : BookHistograms(GroupElectronLHTight,"LHTight","Electron","WithTrigger")
    '''

    if particletype in ["Electron","Photon"] :
        if particletype in ["Electron"] :
            prefix = "electron"
            tlabel = "_{e}"
        if particletype in ["Photon"] :
            prefix = "photon"
            tlabel = "_{#gamma}"

        hname = "N"
        htitle = "Number of " + prefix + "s (" + reconame + " " + withTrigger+ ") ; N" + tlabel + " ; N_{event}"
        groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

        if particletype in ["Photon"] :
            hname = "NConv"
            htitle = "Number of converted photons (" + reconame + " " + withTrigger+ ") ; N_{Conv. #gamma} ; N_{event}"
            groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

            hname = "NUnconv"
            htitle = "Number of unconverted photons (" + reconame + " " + withTrigger+ ") ; N_{Unconv. #gamma} ; N_{event}"
            groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

        hname= "Et"
        htitle= particletype + " transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ")" + " ; Et"+ tlabel +" ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=2.,xmax=102000.0)

        if particletype in ["Photon"] :

            hname= "EtConv"
            htitle= "Converted #gamma transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ")" + " ; Et" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=2.,xmax=102000.0, cutmask = 'is_pt_gt_4gevandconv')

            hname= "EtUnconv"
            htitle= "Unconverted #gamma transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ")" + " ; Et" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=2.,xmax=102000.0, cutmask = 'is_pt_gt_4gevandunconv')

        hname= "Eta"
        htitle= particletype + " #eta" + " (" + reconame + " " + withTrigger + ")"+ " ; #eta" + tlabel + " ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',type = "TH1D", xbins=64,xmin=-3.2,xmax=3.2)

        if particletype in ["Photon"] :
            hname= "EtaConv"
            htitle= "Converted photons #eta" + " (" + reconame + " " + withTrigger + ")"+ " ; #eta" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2,cutmask = 'is_pt_gt_4gevandconv')

            hname= "EtaUnconv"
            htitle= "Unconverted photons #eta" + " (" + reconame + " " + withTrigger + ")"+ " ; #eta" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2,cutmask = 'is_pt_gt_4gevandunconv')

        hname= "Phi"
        htitle= particletype + " #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #phi" + tlabel + " ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2)

        if particletype in ["Photon"] :
            hname= "PhiConv"
            htitle= "Converted photon #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #phi" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2,cutmask = 'is_pt_gt_4gevandconv')

            hname= "PhiUnconv"
            htitle= "Unconverted photon #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #phi" + tlabel + " ; N" + tlabel
            groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2, cutmask = 'is_pt_gt_4gevandunconv')

        hname= "Eta,Phi;Eta_Phi_with_Pt_gt_2.5GeV"
        htitle= particletype + " #eta,#phi map (candidates with Pt>2.5GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
        groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask = 'is_pt_gt_2_5gev')

        if particletype in ["Photon"] :
            hname= "Eta,Phi;Eta_Phi_Conv_with_Pt_gt_2.5GeV"
            htitle= "(#eta,#phi) map of Conv. #gamma with Pt>2.5GeV" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2, cutmask = 'is_pt_gt_2_5gevandconv')

            hname= "Eta,Phi;Eta_Phi_Unconv_with_Pt.gt.2.5GeV"
            htitle= "(#eta,#phi) map of Unconv. #gamma with Pt>2.5GeV" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2, cutmask = 'is_pt_gt_2_5gevandunconv')

        hname= "Eta,Phi;Eta_Phi_distribution_with_Pt.gt.4GeV"
        htitle= particletype + " #eta,#phi map (candidates with Pt>4GeV)" + " (" + reconame + " " + prefix + " " + withTrigger + ") ; #eta ; #phi"
        groupe.defineHistogram(hname,title=htitle,path='',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_4gev')

        if particletype in ["Photon"] :
            hname= "Eta,Phi;Eta_Phi_Conv_with_Pt.gt.4GeV"
            htitle= particletype + " #eta,#phi map (candidates with Pt>4GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_4gevandconv')

            hname= "Eta,Phi;Eta_Phi_Unconv_with_Pt.gt.4GeV"
            htitle= particletype + " #eta,#phi map (candidates with Pt>4GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_4gevandunconv')

        hname= "Eta,Phi;Eta_Phi_with_Pt.gt.20GeV"
        htitle= particletype + " #eta,#phi map (candidates with Pt>20GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi "
        groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_20gev')

        if particletype in ["Photon"] :
            hname= "Eta,Phi;Eta_Phi_distribution_of_Conv._photons_with_Pt.gt.20GeV"
            htitle= particletype + " #eta,#phi map (candidates with Pt>20GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_20gevandconv')

            hname= "Eta,Phi;Eta_Phi_distribution_of_Unconv._photons_with_Pt.gt.20GeV"
            htitle= particletype + " #eta,#phi map (candidates with Pt>20GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
            groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_20gevandunconv')

        hname= "TopoEtCone40"
        htitle = particletype + " Topocluster Isolation Energy" + " (" + reconame + " " + withTrigger + ") ; TopoEtCone40 ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=-10000.,xmax=40000.)

        hname= "PtCone20"
        htitle = particletype + " Track Isolation Pt" + " (" + reconame + " " + withTrigger + ") ; PtCone20 ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-10000.,xmax=40000.)

        hname= "Time"
        htitle = "Time associated with " + particletype + " cluster [ns]" + " (" + reconame + " " + withTrigger + ") ; Time [ns] ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle, path='',xbins=90,xmin=-30.,xmax=60.)

        if particletype in ["Photon"] :
            hname= "RConv"
            htitle = "Photon Conversion radius [mm] " + " (" + reconame + " " + withTrigger + ") ; Rconv ; N_{#gamma}"
            groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=0.,xmax=800.)

        # histograms per region
        hname= "Et"
        htitle= particletype + " transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ") ; Et" + tlabel + " ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=100,xmin=2.,xmax=102000.0,thecut='is_pt_gt_2_5gev')

        hname= "Eta"
        htitle= particletype + " #eta" + " (" + reconame + " " + withTrigger + ")" + " ; #eta" + tlabel + " ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=64,xmin=-3.2,xmax=3.2)

        hname= "Phi"
        htitle= particletype + " #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #eta" + tlabel + " ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=64,xmin=-3.2,xmax=3.2)

        hname= "Time"
        htitle = "Time associated with " + particletype + " cluster [ns]" + " (" + reconame + " " + withTrigger + ") ; Time [ns] ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=90,xmin=-30.,xmax=60.)

        if particletype in ["Photon"] :
            hname= "RConv"
            htitle = "Photon Conversion radius [mm] " + " (" + reconame + " " + withTrigger + ") ; Rconv [mm] ; N_{#gamma} "
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=100,xmin=0.,xmax=800.)

            hname= "ConvType"
            htitle = "Photon Conversion Type " + " (" + reconame + " " + withTrigger + ") ; ConvType ; N_{#gamma}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=6,xmin=-0.5,xmax=5.5)

            hname= "ConvTrkMatch1"
            htitle = "Photon Conversion Trak1 Match " + " (" + reconame + " " + withTrigger + ") ; ConTrkMatch1 ; N_{#gamma}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=4,xmin=-0.5,xmax=0.)

            hname= "ConvTrkMatch2"
            htitle = "Photon Conversion Trak1 Match " + " (" + reconame + " " + withTrigger + ") ; ConTrkMatch2 ; N_{#gamma}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=4,xmin=0.,xmax=4.)

        # ID Panel

        hname= "Ehad1"
        htitle = particletype + " energy leakage in 1st hadronic sampling " + " (" + reconame + " " + withTrigger + ")" + "; Ehad1 (MeV) ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-1000.,xmax=10000.)

        if particletype in ["Electron"] :
            hname= "EoverP"
            htitle = particletype + " matched track E over P " + " (" + reconame + " " + withTrigger + ")" + "; EoverP ; N" + tlabel
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=5.)

        hname= "CoreEM"
        htitle = particletype + " core energy in EM calorimeter " + " (" + reconame + " " + withTrigger + ")" + "; E (MeV) ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-5000.,xmax=250000.)

        hname= "F0"
        htitle = particletype + " fractional energy in Presampler " + " (" + reconame + " " + withTrigger + ")" + "; F0 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

        hname= "F1"
        htitle = particletype + " fractional energy in 1st sampling " + " (" + reconame + " " + withTrigger + ")" + "; F1 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

        hname= "F2"
        htitle = particletype + " fractional energy in 2nd sampling " + " (" + reconame + " " + withTrigger + ")" + "; F2 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

        hname= "F3"
        htitle = particletype + " fractional energy in 3rd sampling " + " (" + reconame + " " + withTrigger + ")" + "; F3 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

        hname= "Re233e237"
        htitle = particletype + "  uncor. energy fraction in 3x3/3x7 cells in em sampling 2 " + " (" + reconame + " " + withTrigger + ")" + "; R 3x3/3x7 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=2.)

        hname= "Re237e277"
        htitle = particletype + "  uncor. energy fraction in 3x7/7x7 cells in em sampling 2 " + " (" + reconame + " " + withTrigger + ")" + "; R 3x7/7x7 ; N" + tlabel
        BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=2.)

        # Specific plots for electrons (related to electron track) in Expert/Tracks Panel

        if particletype in ["Electron"] :
            hname= "NOfBLayerHits"
            htitle = particletype + " NOfBLayerHits (" + reconame + " " + withTrigger + ") ; N_{BlayerHits} ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=6,xmin=-0.5,xmax=5.5)

            hname= "NOfPixelHits"
            htitle = particletype + " NOfPixelHits (" + reconame + " " + withTrigger + ") ; N_{PixelHits} ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=6,xmin=-0.5,xmax=5.5)

            hname= "NOfSCTHits"
            htitle = particletype + " NOfSCTHits (" + reconame + " " + withTrigger + ") ; N_{SCTHits} ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=25.5)

            hname= "NOfTRTHits"
            htitle = particletype + " NOfTRTHits (" + reconame + " " + withTrigger + ") ; N_{TRTHits} ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=50.5)

            hname= "NOfTRTHighThresholdHits"
            htitle = particletype + " NOfTRTHighThresholdHits (" + reconame + " " + withTrigger + ") ; N_{TRT HighThres. Hits} ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=50.5)

            hname= "DeltaEta1"
            htitle = particletype + " track match #Delta #eta (1st sampling) " + " (" + reconame + " " + withTrigger + ")  ; #Delta #eta ; N_{e} "
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=50,xmin=-0.05,xmax=0.05)

            hname= "DeltaPhi2"
            htitle = particletype + " track match #Delta #Phi (2st sampling) " + " (" + reconame + " " + withTrigger + ")  ; #Delta #phi ; N_{e} "
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=50,xmin=-0.15,xmax=0.15)

            hname= "Trackd0"
            htitle = particletype + " track d0 " + " (" + reconame + " " + withTrigger + ")  ; d0 ; N_{e}"
            BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=100,xmin=-5.,xmax=5.)

        # LumiBlock dependant histogram

        hname= "LBEvoN"
        htitle= "Number of " + particletype + "s per LB (" + reconame + " " + withTrigger + ") ; LB ; N" + tlabel
        groupe.defineHistogram(hname,title=htitle,path='byLB',type='TH1F',xbins=2000,xmin=-0.5,xmax=1999.5)

        if particletype in ["Photon"] :
            hname= "LBEvoNPhotonsUnconv;Number_of_Unconv._photons_vs_LB"
            htitle = "Unconverted photon number versus LB " + " (" + reconame + " " + withTrigger + ") ; LB ; N_{Unconv. #gamma}"
            groupe.defineHistogram(hname,title=htitle, path='byLB',xbins=2000,xmin=-0.5,xmax=1999.)

            hname= "LBEvoNPhotonsConv;Number_of_Conv._photons_vs_LB"
            htitle = "Converted photon number versus LB " + " (" + reconame + " " + withTrigger + ") ; LB ; N_{Conv. #gamma}"
            groupe.defineHistogram(hname,title=htitle, path='byLB',xbins=2000,xmin=-0.5,xmax=1999.)

            hname= "is_pt_gt_4gevandconv,LB;Conversion_fraction_vs_LB"
            htitle = "Converted photon fraction versus LB " + " (" + reconame + " " + withTrigger + ") ; LB ; Conv. #gamma fraction"
            groupe.defineHistogram(hname,title=htitle, path='byLB',type = 'TEfficiency', xbins=2000,xmin=-0.5,xmax=1999.)

    return True

def BookTnPHistograms(groupe,reconame,TnPType,MassLowerCut,MassUpperCut):
    '''
    Function to configure the egamma histograms for a Tag and Probe resonance decayong to ee.
    usage example : BookHistogramsTnP(GroupElectronLHTight,"LHTight","Electron","WithTrigger")
    '''

    tlabel = "_{" + TnPType + "}"
    particletype = 'Electron'

    hname = "LB;Number_of_" + TnPType + "_candidates_vs_LB"
    htitle = "Number of " + TnPType +" candidates vs LB ; LB ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=2000,xmin=-0.5,xmax=1999.5)

    hname = "MassZ;" + TnPType + "_candidate_mass"
    htitle = TnPType + " candidate mass ; M_{ee} [MeV/ c^{2}]; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=100,xmin=MassLowerCut,xmax=MassUpperCut)

    hname = "MassZ_BARREL;" + TnPType + "_candidate_mass_distribution_Lead_in_Barrel"
    htitle = TnPType + " candidate mass in Barrel ; M_{ee} [MeV/ c^{2}]; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=100,xmin=MassLowerCut,xmax=MassUpperCut)

    hname = "MassZ_ENDCAP;" + TnPType + "_candidate_mass_distribution_Lead_in_EndCap"
    htitle = TnPType + " candidate mass in EndCap ; M_{ee} [MeV/ c^{2}]; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=100,xmin=MassLowerCut,xmax=MassUpperCut)

    hname = "MassZ_CRACK;" + TnPType + "_candidate_mass_distribution_Lead_in_Crack"
    htitle = TnPType + " candidate mass in Crack Regionl ; M_{ee} [MeV/ c^{2}]; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=100,xmin=MassLowerCut,xmax=MassUpperCut)

    # The next two histos receive 1 entry per leading electron probe candidate

    hname = "EtaZ;Number_of_"+ TnPType + "_candidates_vs_eta_leading_e"
    htitle = "Number of "+ TnPType +" candidates vs #eta of leading e" + " ; #eta_{e_{lead}} ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=64,xmin=-3.2,xmax=3.2)

    hname = "PhiZ;Number_of_"+ TnPType + "_candidates_vs_phi_leading_e"
    htitle = "Number of "+ TnPType +" candidates vs #phi of e" + " ; #phi_{e_{lead}} ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='TnPCandidate',xbins=64,xmin=-3.2,xmax=3.2)

    # EFFICIENCIES IN EFFICIENCIES PANEL

    hname = "is_Tight,Etprobe;TnP_ID_efficiency_vs_etprobe"
    htitle = "Tag & Probe  ID efficiency vs Et_{e} [MeV] ; Et_{e} [MeV] ; LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=100,xmin=0.,xmax=250000.)

    hname = "is_Tight,Etaprobe;TnP_ID_efficiency_vs_etaprobe"
    htitle = "Tag & Probe  ID efficiency vs #eta ; #eta _{e} ; LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "is_Tight,Phiprobe;TnP_ID_efficiency_vs_phiprobe"
    htitle = "Tag & Probe  ID efficiency vs #eta ; #phi _{e} ; LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "is_Iso,Etprobe;TnP_Calo_Iso_efficiency_vs_etprobe"
    htitle = "Tag & Probe  Calo. Isolation efficiency vs Et _{e} [MeV] ; Et_{e} [MeV] ; Calo Iso Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=100,xmin=0.,xmax=250000.)

    hname = "is_Iso,Etaprobe;TnP_Calo_Iso_efficiency_vs_etaprobe"
    htitle = "Tag & Probe  Calo. Isolation efficiency vs #eta_{e} ; #eta_{e} ; Calo. Iso Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "is_Iso,Phiprobe;TnP_Calo_Iso_efficiency_vs_phiprobe"
    htitle = "Tag & Probe Calo. Isolation efficiency vs #phi_{e} ; #phi_{e}; Calo Iso. Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "is_IsoandTight,Etprobe;TnP_Calo_Iso_and_LHTight_efficiency_vs_etprobe"
    htitle = "Tag & Probe (Isolation & LHTight ID) efficiency vs Et_{e} [MeV] ; Et_{e} [MeV] ; Calo Iso & LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=100,xmin=0.,xmax=250000.)

    hname = "is_IsoandTight,Etaprobe;TnP_Calo_Iso_and_LHTight_efficiency_vs_etaprobe"
    htitle = "Tag & Probe (Isolation & LHTight ID) efficiency vs #eta ; #eta _{e} ; Calo. Iso & LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "is_IsoandTight,Phiprobe;TnP_Calo_Iso_and_LHTight_ID_efficiency_vs_phiprobe"
    htitle = "Tag & Probe (Isolation & LHTight ID) efficiency vs #phi_{e} ; #phi_{e} ; Calo Iso & LHTight ID Eff. "
    groupe.defineHistogram(hname,title=htitle, path='Efficiencies',type = 'TEfficiency', xbins=64,xmin=-3.2,xmax=3.2)


    # Probe electron candidate distributions

    hname = "Etprobe;Etprobe_distribution"
    htitle = "Number of "+ TnPType +" candidates vs of leading e ; Et_{probe} ; N_{e_{probe}}"
    groupe.defineHistogram(hname,title=htitle, path='ElectronProbes',xbins=100,xmin=-1000.0,xmax=200000.0)

    hname = "Etaprobe;Etaprobe_distribution"
    htitle = "#eta distribution of probe candidates ; #eta_{e_{probe}} ; N_{e_{probe}}"
    groupe.defineHistogram(hname,title=htitle,path='ElectronProbes', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "Phiprobe;Phiprobe_distribution"
    htitle = "#phi distribution of probe candidates ; #phi_{e_{probe}} ; N_{e_{probe}}"
    groupe.defineHistogram(hname,title=htitle,path='ElectronProbes', xbins=64,xmin=-3.2,xmax=3.2)

    hname = "Etaprobe,Phiprobe;Eta_Phi_map_of_electron_probes_for_T_and_P_"+ TnPType
    htitle= "Concerted photons #eta,#phi map (candidates with Pt>4GeV)" + " ; #eta_{e} ; #phi_{e}"
    groupe.defineHistogram(hname,title=htitle,path='ElectronProbes',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2)

    hname = "Timeprobe;Time_of_electron_probe"
    htitle = "Time distribution of probe candidates ; Time_{e_{probe}} [ns] ; N_{e_{probe}}"
    groupe.defineHistogram(hname,title=htitle,path='ElectronProbes',xbins=90,xmin=-30.,xmax=60.)

    hname = "TopoEtCone40probe;TopoEtCone40_of_electron_probe"
    htitle = "Electron probe Topocluster Isolation Energy ; TopoEtCone40 [MeV] ; N_{e_{probe}} "
    groupe.defineHistogram(hname,title=htitle, path='ElectronProbes',xbins=100,xmin=-10000.,xmax=40000.)

    hname = "PtCone20probe;PtCone20_of_electron_probe"
    htitle = "Electron probe Track Isolation Pt [MeV]; PtCone20 [MeV] ; N_{e_{probe}} "
    groupe.defineHistogram(hname,title=htitle, path='ElectronProbes',xbins=64,xmin=-10000.,xmax=40000.)

    # Track distributions

    hname = "NOfBLayerHitsProbe"
    htitle = "Electron probe NOfBLayerHits ; N_{BlayerHits} ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=6,xmin=-0.5,xmax=5.5)

    hname = "NOfPixelHitsProbe"
    htitle = "Electron probe NOfPixelHits ; N_{PixelHits} ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=6,xmin=-0.5,xmax=5.5)

    hname = "NOfSCTHitsProbe"
    htitle = "Electron probe NOfSCTHits ; N_{SCTHits} ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=25.5)

    hname = "NOfTRTHitsProbe"
    htitle = "Electron probe NOfTRTHits ; N_{TRTHits} ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=50.5)

    hname = "NOfTRTHighThresholdHitsProbe"
    htitle = "Electron probe NOfTRTHighThresholdHits ; N_{TRT HighThres. Hits} ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=26,xmin=-0.5,xmax=50.5)

    hname = "DeltaEta1Probe"
    htitle = "Electron probe track match #Delta #eta (1st sampling) ; #Delta #eta ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=50,xmin=-0.05,xmax=0.05)

    hname = "DeltaPhi2Probe"
    htitle = "Electron probe track match #Delta #Phi (2st sampling) ; #Delta #phi ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=50,xmin=-0.15,xmax=0.15)

    hname = "Trackd0Probe"
    htitle = "Electron probe track d0 ; d0 ; N_{e_{probe}}"
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Tracks',xbins=100,xmin=-5.,xmax=5.)

    # ID distributions

    hname = "Ehad1Probe"
    htitle = "Electron probe energy leakage in 1st hadronic sampling ; Ehad1 (MeV) ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-1000.,xmax=10000.)

    hname = "EoverPProbe"
    htitle = "Electron probe matched track E over P ; EoverP ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=5.)

    hname = "CoreEMProbe"
    htitle = "Electron probe core energy in EM calorimeter ; E (MeV) ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-5000.,xmax=250000.)

    hname = "F0Probe"
    htitle = "Electron probe fractional energy in Presampler ; F0 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

    hname = "F1Probe"
    htitle = "Electron probe fractional energy in 1st sampling ; F1 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

    hname = "F2Probe"
    htitle = "Electron probe fractional energy in 2nd sampling ; F2 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

    hname = "F3Probe"
    htitle = "Electron probe fractional energy in 3rd sampling ; F3 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=-0.2,xmax=1.)

    hname = "Re233e237Probe"
    htitle = "Electron probe  uncor. energy fraction in 3x3/3x7 cells in em sampling ; R 3x3/3x7 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=2.)

    hname = "Re237e277Probe"
    htitle = "Electron probe  uncor. energy fraction in 3x7/7x7 cells in em sampling 2 ; R 3x7/7x7 ; N_{e_{probe}} "
    BookHistogramsPerRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='ID',xbins=50,xmin=0.,xmax=2.)

    hname = "LBEvoNProbe"
    htitle= "Number of " + TnPType + " electron probe per LB ; LB ; N_{e_{probe}} "
    groupe.defineHistogram(hname,title=htitle,path='byLB',type='TH1F',xbins=2000,xmin=-0.5,xmax=1999.5)

    return True


def BookFwdElectronHistograms(groupe,reconame,particletype,withTrigger=""):
    '''
    Function to configure the forward electrons histograms
    usage example : BookHistograms(GroupElectronLHTight,"LHTight","Electron","WithTrigger")
    '''

    prefix = "fwdElectron"
    tlabel = "_{fwd. electron}"

    hname = "N"
    htitle = "Number of " + prefix + "s (" + reconame + " " + withTrigger+ ") ; N" + tlabel + " ; N_{event}"
    groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

    hname = "NinENDCAP"
    htitle = "Number of " + prefix + "s in ENDCAP (" + reconame + " " + withTrigger+ ") ; N" + tlabel + " ; N_{event}"
    groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

    hname = "NinFORWARD"
    htitle = "Number of " + prefix + "s in FORWARD (" + reconame + " " + withTrigger+ ") ; N" + tlabel + " ; N_{event}"
    groupe.defineHistogram(hname,title=htitle, path='',xbins=20,xmin=-0.5,xmax=19.5)

    hname= "Et" 
    htitle= particletype + " transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ")" + " ; Et"+ tlabel +" ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=2.,xmax=102000.0)

    hname= "Eta"
    htitle= particletype + " #eta" + " (" + reconame + " " + withTrigger + ")"+ " ; #eta" + tlabel + " ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2)

    hname= "Phi"
    htitle= particletype + " #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #phi" + tlabel + " ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=64,xmin=-3.2,xmax=3.2)

    hname= "Eta,Phi;Eta_Phi_distribution_Pt_gt_2.5GeV"
    htitle= particletype + " #eta,#phi map (candidates with Pt>2.5GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi"
    groupe.defineHistogram(hname,title=htitle,path='Expert',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask = 'is_pt_gt_2_5gev')

    hname= "Eta,Phi;Eta_Phi_distributionf_Pt_gt_10GeV"
    htitle= particletype + " #eta,#phi map (candidates with Pt>10GeV)" + " (" + reconame + " " + withTrigger + ") ; #eta ; #phi "
    groupe.defineHistogram(hname,title=htitle,path='',type='TH2F',xbins=64,xmin=-3.2,xmax=3.2,ybins=64,ymin=-3.2,ymax=3.2,cutmask='is_pt_gt_10gev')

    # shower variables

    hname= "EnergyDensity"
    htitle= "Fwd electron 1st Moment Energy Density" + " (" + reconame + " " + withTrigger + ")" + " ; firstENGdensity ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=200,xmin=0.,xmax=2.0)

    hname= "FracMax"
    htitle= "Fwd electron fraction of most energetic cell distribution" + " (" + reconame + " " + withTrigger + ")" + " ; lateral moment ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=50,xmin=0.,xmax=1.0)

    hname= "Lateral"
    htitle= "Fwd electron lateral moment distribution" + " (" + reconame + " " + withTrigger + ")" + " ; lateral moment ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=0.,xmax=1.0)

    hname= "Longitudinal"
    htitle= "Fwd electron longitudinal moment distribution" + " (" + reconame + " " + withTrigger + ")" + " ; longitudinal moment ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=100,xmin=0.,xmax=1.0)

    hname= "SecondLambda"
    htitle= "Fwd electron lambda second moment distribution" + " (" + reconame + " " + withTrigger + ")" + " ; Second#Lambda ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=500,xmin=0.,xmax=10000.0)

    hname= "SecondR"
    htitle= "Fwd electron SecondR lateral moment distribution" + " (" + reconame + " " + withTrigger + ")" + " ; SecondR ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=500,xmin=0.,xmax=20000.0)

    hname= "CenterLambda"
    htitle= "Fwd Electron shower center from calo front face distribution" + " (" + reconame + " " + withTrigger + ")" + " ; Center lambda ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=500,xmin=0.,xmax=2000.0)

    hname= "Time"
    htitle = "Time associated with " + particletype + " cluster [ns]" + " (" + reconame + " " + withTrigger + ") ; Time [ns] ; N" + tlabel
    groupe.defineHistogram(hname,title=htitle, path='',xbins=90,xmin=-30.,xmax=60.)

    # info per forward region

    # histograms per region
    hname= "Et"
    htitle= particletype + " transverse energy [MeV]" + " (" + reconame + " " + withTrigger + ") ; Et" + tlabel + " ; N" + tlabel
    BookHistogramsPerForwardRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=100,xmin=2.,xmax=102000.0,thecut='is_pt_gt_2_5gev')

    hname= "Eta"
    htitle= particletype + " #eta" + " (" + reconame + " " + withTrigger + ")" + " ; #eta" + tlabel + " ; N" + tlabel
    BookHistogramsPerForwardRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=64,xmin=-3.2,xmax=3.2)

    hname= "Phi"
    htitle= particletype + " #phi" + " (" + reconame + " " + withTrigger + ")" + " ; #eta" + tlabel + " ; N" + tlabel
    BookHistogramsPerForwardRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=64,xmin=-3.2,xmax=3.2)

    hname= "Time"
    htitle = "Time associated with " + particletype + " cluster [ns]" + " (" + reconame + " " + withTrigger + ") ; Time [ns] ; N" + tlabel
    BookHistogramsPerForwardRegions(thegroupe = groupe, theparttype = particletype, thename = hname, title=htitle, path='Expert',xbins=90,xmin=-30.,xmax=60.)

    return True

def MonitorElectronConfig(inputFlags):
    '''Function to configures some algorithms in the monitoring system.'''

    ### STEP 1 ###

    from AthenaMonitoring import AthMonitorCfgHelper
    helper = AthMonitorCfgHelper(inputFlags,'ElectronAthMonitorCfg')

    ### STEP 2 ###

    SpareElectronMonitoringGroups = [ ]

    MonitorElectronAlgorithm=CompFactory.MonitorElectronAlgorithm
    elLHTightMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elLHTightMonAlg')
    elLHTightTrigMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elLHTightTrigMonAlg')

    elLHLooseMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elLHLooseMonAlg')
    elLHLooseTrigMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elLHLooseTrigMonAlg')

    if ("CBTight" in SpareElectronMonitoringGroups) :
        elCBTightMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elCBTightMonAlg')
    if ("CBTightTrig" in SpareElectronMonitoringGroups) :
        elCBTightTrigMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elCBTightTrigMonAlg')
    if ("CBLoose" in SpareElectronMonitoringGroups) :
        elCBLooseMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elCBLooseMonAlg')
    if ("CBLooseTrig" in SpareElectronMonitoringGroups) :
        elCBLooseTrigMonAlg = helper.addAlgorithm(MonitorElectronAlgorithm,'elCBLooseTrigMonAlg')

    ### STEP 3 ###
    # Edit properties of algorithms

    elLHTightMonAlg.ParticleContainerName = "Electrons"
    elLHTightMonAlg.RecoName = "LHTight"
    elLHTightMonAlg.ParticlePrefix = "electron"
    elLHTightMonAlg.WithTrigger = ""

    elLHTightTrigMonAlg.ParticleContainerName = "Electrons"
    elLHTightTrigMonAlg.RecoName = "LHTight"
    elLHTightTrigMonAlg.ParticlePrefix = "electron"
    elLHTightTrigMonAlg.WithTrigger = "WithTrigger"
    # to enable a trigger filter, for example:
    elLHTightTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    elLHLooseMonAlg.ParticleContainerName = "Electrons"
    elLHLooseMonAlg.RecoName = "LHLoose"
    elLHLooseMonAlg.ParticlePrefix = "electron"

    elLHLooseTrigMonAlg.ParticleContainerName = "Electrons"
    elLHLooseTrigMonAlg.RecoName = "LHLoose"
    elLHLooseTrigMonAlg.ParticlePrefix = "electron"
    elLHLooseTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    if ("CBTight" in SpareElectronMonitoringGroups) :
        elCBTightMonAlg.ParticleContainerName = "Electrons"
        elCBTightMonAlg.RecoName = "Tight"
        elCBTightMonAlg.ParticlePrefix = "electron"

    if ("CBTightTrig" in SpareElectronMonitoringGroups) :
        elCBTightTrigMonAlg.ParticleContainerName = "Electrons"
        elCBTightTrigMonAlg.RecoName = "Tight"
        elCBTightTrigMonAlg.ParticlePrefix = "electron"
        elCBTightTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    if ("CBLoose" in SpareElectronMonitoringGroups) :
        elCBLooseMonAlg.ParticleContainerName = "Electrons"
        elCBLooseMonAlg.RecoName = "Loose"
        elCBLooseMonAlg.ParticlePrefix = "electron"

    if ("CBLooseTrig" in SpareElectronMonitoringGroups) :
        elCBLooseTrigMonAlg.ParticleContainerName = "Electrons"
        elCBLooseTrigMonAlg.RecoName = "Loose"
        elCBLooseTrigMonAlg.ParticlePrefix = "electron"
        elCBLooseTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    ### STEP 4 ###

    GroupElectronLHTight = helper.addGroup(elLHTightMonAlg, 'MonitorElectron', 'egamma/LHTightElectrons/')
    GroupElectronLHTightTriggered = helper.addGroup(elLHTightTrigMonAlg, 'MonitorElectron', 'egamma/LHTightElectronsWithTrigger/')
    GroupElectronLHLoose = helper.addGroup(elLHLooseMonAlg, 'MonitorElectron', 'egamma/LHLooseElectrons/')
    GroupElectronLHLooseTriggered = helper.addGroup(elLHLooseTrigMonAlg, 'MonitorElectron', 'egamma/LHLooseElectronsWithTrigger/')
    if ("CBTight" in SpareElectronMonitoringGroups) :
        GroupElectronCBTight = helper.addGroup(elCBTightMonAlg, 'MonitorElectron', 'egamma/CBTightElectrons/')
    if ("CBTightTrig" in SpareElectronMonitoringGroups) :
        GroupElectronCBTightTriggered = helper.addGroup(elCBTightTrigMonAlg, 'MonitorElectron', 'egamma/CBTightTriggeredElectrons/')
    if ("CBLoose" in SpareElectronMonitoringGroups) :
        GroupElectronCBLoose = helper.addGroup(elCBLooseMonAlg, 'MonitorElectron', 'egamma/CBLooseElectrons/')
    if ("CBLooseTrig" in SpareElectronMonitoringGroups) :
        GroupElectronCBLooseTriggered = helper.addGroup(elCBLooseTrigMonAlg, 'MonitorElectron', 'egamma/CBLooseTriggeredElectrons/')

    ### STEP 5 ###
    # Configure histograms

    # LHTight Electrons
    BookHistograms(GroupElectronLHTight,"LHTight","Electron")
    BookHistograms(GroupElectronLHTightTriggered,"LHTightTrig","Electron","WithTrigger")
    # LHLoose Electrons
    BookHistograms(GroupElectronLHLoose,"LHLoose","Electron")
    BookHistograms(GroupElectronLHLooseTriggered,"LHLooseTrig","Electron","WithTrigger")
    # Cut Based Tight Electrons
    if ("CBTight" in SpareElectronMonitoringGroups) :
        BookHistograms(GroupElectronCBTight,"CBTight","Electron")
    if ("CBTightTrig" in SpareElectronMonitoringGroups) :
        BookHistograms(GroupElectronCBTightTriggered,"CBTightTrig","Electron","WithTrigger")
    # Cut Based Loose Electrons
    if ("CBLoose" in SpareElectronMonitoringGroups) :
        BookHistograms(GroupElectronCBLoose,"CBLoose","Electron")
    if ("CBLooseTrig" in SpareElectronMonitoringGroups) :
        BookHistograms(GroupElectronCBLooseTriggered,"CBLooseTrig","Electron","WithTrigger")

    ### STEP 6 ###
    return helper.result()


def MonitorPhotonConfig(inputFlags):
    '''Function to configures some algorithms in the monitoring system.'''

    ### STEP 1 ###

    from AthenaMonitoring import AthMonitorCfgHelper
    helper = AthMonitorCfgHelper(inputFlags,'PhotonAthMonitorCfg')

    ### STEP 2 ###

    MonitorPhotonAlgorithm=CompFactory.MonitorPhotonAlgorithm
    phCBTightMonAlg = helper.addAlgorithm(MonitorPhotonAlgorithm,'phCBTightMonAlg')
    phCBTightTrigMonAlg = helper.addAlgorithm(MonitorPhotonAlgorithm,'phCBTightTrigMonAlg')

    phCBLooseMonAlg = helper.addAlgorithm(MonitorPhotonAlgorithm,'phCBLooseMonAlg')
    phCBLooseTrigMonAlg = helper.addAlgorithm(MonitorPhotonAlgorithm,'phCBLooseTrigMonAlg')

    ### STEP 3 ###
    # Edit properties of algorithms

    phCBTightMonAlg.ParticleContainerName = "Photons"
    phCBTightMonAlg.RecoName = "Tight"
    phCBTightMonAlg.ParticlePrefix = "photon"

    phCBTightTrigMonAlg.ParticleContainerName = "Photons"
    phCBTightTrigMonAlg.RecoName = "Tight"
    phCBTightTrigMonAlg.ParticlePrefix = "photon"
    phCBTightTrigMonAlg.WithTrigger = "WithTrigger"
    # to enable a trigger filter, for example:
    phCBTightTrigMonAlg.TriggerChain = egammaConf.primary_double_pho[0]

    phCBLooseMonAlg.ParticleContainerName = "Photons"
    phCBLooseMonAlg.RecoName = "Loose"
    phCBLooseMonAlg.ParticlePrefix = "photon"

    phCBLooseTrigMonAlg.ParticleContainerName = "Photons"
    phCBLooseTrigMonAlg.RecoName = "Loose"
    phCBLooseTrigMonAlg.ParticlePrefix = "photon"
    phCBLooseTrigMonAlg.WithTrigger = "WithTrigger"
    # to enable a trigger filter, for example:
    phCBLooseTrigMonAlg.TriggerChain = egammaConf.primary_double_pho[0]


    ### STEP 4 ###

    GroupPhotonCBTight = helper.addGroup(phCBTightMonAlg, 'MonitorPhoton', 'egamma/CBTightPhotons/')
    GroupPhotonCBTightTriggered = helper.addGroup(phCBTightTrigMonAlg, 'MonitorPhoton', 'egamma/CBTightTriggeredPhotons/')
    GroupPhotonCBLoose = helper.addGroup(phCBLooseMonAlg, 'MonitorPhoton', 'egamma/CBLoosePhotons/')
    GroupPhotonCBLooseTriggered = helper.addGroup(phCBLooseTrigMonAlg, 'MonitorPhoton', 'egamma/CBLooseTriggeredPhotons/')

    ### STEP 5 ###
    # Configure histograms

    # Cut Based Tight Photons
    BookHistograms(GroupPhotonCBTight,"CBTight","Photon")
    BookHistograms(GroupPhotonCBTightTriggered,"CBTightTrig","Photon")
    # Cut Based Loose Photons
    BookHistograms(GroupPhotonCBLoose,"CBLoose","Photon")
    BookHistograms(GroupPhotonCBLooseTriggered,"CBLooseTrig","Photon")

    ### STEP 6 ###
    return helper.result()


def MonitorTnPConfig(inputFlags):
    '''Function to configures some algorithms in the monitoring system for Tag and Probe candidates (Z, J/Psi).'''

    ### STEP 1 ###

    from AthenaMonitoring import AthMonitorCfgHelper
    helper = AthMonitorCfgHelper(inputFlags,'TnPPAthMonitorCfg')


    ### STEP 2 ###

    MonitorTnPAlgorithm = CompFactory.MonitorTnPAlgorithm
    ZeeMonAlg = helper.addAlgorithm(MonitorTnPAlgorithm,'TnPZeeMonAlg')
    JPsiMonAlg = helper.addAlgorithm(MonitorTnPAlgorithm,'TnPJpsiMonAlg')


    ### STEP 3 ###
    # Edit properties of algorithms

    ZeeMonAlg.ParticleContainerName = "Electrons"
    ZeeMonAlg.RecoName = "LHLoose"
    ZeeMonAlg.ParticlePrefix = "electron"
    ZeeMonAlg.MassPeak = 91188.
    ZeeMonAlg.ElectronEtCut = 15000.
    ZeeMonAlg.MassLowerCut = 70000.
    ZeeMonAlg.MassUpperCut = 110000.
    ZeeMonAlg.TnPType = "Z"

    # get trigger chain from egammaConf.monitoring_Zee
    chain=""
    for el in egammaConf.monitoring_Zee:
        if chain != "":
            chain = chain + ", " + el
        else :
            chain = el
    #ZeeMonAlg.TriggerChain = chain

    JPsiMonAlg.ParticleContainerName = "Electrons"
    JPsiMonAlg.RecoName = "LHLoose"
    JPsiMonAlg.ParticlePrefix = "electron"
    JPsiMonAlg.MassPeak = 3097.
    JPsiMonAlg.ElectronEtCut = 3000.
    JPsiMonAlg.MassLowerCut = 2000.
    JPsiMonAlg.MassUpperCut = 5000.
    JPsiMonAlg.TnPType = "JPsi"

    # get trigger chain from egammaConf.monitoring_Jpsiee
    chain =""
    for el in egammaConf.monitoring_Jpsiee:
        if chain != "":
            chain = chain + ", " + el
        else :
            chain = el
    JPsiMonAlg.TriggerChain = chain

    ### STEP 4 ###

    GroupZee = helper.addGroup(ZeeMonAlg, 'MonitorTnP', 'egamma/TnPZ/')
    GroupJPsi = helper.addGroup(JPsiMonAlg, 'MonitorTnP', 'egamma/TnPJPsi/')

    ### STEP 5 ###
    # Configure histograms

    # Zee
    BookTnPHistograms(GroupZee, ZeeMonAlg.RecoName, ZeeMonAlg.TnPType, ZeeMonAlg.MassLowerCut, ZeeMonAlg.MassUpperCut)
    # JPsi
    BookTnPHistograms(GroupJPsi, JPsiMonAlg.RecoName, JPsiMonAlg.TnPType, JPsiMonAlg.MassLowerCut, JPsiMonAlg.MassUpperCut)

    ### STEP 6 ###
    return helper.result()


def MonitorForwardElectronConfig(inputFlags):
    '''Function to configures some algorithms in the monitoring system.'''

    ### STEP 1 ###

    from AthenaMonitoring import AthMonitorCfgHelper
    helper = AthMonitorCfgHelper(inputFlags,'FwdElectronAthMonitorCfg')

    ### STEP 2 ###

    SpareFwdElectronMonitoringGroups = [ "CBFwdTight" ]

    MonitorFwdElectronAlgorithm=CompFactory.MonitorFwdElectronAlgorithm
    # fwdelLHTightMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelLHTightMonAlg')
    # fwdelLHTightTrigMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelLHTightTrigMonAlg')
    if ("CBFwdTight" in SpareFwdElectronMonitoringGroups) :
        fwdelCBTightMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelCBTightMonAlg')
    if ("CBfwdTightTrig" in SpareFwdElectronMonitoringGroups) :
        fwdelCBTightTrigMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelCBTightTrigMonAlg')
    if ("CBLoose" in SpareFwdElectronMonitoringGroups) :
        fwdelCBLooseMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelCBLooseMonAlg')
    if ("CBLooseTrig" in SpareFwdElectronMonitoringGroups) :
        fwdelCBLooseTrigMonAlg = helper.addAlgorithm(MonitorFwdElectronAlgorithm,'fwdelCBLooseTrigMonAlg')

    ### STEP 3 ###
    # Edit properties of algorithms

    # fwdelLHTightMonAlg.ParticleContainerName = "ForwardElectrons"
    # fwdelLHTightMonAlg.RecoName = "LHTight"
    # fwdelLHTightMonAlg.ParticlePrefix = "FwdElectron"

    # fwdelLHTightTrigMonAlg.ParticleContainerName = "ForwardElectrons"
    # fwdelLHTightTrigMonAlg.RecoName = "LHTight"
    # fwdelLHTightTrigMonAlg.ParticlePrefix = "FwdElectron"
    # fwdelLHTightTrigMonAlg.WithTrigger = "WithTrigger"
    # to enable a trigger filter, for example:
    # fwdelLHTightTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    # fwdelLHLooseMonAlg.ParticleContainerName = "ForwardElectrons"
    # fwdelLHLooseMonAlg.RecoName = "LHLoose"
    # fwdelLHLooseMonAlg.ParticlePrefix = "FwdElectron"

    # fwdelLHLooseTrigMonAlg.ParticleContainerName = "ForwardElectrons"
    # fwdelLHLooseTrigMonAlg.RecoName = "LHLoose"
    # fwdelLHLooseTrigMonAlg.ParticlePrefix = "FwdElectron"
    # to enable a trigger filter, for example:
    # fwdelLHLooseTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    if ("CBFwdTight" in SpareFwdElectronMonitoringGroups) :
        fwdelCBTightMonAlg.ParticleContainerName = "ForwardElectrons"
        fwdelCBTightMonAlg.RecoName = "Tight"
        fwdelCBTightMonAlg.ParticlePrefix = "FwdElectron"

    if ("CBFwdTightTrig" in SpareFwdElectronMonitoringGroups) :
        fwdelCBTightTrigMonAlg.ParticleContainerName = "ForwardElectrons"
        fwdelCBTightTrigMonAlg.RecoName = "Tight"
        fwdelCBTightTrigMonAlg.ParticlePrefix = "FwdElectron"
        # to enable a trigger filter, for example:
        fwdelCBTightTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    if ("CBFwdLoose" in SpareFwdElectronMonitoringGroups) :
        fwdelCBLooseMonAlg.ParticleContainerName = "ForwardElectrons"
        fwdelCBLooseMonAlg.RecoName = "Loose"
        fwdelCBLooseMonAlg.ParticlePrefix = "FwdElectron"

    if ("CBFwdLooseTrig" in SpareFwdElectronMonitoringGroups) :
        fwdelCBLooseTrigMonAlg.ParticleContainerName = "ForwardElectrons"
        fwdelCBLooseTrigMonAlg.RecoName = "Loose"
        fwdelCBLooseTrigMonAlg.ParticlePrefix = "FwdElectron"
        # to enable a trigger filter, for example:
        fwdelCBLooseTrigMonAlg.TriggerChain = egammaConf.primary_single_ele[0]

    ### STEP 4 ###

    #GroupFwdElectronLHTight = helper.addGroup(fwdelLHTightMonAlg, 'MonitorFwdElectron', 'egamma/LHTightFwdElectrons/')
    #GroupFwdElectronLHTightTriggered = helper.addGroup(fwdelLHTightTrigMonAlg, 'MonitorFwdElectron', 'egamma/LHTightFwdElectronsWithTrigger/')
    #GroupFwdElectronLHLoose = helper.addGroup(fwdelLHLooseMonAlg, 'MonitorFwdElectron', 'egamma/LHLooseFwdElectrons/')
    #GroupFwdElectronLHLooseTriggered = helper.addGroup(fwdelLHLooseTrigMonAlg, 'MonitorFwdElectron', 'egamma/LHLooseFwdElectronsWithTrigger/')
    if ("CBFwdTight" in SpareFwdElectronMonitoringGroups) :
        GroupFwdElectronCBTight = helper.addGroup(fwdelCBTightMonAlg, 'MonitorFwdElectron', 'egamma/CBTightFwdElectrons/')
    if ("CBFwdTightTrig" in SpareFwdElectronMonitoringGroups) :
        GroupFwdElectronCBTightTriggered = helper.addGroup(fwdelCBTightTrigMonAlg, 'MonitorFwdElectron', 'egamma/CBTightTriggeredFwdElectrons/')
    if ("CBFwdLoose" in SpareFwdElectronMonitoringGroups) :
        GroupFwdElectronCBLoose = helper.addGroup(fwdelCBLooseMonAlg, 'MonitorFwdElectron', 'egamma/CBLooseFwdElectrons/')
    if ("CBFwdLooseTrig" in SpareFwdElectronMonitoringGroups) :
        GroupFwdElectronCBLooseTriggered = helper.addGroup(fwdelCBLooseTrigMonAlg, 'MonitorFwdElectron', 'egamma/CBLooseTriggeredFwdElectrons/')

    ### STEP 5 ###
    # Configure histograms

    # LHTight FwdElectrons
    #BookFwdHistograms(GroupFwdElectronLHTight,"LHTight","FwdElectron")
    #BookFwdHistograms(GroupFwdElectronLHTightTriggered,"LHTightTrig","FwdElectron","WithTrigger")
    # LHLoose FwdElectrons
    #BookFwdHistograms(GroupFwdElectronLHLoose,"LHLoose","FwdElectron")
    #BookFwdHistograms(GroupFwdElectronLHLooseTriggered,"LHLooseTrig","FwdElectron","WithTrigger")
    # Cut Based Tight FwdElectrons
    if ("CBFwdTight" in SpareFwdElectronMonitoringGroups) :
        BookFwdElectronHistograms(GroupFwdElectronCBTight,"CBTight","FwdElectron")
    if ("CBFwdTightTrig" in SpareFwdElectronMonitoringGroups) :
        BookFwdElectronHistograms(GroupFwdElectronCBTightTriggered,"CBTightTrig","FwdElectron","WithTrigger")
    # Cut Based Loose FwdElectrons
    if ("CBFwdLoose" in SpareFwdElectronMonitoringGroups) :
        BookFwdElectronHistograms(GroupFwdElectronCBLoose,"CBLoose","FwdElectron")
    if ("CBFwdLooseTrig" in SpareFwdElectronMonitoringGroups) :
        BookFwdElectronHistograms(GroupFwdElectronCBLooseTriggered,"CBLooseTrig","FwdElectron","WithTrigger")

    ### STEP 6 ###
    return helper.result()


if __name__=='__main__':
    # Setup logs
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    log.setLevel(INFO)
    
    # Set the Athena configuration flags
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    nightly = '/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CommonInputs/'
    file = 'data16_13TeV.00311321.physics_Main.recon.AOD.r9264/AOD.11038520._000001.pool.root.1'

    ConfigFlags.Input.Files = [nightly+file]
    # ConfigFlags.Input.Files = ['AOD.11038520._000001.pool.root.1']
    ConfigFlags.Input.isMC = False
    # ConfigFlags.Output.HISTFileName = 'MonitorEgammaOutput.root'

# To produce WebDisplay, filename must follow a certain format
    ConfigFlags.Output.HISTFileName = 'data16_13TeV.00311321.physics_Main.merge.HIST.f1156_h347._0002.1.root'

    ConfigFlags.lock()

    # Initialize configuration object, add accumulator, merge, and run.
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg = MainServicesCfg(ConfigFlags)
    cfg.merge(PoolReadCfg(ConfigFlags))

    MonitorElectronAcc = MonitorElectronConfig(ConfigFlags)
    cfg.merge(MonitorElectronAcc)

    MonitorPhotonAcc = MonitorPhotonConfig(ConfigFlags)
    cfg.merge(MonitorPhotonAcc)

    MonitorTnPAcc = MonitorTnPConfig(ConfigFlags)
    cfg.merge(MonitorTnPAcc)

    MonitorFwdElectronAcc = MonitorForwardElectronConfig(ConfigFlags)
    cfg.merge(MonitorFwdElectronAcc)

    # If you want to turn on more detailed messages ...

    #MonitorElectronAcc.getEventAlgo('elLHLooseMonAlg').OutputLevel = 2 # 2 = DEBUG
    #MonitorPhotonAcc.getEventAlgo('phCBLooseMonAlg').OutputLevel = 2 # 2 = DEBUG
    #MonitorTnPAcc.getEventAlgo('TnPZeeMonAlg').OutputLevel = 2 # 2 = DEBUG
    #MonitorTnPAcc.getEventAlgo('TnPJpsiMonAlg').OutputLevel = 2 # 2 = DEBUG
    #MonitorFwdElectronAcc.getEventAlgo('fwdelCBTightMonAlg').OutputLevel = 2 # 2 = DEBUG

    MonitorFwdElectronAcc = MonitorForwardElectronConfig(ConfigFlags)
    cfg.merge(MonitorFwdElectronAcc)

    cfg.printConfig(withDetails=False) # set True for exhaustive info
    #cfg.printConfig(withDetails=True) # set True for exhaustive info

    #cfg.run(10)
    cfg.run(-1)  # the one to process all events  
    #cfg.run(600) # use cfg.run(20) to only run on first 20 events
