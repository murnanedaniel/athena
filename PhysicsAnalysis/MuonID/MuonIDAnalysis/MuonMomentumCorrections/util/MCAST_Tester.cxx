/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

// System include(s):
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>

// ROOT include(s):
#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

// EDM include(s):
#include "AsgMessaging/Check.h"
#include "AsgMessaging/MessageCheck.h"
#include "AsgMessaging/StatusCode.h"
#include "AsgTools/StandaloneToolHandle.h"
#include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicVariation.h"
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/tools/IOStats.h"
#include "xAODCore/tools/ReadStats.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"

ANA_MSG_HEADER(msgMMC)

#include <MuonMomentumCorrections/CalibContainer.h>
#include <MuonMomentumCorrections/CalibInitializer.h>


int main(int argc, char* argv[]) {

    const char* APP_NAME = argv[0];

    // setup for ANA_CHECK()
     using namespace msgMMC;
    ANA_CHECK_SET_TYPE(int);

    bool useCorrectedCopy = true;

    // The application's name:

    // Check if we received a file name:
    if (argc < 2) {
        Error(APP_NAME, "==============================================");
        Error(APP_NAME, "No file name received!");
        Error(APP_NAME, "Usage: $> %s [xAOD file name]", APP_NAME);
        Error(APP_NAME, " $> %s [xAOD file name]", APP_NAME);
        Error(APP_NAME, " $> %s [xAOD file name] -n X  |  X = number of events you want to run on", APP_NAME);
        Error(APP_NAME, " $> %s [xAOD file name] -event X  |  X = specific number of the event to run on - for debugging", APP_NAME);
        Error(APP_NAME, "==============================================");
        return 1;
    }

    ////////////////////////////////////////////////////
    //:::  parse the options
    ////////////////////////////////////////////////////
    std::string options;
    for (int i = 0; i < argc; i++) { options += (argv[i]); }

    int Ievent = -1;
    
    if (options.find("-event") != std::string::npos) {
        for (int ipos = 0; ipos < argc; ipos++) {
            if (std::string(argv[ipos]).compare("-event") == 0) {
                Ievent = atoi(argv[ipos + 1]);
                Info(APP_NAME, "Argument (-event) : Running only on event # %i", Ievent);
                break;
            }
        }
    }

    int nEvents = -1;
    if (options.find("-n") != std::string::npos) {
        for (int ipos = 0; ipos < argc; ipos++) {
            if (std::string(argv[ipos]).compare("-n") == 0) {
                nEvents = atoi(argv[ipos + 1]);
                Info(APP_NAME, "Argument (-n) : Running on NEvents = %i", nEvents);
                break;
            }
        }
    }

    TString limitSys = "";
    if (options.find("-s") != std::string::npos) {
        for (int ipos = 0; ipos < argc; ipos++) {
            if (std::string(argv[ipos]).compare("-s") == 0) {
                limitSys = TString(argv[ipos + 1]);
                break;
            }
        }
    }

    bool doSys = false;
    if (options.find("-doSys") != std::string::npos) {
       doSys = true;
    }
    ////////////////////////////////////////////////////
    //:::  initialize the application and get the event
    ////////////////////////////////////////////////////
    ANA_CHECK(xAOD::Init(APP_NAME));

    //::: Open the input file:
    std::string fileName = argv[1];
    Info(APP_NAME, "Opening file: %s", fileName.c_str());
    std::unique_ptr<TFile> ifile(TFile::Open(fileName.c_str(), "READ"));
    if (!ifile) Error(APP_NAME, "Cannot find file %s", fileName.c_str());

    std::unique_ptr<TChain> chain = std::make_unique<TChain>("CollectionTree", "CollectionTree");
    chain->Add(fileName.c_str());

    //::: Create a TEvent object:
    xAOD::TEvent event((TTree*)chain.get(), xAOD::TEvent::kAthenaAccess);
    Info(APP_NAME, "Number of events in the file: %i", static_cast<int>(event.getEntries()));

    //::: Create a transient object store. Needed for the tools.
    xAOD::TStore store;

    //::: Decide how many events to run over:
    Long64_t entries = event.getEntries();
    if (fileName.find("data") != std::string::npos) entries = 2000;
    // if (fileName.find("mc20") != std::string::npos) entries = 2000;

    if(doSys) entries = 200;

    ////////////////////////////////////////////////////
    //::: MuonCalibrationAndSmearingTool
    // setup the tool handle as per the
    // recommendation by ASG - https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AthAnalysisBase#How_to_use_AnaToolHandle
    ////////////////////////////////////////////////////
    //::: create the tool handle
    asg::StandaloneToolHandle<CP::IMuonCalibrationAndSmearingTool> corrTool;  //!
    corrTool.setTypeAndName("CP::MuonCalibrationAndSmearingTool/MCT");
    //::: set the properties
    if (fileName.find("data18") != std::string::npos) corrTool.setProperty("Year", "Data18").ignore();
    if (fileName.find("data17") != std::string::npos) corrTool.setProperty("Year", "Data17").ignore();
    if (fileName.find("data16") != std::string::npos) corrTool.setProperty("Year", "Data16").ignore();
    if (fileName.find("data15") != std::string::npos) corrTool.setProperty("Year", "Data16").ignore();

    if (fileName.find("r13145") != std::string::npos) corrTool.setProperty("Year", "Data18").ignore();
    if (fileName.find("r13144") != std::string::npos) corrTool.setProperty("Year", "Data17").ignore();
    if (fileName.find("r13167") != std::string::npos) corrTool.setProperty("Year", "Data16").ignore();
    corrTool.setProperty("systematicCorrelationScheme", "Corr_Scale").ignore();
    corrTool.setProperty("SagittaCorr", true).ignore();
    corrTool.setProperty("doSagittaMCDistortion", false).ignore();
    corrTool.setProperty("doDirectCBCalib", false).ignore();
    // corrTool.setProperty("doExtraSmearing", true).ignore();
    // corrTool.setProperty("do2StationsHighPt", true).ignore();

    asg::StandaloneToolHandle<CP::IMuonCalibrationAndSmearingTool> newcorrTool;  //!
    newcorrTool.setProperty("calibMode", 1).ignore();
    // newcorrTool.setProperty("doExtraSmearing", true).ignore();
    // newcorrTool.setProperty("do2StationsHighPt", false).ignore();

    newcorrTool.setTypeAndName("CP::MuonCalibTool/NewMCT");



    bool isDebug = false;
    if (nEvents >= 0 || Ievent >= 0) isDebug = true;

    if (isDebug) corrTool.setProperty("OutputLevel", MSG::VERBOSE).ignore();
    if (isDebug) newcorrTool.setProperty("OutputLevel", MSG::VERBOSE).ignore();
    //::: retrieve the tool
    StatusCode sc = newcorrTool.retrieve();
    if (sc.isFailure()) {
        Error(APP_NAME, "Cannot retrieve MuonCorrectionTool");
        return 1;
    }

    //::: retrieve the tool
    sc = corrTool.retrieve();
    if (sc.isFailure()) {
        Error(APP_NAME, "Cannot retrieve MuonCorrectionTool");
        return 1;
    }


    ////////////////////////////////////////////////////
    //::: MuonSelectionTool
    // setup the tool handle as per the
    // recommendation by ASG - https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AthAnalysisBase#How_to_use_AnaToolHandle
    ////////////////////////////////////////////////////
    //::: create the tool handle
    asg::StandaloneToolHandle<CP::IMuonSelectionTool> selTool;  //!
    selTool.setTypeAndName("CP::MuonSelectionTool/MuonSelectionTool");

    //::: set the properties
    selTool.setProperty("MaxEta", 2.5).ignore();
    selTool.setProperty("MuQuality", (int)xAOD::Muon::Loose).ignore();  // corresponds to 0=Tight, 1=Medium, 2=Loose, 3=VeryLoose, 4=HighPt, 5=LowPtEfficiency

    //::: retrieve the tool
    if (selTool.retrieve().isFailure() || sc.isFailure()) {
        Error(APP_NAME, "Cannot retrieve MuonSelectionTool");
        return 1;
    }

    ////////////////////////////////////////////////////
    //::: Systematics initialization
    ////////////////////////////////////////////////////
    std::vector<CP::SystematicSet> sysList;
    if(limitSys.Length() == 0)sysList.push_back(CP::SystematicSet());

    if(doSys)
    {
        const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
        const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics();
        for (CP::SystematicSet::const_iterator sysItr = recommendedSystematics.begin(); sysItr != recommendedSystematics.end(); ++sysItr) {
            if(!TString(sysItr->name()).Contains(limitSys) && limitSys.Length() > 0) continue;
            sysList.push_back(CP::SystematicSet());
            sysList.back().insert(*sysItr);
        }
    }

    std::vector<CP::SystematicSet>::const_iterator sysListItr;

    msgMMC::ANA_MSG_INFO("main() - Systematics are ");
    for (sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr) msgMMC::ANA_MSG_INFO(sysListItr->name());

    // branches to be stored
    Float_t InitPtCB(0.), InitPtID(0.), InitPtMS(0.);
    Float_t CorrPtCB(0.), CorrPtID(0.), CorrPtMS(0.);
    Float_t Eta(0.), Phi(0.), Charge(0.);
    Float_t ExpResoCB(0.), ExpResoID(0.), ExpResoMS(0.);
    long long unsigned int eventNum(0);

    // output file
    TFile* outputFile = TFile::Open("output.root", "recreate");

    // it contains a single TTree per systematic variation
    std::unordered_map<CP::SystematicSet, TTree*> sysTreeMap;
    for (sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr) {
        // create new tree for the systematic in question
        std::string treeName = "test_tree_" + (sysListItr->name().size() == 0 ? "NOMINAL" : sysListItr->name());
        TTree* sysTree = new TTree(treeName.c_str(), "test tree for MCAST");

        // add branches
        sysTree->Branch("InitPtCB", &InitPtCB, "InitPtCB/F");
        sysTree->Branch("InitPtID", &InitPtID, "InitPtID/F");
        sysTree->Branch("InitPtMS", &InitPtMS, "InitPtMS/F");
        sysTree->Branch("CorrPtCB", &CorrPtCB, "CorrPtCB/F");
        sysTree->Branch("CorrPtID", &CorrPtID, "CorrPtID/F");
        sysTree->Branch("CorrPtMS", &CorrPtMS, "CorrPtMS/F");
        sysTree->Branch("Eta", &Eta, "Eta/F");
        sysTree->Branch("Phi", &Phi, "Phi/F");
        sysTree->Branch("Charge", &Charge, "Charge/F");
        sysTree->Branch("ExpResoCB", &ExpResoCB, "ExpResoCB/F");
        sysTree->Branch("ExpResoID", &ExpResoID, "ExpResoID/F");
        sysTree->Branch("ExpResoMS", &ExpResoMS, "ExpResoMS/F");
        sysTree->Branch("eventNum", &eventNum);

        sysTreeMap[*sysListItr] = sysTree;
    }

    ////////////////////////////////////////////////////
    //::: Loop over the events
    ////////////////////////////////////////////////////
    for (Long64_t entry = 0; entry < entries; ++entry) {
        if (nEvents != -1 && entry > nEvents) break;
        // Tell the object which entry to look at:
        event.getEntry(entry);

        // Print some event information
        const xAOD::EventInfo* evtInfo = 0;
        ANA_CHECK(event.retrieve(evtInfo, "EventInfo"));
        if (Ievent != -1 && static_cast<int>(evtInfo->eventNumber()) != Ievent) { continue; }

        eventNum = evtInfo->eventNumber();

        //::: Get the Muons from the event:
        const xAOD::MuonContainer* muons = 0;
        ANA_CHECK(event.retrieve(muons, "Muons"));

        //::: Loop over systematics
        for (sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr) {
            // create a shallow copy of the muons container
            std::pair<xAOD::MuonContainer*, xAOD::ShallowAuxContainer*> muons_shallowCopy = xAOD::shallowCopyContainer(*muons);

            xAOD::MuonContainer* muonsCorr = muons_shallowCopy.first;

            if (isDebug) {
                Info(APP_NAME, "-----------------------------------------------------------");
                Info(APP_NAME, "Looking at %s systematic", (sysListItr->name()).c_str());
            }
            //::: Check if systematic is applicable to the correction tool
            if (corrTool->applySystematicVariation(*sysListItr) != StatusCode::SUCCESS) {
                Error(APP_NAME, "Cannot configure muon calibration tool for systematic");
            }
            if (newcorrTool->applySystematicVariation(*sysListItr) != StatusCode::SUCCESS) {
                Error(APP_NAME, "Cannot configure muon calibration tool for systematic");
            }

            //::: Loop over muon container
            for (auto muon : *muonsCorr) {
                //::: Select "good" muons:
                if (!selTool->accept(*muon)) {
                    if (isDebug) Info(APP_NAME, "This muon doesn't pass the ID hits quality cuts");
                    continue;
                }

                //::: Should be using correctedCopy here, testing behaviour of applyCorrection though
                InitPtCB = muon->pt();
                InitPtID = -999;
                if (muon->inDetTrackParticleLink().isValid()) {
                    const ElementLink<xAOD::TrackParticleContainer>& id_track = muon->inDetTrackParticleLink();
                    InitPtID = (!id_track) ? 0 : (*id_track)->pt();
                }
                InitPtMS = -999;
                if (muon->extrapolatedMuonSpectrometerTrackParticleLink().isValid()) {
                    const ElementLink<xAOD::TrackParticleContainer>& ms_track = muon->extrapolatedMuonSpectrometerTrackParticleLink();
                    InitPtMS = (!ms_track) ? 0 : (*ms_track)->pt();
                }

                Eta = muon->eta();
                Phi = muon->phi();
                Charge = muon->charge();

                //::: Print some info about the selected muon:
                if (isDebug) Info(APP_NAME, "Selected muon: eta = %g, phi = %g, pt = %g", muon->eta(), muon->phi(), muon->pt() / 1e3);

                float ptCB = 0;
                if (muon->primaryTrackParticleLink().isValid()) {
                    const ElementLink<xAOD::TrackParticleContainer>& cb_track = muon->primaryTrackParticleLink();
                    ptCB = (!cb_track) ? 0 : (*cb_track)->pt();
                } else {
                    if (isDebug)
                        Info(APP_NAME, "Missing primary track particle link for --> CB %g, author: %d, type: %d", ptCB, muon->author(),
                             muon->muonType());
                }
                float ptID = 0;
                if (muon->inDetTrackParticleLink().isValid()) {
                    const ElementLink<xAOD::TrackParticleContainer>& id_track = muon->inDetTrackParticleLink();
                    ptID = (!id_track) ? 0 : (*id_track)->pt();
                }
                float ptME = 0;
                if (muon->extrapolatedMuonSpectrometerTrackParticleLink().isValid()) {
                    const ElementLink<xAOD::TrackParticleContainer>& ms_track = muon->extrapolatedMuonSpectrometerTrackParticleLink();
                    ptME = (!ms_track) ? 0 : (*ms_track)->pt();
                }

                if (isDebug)
                    Info(APP_NAME, "--> CB %g, ID %g, ME %g, author: %d, type: %d", ptCB / 1e3, ptID / 1e3, ptME / 1e3, muon->author(),
                         muon->muonType());

                // either use the correctedCopy call or correct the muon object itself
                if (useCorrectedCopy) {
                    //::: Create a calibrated muon:
                    xAOD::Muon* mu = 0;
                    if (!corrTool->correctedCopy(*muon, mu)) {
                        Error(APP_NAME, "Cannot really apply calibration nor smearing");
                        continue;
                    }
                    CorrPtCB = mu->pt();
                    CorrPtID = mu->auxdata<float>("InnerDetectorPt");
                    CorrPtMS = mu->auxdata<float>("MuonSpectrometerPt");


                    xAOD::Muon* muNew = 0;
                    if (!newcorrTool->correctedCopy(*muon, muNew)) {
                        Error(APP_NAME, "Cannot really apply calibration nor smearing");
                        continue;
                    }
                    double NewCorrPtCB = muNew->pt();
                    double NewCorrPtID = muNew->auxdata<float>("InnerDetectorPt");
                    double NewCorrPtMS = muNew->auxdata<float>("MuonSpectrometerPt");

                    if (isDebug)
                    {
                        Info(APP_NAME, "Old Calibrated muon: eta = %g, phi = %g, pt(CB) = %g, pt(ID) = %g, pt(MS) = %g", mu->eta(), mu->phi(), mu->pt() / 1e3, mu->auxdata<float>("InnerDetectorPt") / 1e3, mu->auxdata<float>("MuonSpectrometerPt") / 1e3);
                        Info(APP_NAME, "new Calibrated muon: eta = %g, phi = %g, pt(CB) = %g, pt(ID) = %g, pt(MS) = %g", muNew->eta(), muNew->phi(), muNew->pt() / 1e3, muNew->auxdata<float>("InnerDetectorPt") / 1e3, muNew->auxdata<float>("MuonSpectrometerPt") / 1e3);
                    }
                    if(std::abs(NewCorrPtCB-CorrPtCB) > 0.5) Warning(APP_NAME, "CB pt not matching: old = %g, new = %g, event=%lli, eta = %g, phi = %g, for %s", CorrPtCB, NewCorrPtCB, eventNum, mu->eta(), mu->phi(),  sysListItr->name().c_str());
                    if(std::abs(NewCorrPtID-CorrPtID) > 0.5) Warning(APP_NAME, "ID pt not matching: old = %g, new = %g, event=%lli, eta = %g, phi = %g, for %s", CorrPtID, NewCorrPtID, eventNum, mu->eta(), mu->phi(),  sysListItr->name().c_str());
                    if(std::abs(NewCorrPtMS-CorrPtMS) > 0.5) Warning(APP_NAME, "MS pt not matching: old = %g, new = %g, event=%lli, eta = %g, phi = %g, for %s", CorrPtMS, NewCorrPtMS, eventNum, mu->eta(), mu->phi(),  sysListItr->name().c_str());

                    
                    sysTreeMap[*sysListItr]->Fill();
                    //::: Delete the calibrated muon:
                    delete mu;
                } else {
                    if (!corrTool->applyCorrection(*muon)) {
                        Error(APP_NAME, "Cannot really apply calibration nor smearing");
                        continue;
                    }
                    CorrPtCB = muon->pt();
                    CorrPtID = muon->auxdata<float>("InnerDetectorPt");
                    CorrPtMS = muon->auxdata<float>("MuonSpectrometerPt");
                    ExpResoCB = corrTool->expectedResolution("CB", *muon, true);
                    ExpResoID = corrTool->expectedResolution("ID", *muon, true);
                    ExpResoMS = corrTool->expectedResolution("MS", *muon, true);
                    if (isDebug)
                        Info(APP_NAME, "Calibrated muon: eta = %g, phi = %g, pt(CB) = %g, pt(ID) = %g, pt(MS) = %g", muon->eta(),
                             muon->phi(), muon->pt() / 1e3, CorrPtID / 1e3, CorrPtMS / 1e3);
                    if (isDebug)
                        Info(APP_NAME, " expReso : ExpResoCB = %g , ExpResoID = %g , ExpResoMS = %g", ExpResoCB, ExpResoID, ExpResoMS);
                    sysTreeMap[*sysListItr]->Fill();
                }
            }
            if (isDebug) Info(APP_NAME, "-----------------------------------------------------------");

            delete muons_shallowCopy.first;
            delete muons_shallowCopy.second;
        }
        //::: Close with a message:
        if (entry % 5 == 0)
            Info(APP_NAME,
                 "===>>>  done processing event #%i, run #%i %i events processed so far  <<<===", static_cast<int>(evtInfo->eventNumber()),
                 static_cast<int>(evtInfo->runNumber()), static_cast<int>(entry + 1));
    }

    for (sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr) { sysTreeMap[*sysListItr]->Write(); }

    //::: Close output file
    outputFile->Close();

    xAOD::IOStats::instance().stats().printSmartSlimmingBranchList(); 

    //::: Return gracefully:
    return 0;
}
