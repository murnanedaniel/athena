/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

// MuonPhysValMonitoringTool.cxx
// Implementation file for class MuonPhysValMonitoringTool
// Author: Felix Socher <Felix.Socher@cern.ch>
///////////////////////////////////////////////////////////////////

#ifndef MUONPHYSVALMONITORING_MUONPHYSVALMONITORINGTOOL_H
#define MUONPHYSVALMONITORING_MUONPHYSVALMONITORINGTOOL_H



#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/ServiceHandle.h"
#include "IsolationSelection/IIsolationSelectionTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "MuonCombinedToolInterfaces/IMuonPrintingTool.h"


#include "StoreGate/ReadHandleKey.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

#include "TrkToolInterfaces/ITrackSelectorTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/Muon.h" //typedef
#include "xAODMuon/MuonSegment.h" //typedef
#include "xAODMuon/SlowMuon.h" //typedef
#include "xAODTrigMuon/L2CombinedMuon.h" //typedef

#include "xAODTrigMuon/L2StandAloneMuon.h" //typedef
#include "xAODTruth/TruthParticle.h" //typedef
#include "xAODTruth/TruthParticleContainer.h" //typedef

#include "MuonSegmentValidationPlots.h" //needed for unique_ptr  access to deleter
#include "MuonHistUtils/MuonSegmentPlots.h" //needed for unique_ptr  access to deleter
#include "MuonValidationPlots.h" //needed for unique_ptr  access to deleter
#include "MuonHistUtils/RecoMuonPlotOrganizer.h" //needed for unique_ptr  access to deleter
#include "MuonHistUtils/TruthMuonPlotOrganizer.h" //needed for unique_ptr  access to deleter
#include "TrkValHistUtils/PlotBase.h" //for PlotBase, also for typedef of HistData
#include "SlowMuonValidationPlots.h" //needed for unique_ptr  access to deleter
#include "MuonTrackValidationPlots.h" //needed for unique_ptr  access to deleter
#include "TriggerMuonValidationPlots.h" //needed for unique_ptr  access to deleter

#include <string>
#include <vector>
#include <memory>

class TString;
class TH1F;



namespace MuonPhysValMonitoring {

    class MuonPhysValMonitoringTool : public ManagedMonitorToolBase {
        ///////////////////////////////////////////////////////////////////
        // Public methods:
        ///////////////////////////////////////////////////////////////////
    public:
        // Copy constructor:

        /// Constructor with parameters:
        MuonPhysValMonitoringTool(const std::string& type, const std::string& name, const IInterface* parent);

        /// Destructor:
        virtual ~MuonPhysValMonitoringTool() = default;

        // Athena algtool's Hooks
        virtual StatusCode initialize() override;       
        virtual StatusCode bookHistograms() override;
        virtual StatusCode fillHistograms() override;
        virtual StatusCode procHistograms() override;

        ///////////////////////////////////////////////////////////////////
        // Const methods:
        ///////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////
        // Non-const methods:
        ///////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////
        // Private data:
        ///////////////////////////////////////////////////////////////////
    private:
        enum MUCATEGORY { ALL = 0, PROMPT, INFLIGHT, NONISO, REST };

        void handleMuon(const xAOD::Muon* mu, const xAOD::SlowMuon* smu = nullptr, float weight = 1.0);
        void handleSlowMuon(const xAOD::SlowMuon* smu, float weight = 1.0);
        void handleTruthMuon(const xAOD::TruthParticle* truthMu, float weight = 1.0);
        void handleMuonTrack(const xAOD::TrackParticle* tp, xAOD::Muon::TrackParticleType type, float weight = 1.0);
        void handleMuonSegment(const xAOD::MuonSegment* muSeg, float weight = 1.0);
        void handleTruthMuonSegment(const xAOD::MuonSegment* truthMuSeg, const xAOD::TruthParticleContainer* muonTruthContainer,
                                    float weight = 1.0);

        void handleMuonTrees(const xAOD::EventInfo* eventInfo, bool isData);

        void handleMuonL1Trigger(const xAOD::MuonRoI* TrigL1mu);
        void handleMuonL2Trigger(const xAOD::L2StandAloneMuon* L2SAMu);
        void handleMuonL2Trigger(const xAOD::L2CombinedMuon* L2CBMu);
        void handleMuonTrigger(const xAOD::Muon* mu);
        void L2SATriggerResolution();
        void L2CBTriggerResolution();
        void EFTriggerResolution();

        void printMuonDebug(const xAOD::Muon* mu);
        void printTruthMuonDebug(const xAOD::TruthParticle* truthMu, const xAOD::Muon* mu);

        StatusCode bookValidationPlots(PlotBase& valPlots);
        const xAOD::Muon* findRecoMuon(const xAOD::TruthParticle* truthMu);
        const xAOD::SlowMuon* findRecoSlowMuon(const xAOD::TruthParticle* truthMu);
        const xAOD::MuonSegment* findRecoMuonSegment(const xAOD::MuonSegment* truthMuSeg);
        std::unique_ptr<xAOD::Muon> getCorrectedMuon(const xAOD::Muon& mu);

        const xAOD::TrackParticleContainer* m_MSTracks{nullptr};
        std::map<std::string, int> m_counterBits;
        std::vector<std::string> m_muonItems;
        std::vector<std::string> m_L1Seed;
        int m_SelectedAuthor{0};

        TH1F* findHistogram(const std::vector<HistData>& hists, const std::string& hnameTag, const std::string& hdirTag,
                            const std::string& hNewName);
        void modifyHistogram(TH1* hist);

        Gaudi::Property<std::string> m_tracksName{this, "TrackContainerName", "InDetTrackParticles"};
        Gaudi::Property<std::string> m_fwdtracksName{this, "FwdTrackContainerName", ""};
        Gaudi::Property<std::string> m_muonsName{this, "MuonContainerName", "Muons"};
        Gaudi::Property<std::string> m_slowMuonsName{this, "SlowMuonContainerName", "SlowMuons"};
        Gaudi::Property<std::string> m_muonsTruthName{this, "MuonTruthParticleContainerName", "MuonTruthParticles"};
        Gaudi::Property<std::string> m_muonTracksName{this, "MuonTrackContainerName", "MuonSpectrometerTrackParticles"};
        Gaudi::Property<std::string> m_muonExtrapolatedTracksName{this, "MuonExtrapolatedTrackContainerName",
                                                                  "ExtrapolatedMuonTrackParticles"};
        Gaudi::Property<std::string> m_muonMSOnlyExtrapolatedTracksName{this, "MuonOnlyExtrapolatedTrackContainerName",
                                                                        "MSOnlyExtrapolatedMuonTrackParticles"};
        Gaudi::Property<std::string> m_muonSegmentsName{this, "MuonSegmentContainerName", "MuonSegments"};
        Gaudi::Property<std::string> m_muonSegmentsTruthName{this, "MuonTruthSegmentContainerName", "MuonTruthSegments"};
        Gaudi::Property<std::string> m_muonL1TrigName{this, "L1TrigMuonContainerName", "LVL1MuonRoIs"};
        Gaudi::Property<std::string> m_muonL2SAName{this, "L2SAMuonContainerName", "HLT_xAOD__L2StandAloneMuonContainer_MuonL2SAInfo"};
        Gaudi::Property<std::string> m_muonL2CBName{this, "L2CBMuonContainerName", "HLT_xAOD__L2CombinedMuonContainer_MuonL2CBInfo"};
        Gaudi::Property<std::string> m_muonEFCombTrigName{this, "EFCombTrigMuonContainerName", "HLT_xAOD__MuonContainer_MuonEFInfo"};

        Gaudi::Property<std::vector<int>> m_selectMuonWPs{
            this, "SelectMuonWorkingPoints", {xAOD::Muon::Loose, xAOD::Muon::Medium, xAOD::Muon::Tight}};
        Gaudi::Property<std::vector<unsigned int>> m_selectMuonAuthors{
            this,
            "SelectMuonAuthors",
            {xAOD::Muon::MuidCo, xAOD::Muon::MuTagIMO, xAOD::Muon::MuidSA, xAOD::Muon::MuGirl, xAOD::Muon::CaloTag, xAOD::Muon::CaloScore}};
        /// Flag to tell whether muons with the comissioning author will be selected or not
        Gaudi::Property<bool> m_selectComissioning{this, "SelectComissioningMuons", false};
        Gaudi::Property<std::vector<std::vector<std::string>>> m_selectHLTMuonItems{this, "SelectHLTMuonItems", {}};
        Gaudi::Property<std::vector<std::string>> m_L1MuonItems{this, "SelectL1MuonItems", {}};
        Gaudi::Property<std::vector<unsigned int>> m_selectMuonCategories{
            this,
            "SelectMuonCategories",
            {MUCATEGORY::ALL, MUCATEGORY::PROMPT, MUCATEGORY::INFLIGHT, MUCATEGORY::NONISO, MUCATEGORY::REST}};

        Gaudi::Property<bool> m_doBinnedResolutionPlots{this, "DoBinnedResolutionPlots", true};
        Gaudi::Property<bool> m_doTrigMuonValidation{this, "DoTrigMuonValidation", false};
        Gaudi::Property<bool> m_doTrigMuonL1Validation{this, "DoTrigMuonL1Validation", false};
        Gaudi::Property<bool> m_doTrigMuonL2Validation{this, "DoTrigMuonL2Validation", false};
        Gaudi::Property<bool> m_doTrigMuonEFValidation{this, "DoTrigMuonEFValidation", false};
        Gaudi::Property<bool> m_doMuonTree{this, "DoMuonTree", false};
        Gaudi::Property<bool> m_isData{this, "IsData", false};

        SG::ReadHandleKey<xAOD::EventInfo> m_eventInfo{this, "EventInfo", "EventInfo", "event info"};

        // Tools
        ToolHandle<CP::IMuonSelectionTool> m_muonSelectionTool{this, "MuonSelector", "CP::MuonSelectionTool/MuonSelectionTool"};
        ToolHandle<Rec::IMuonPrintingTool> m_muonPrinter{this, "MuonPrinter", "Rec::MuonPrintingTool/MuonPrintingTool"};
        ToolHandle<Trig::TrigDecisionTool> m_trigDec{this, "TrigDecTool", "Trig::TrigDecisionTool/TrigDecisionTool"};
        ToolHandle<Trk::ITrackSelectorTool> m_trackSelector{this, "TrackSelector", "InDet::InDetDetailedTrackSelectorTool/MuonCombinedInDetDetailedTrackSelectorTool"};
        ToolHandle<CP::IIsolationSelectionTool> m_isoTool{this, "IsoTool", ""};

        std::vector<std::string> m_selectMuonCategoriesStr;
        MuonPhysValMonitoringTool::MUCATEGORY getMuonSegmentTruthCategory(const xAOD::MuonSegment* truthMuSeg,
                                                                          const xAOD::TruthParticleContainer* muonTruthContainer);
        MuonPhysValMonitoringTool::MUCATEGORY getMuonTruthCategory(const xAOD::IParticle* prt);
        bool passesAcceptanceCuts(const xAOD::IParticle* prt);
        void SplitString(TString x, const TString& delim, std::vector<TString>& v);

        // Hists
        std::vector<std::unique_ptr<MuonValidationPlots>> m_muonValidationPlots;
        std::vector<std::unique_ptr<TriggerMuonValidationPlots>> m_TriggerMuonValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonMSTrackValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonMETrackValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonMSOnlyMETrackValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonIDTrackValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonIDSelectedTrackValidationPlots;
        std::vector<std::unique_ptr<MuonTrackValidationPlots>> m_muonIDForwardTrackValidationPlots;
        std::vector<std::unique_ptr<MuonSegmentValidationPlots>> m_muonSegmentValidationPlots;
        std::unique_ptr<Muon::RecoMuonPlotOrganizer> m_oUnmatchedRecoMuonPlots;
        std::unique_ptr<Muon::TruthMuonPlotOrganizer> m_oUnmatchedTruthMuonPlots;
        std::unique_ptr<Muon::RecoMuonTrackPlotOrganizer> m_oUnmatchedRecoMuonTrackPlots;
        std::unique_ptr<Muon::MuonSegmentPlots> m_oUnmatchedRecoMuonSegmentPlots;

        std::vector<std::unique_ptr<SlowMuonValidationPlots>> m_slowMuonValidationPlots;

        // overview hists
        std::vector<TH1F*> m_h_overview_nObjects;
        TH1F* m_h_overview_reco_category{nullptr};
        std::vector<TH1F*> m_h_overview_reco_authors;

        TH1F* m_h_overview_Z_mass{nullptr};
        TH1F* m_h_overview_Z_mass_ME{nullptr};
        TH1F* m_h_overview_Z_mass_ID{nullptr};

        std::vector<const xAOD::TruthParticle*> m_vMatchedTruthMuons;
        std::vector<const xAOD::Muon*> m_vMatchedMuons;
        std::vector<const xAOD::SlowMuon*> m_vMatchedSlowMuons;
        std::vector<const xAOD::TrackParticle*> m_vMatchedMuonTracks;
        std::vector<const xAOD::MuonSegment*> m_vMatchedMuonSegments;
        std::vector<const xAOD::TrackParticle*> m_vZmumuIDTracks;
        std::vector<const xAOD::TrackParticle*> m_vZmumuMETracks;
        std::vector<const xAOD::Muon*> m_vZmumuMuons;
        std::vector<const xAOD::Muon*> m_vEFMuons;
        std::vector<const xAOD::Muon*> m_vEFMuonsSelected;
        std::vector<const xAOD::L2StandAloneMuon*> m_vL2SAMuons;
        std::vector<const xAOD::L2StandAloneMuon*> m_vL2SAMuonsSelected;
        std::vector<const xAOD::L2CombinedMuon*> m_vL2CBMuons;
        std::vector<const xAOD::L2CombinedMuon*> m_vL2CBMuonsSelected;
        std::vector<const xAOD::Muon*> m_vRecoMuons;
        std::vector<const xAOD::Muon*> m_vRecoMuons_EffDen_CB;
        std::vector<const xAOD::Muon*> m_vRecoMuons_EffDen_MS;
        std::vector<const xAOD::Muon*> m_vRecoMuons_EffDen;

        template <class T> const T* getContainer(const std::string& containerName);
    };

    template <class T> const T* MuonPhysValMonitoringTool::getContainer(const std::string& containerName) {
        const T* ptr = evtStore()->retrieve<const T>(containerName);
        if (!ptr) { ATH_MSG_WARNING("Container '" << containerName << "' could not be retrieved"); }
        return ptr;
    }

    // I/O operators
    //////////////////////

    ///////////////////////////////////////////////////////////////////
    // Inline methods:
    ///////////////////////////////////////////////////////////////////
}  // namespace MuonPhysValMonitoring

#endif  //> !MUONPHYSVALMONITORING_MUONPHYSVALMONITORINGTOOL_H
