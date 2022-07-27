/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonValidationPlots.h"

#include <utility>

#include "MuonHistUtils/MuonEnumDefs.h"

typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
typedef ElementLink<xAOD::MuonContainer> MuonLink;

MuonValidationPlots::MuonValidationPlots(PlotBase* pParent, const std::string& sDir, std::vector<int> wps,
                                         std::vector<unsigned int> authors, bool isData, bool doBinnedResolutionPlots,
                                         bool doSeparateSAFMuons, bool doMuonTree) :
    PlotBase(pParent, sDir),
    m_selectedWPs(std::move(wps)),
    m_selectedAuthors(std::move(authors)),
    m_truthSelections(2, ""),
    m_oTruthRelatedMuonPlots(nullptr),
    m_isData(isData),
    m_doSeparateSAFMuons(doSeparateSAFMuons),
    m_MuonTree(nullptr) {
    if (!m_isData) {
        m_truthSelections[0] = "all";           // no selection on truth muons (minimum selection is |eta|<2.5, pt>5 GeV, defined in
                                                // MuonPhysValMonitoringTool::handleTruthMuon()
        m_truthSelections[1] = "MSAcceptance";  // truth muons in MS acceptance (at least 4 associated hits in the MS)

        // histogram classes for all muons
        for (const auto& truthSelection : m_truthSelections) {
            m_oTruthMuonPlots.emplace_back(new Muon::TruthMuonPlotOrganizer(this, "truth/" + truthSelection));
        }
        m_oTruthRelatedMuonPlots.reset(
            new Muon::TruthRelatedMuonPlotOrganizer(this, "matched/AllMuons", doBinnedResolutionPlots));  
    }

    std::vector<int> allPlotCategories(0);
    std::vector<int> selectedPlotCategories(0);
    for (unsigned int i = 0; i < Muon::MAX_RECOPLOTCLASS; i++) {
        allPlotCategories.emplace_back(i);
        if (i != Muon::MUON_CHARGEPARAM) selectedPlotCategories.emplace_back(i);
    }

    // histogram classes for all muons
    m_oRecoMuonPlots.reset(new Muon::RecoMuonPlotOrganizer(this, "reco/AllMuons", allPlotCategories));

    // define a histogram class for each of the selected muon qualities
    for (unsigned int i = 0; i < m_selectedWPs.size(); i++) {
        std::string sQuality = Muon::EnumDefs::toString((xAOD::Muon::Quality)m_selectedWPs[i]);
        m_oRecoMuonPlots_perQuality.emplace_back(new Muon::RecoMuonPlotOrganizer(
            this, "reco/" + sQuality, (sQuality == "Medium" || sQuality == "Tight") ? allPlotCategories : selectedPlotCategories));

        if (!m_isData) {
            bool doBinnedPlots = false;
            if (sQuality == "Medium") doBinnedPlots = true;
            m_oTruthRelatedMuonPlots_perQuality.emplace_back(
                new Muon::TruthRelatedMuonPlotOrganizer(this, "matched/" + sQuality, doBinnedPlots));
        }
    }

    // define a histogram class for each of the selected muon authors (+one inclusive for all authors)
    for (unsigned int i = 0; i < m_selectedAuthors.size(); i++) {
        std::string sAuthor = Muon::EnumDefs::toString((xAOD::Muon::Author)m_selectedAuthors[i]);
        if (sAuthor == "CaloTag") sAuthor = "CaloTagTight";
        m_oRecoMuonPlots_perAuthor.emplace_back(
            new Muon::RecoMuonPlotOrganizer(this, "reco/" + sAuthor, (sAuthor == "MuidCo") ? allPlotCategories : selectedPlotCategories));
        if (!m_isData)
            m_oTruthRelatedMuonPlots_perAuthor.emplace_back(
                new Muon::TruthRelatedMuonPlotOrganizer(this, "matched/" + sAuthor, doBinnedResolutionPlots));
    }

    // define histogram class for loose CaloTag and append to author plots, not very nice workaround though
    for (unsigned int i = 0; i < m_selectedAuthors.size(); i++) {
        if ((xAOD::Muon::Author)m_selectedAuthors[i] == xAOD::Muon::CaloTag) {  // found CaloTag in list, also do CaloTagLoose
            m_oRecoMuonPlots_perAuthor.emplace_back(new Muon::RecoMuonPlotOrganizer(this, "reco/CaloTagLoose", selectedPlotCategories));
            if (!m_isData)
                m_oTruthRelatedMuonPlots_perAuthor.emplace_back(
                    new Muon::TruthRelatedMuonPlotOrganizer(this, "matched/CaloTagLoose", doBinnedResolutionPlots));
        }
    }

    // define histogram class for SiliconAssociatedForwardMuons
    if (m_doSeparateSAFMuons) {
        m_oRecoMuonPlots_SiAssocFwrdMu.emplace_back(new Muon::RecoMuonPlotOrganizer(this, "reco/SiAssocForward", selectedPlotCategories));
        if (!m_isData)
            m_oTruthRelatedMuonPlots_SiAssocFwrdMu.emplace_back(
                new Muon::TruthRelatedMuonPlotOrganizer(this, "matched/SiAssocForward", doBinnedResolutionPlots));
    }

    if (doMuonTree) m_MuonTree.reset(new Muon::MuonTree(this, "", !m_isData));
}

MuonValidationPlots::~MuonValidationPlots() =default;

void MuonValidationPlots::fillRecoMuonPlots(const xAOD::Muon& mu, float weight) {
    // fill hists for all muons
    m_oRecoMuonPlots->fill(mu, weight);
    // fill separate hists for each muon quality
    xAOD::Muon::Quality muqual = mu.quality();
    for (unsigned int i = 0; i < m_selectedWPs.size(); i++) {
        if (muqual <= (xAOD::Muon::Quality)m_selectedWPs[i]) { m_oRecoMuonPlots_perQuality[i]->fill(mu, weight); }
    }
    // fill separate hists for each author
    for (unsigned int i = 0; i < m_selectedAuthors.size(); i++) {
        if (mu.isAuthor((xAOD::Muon::Author)m_selectedAuthors[i])) {
            if ((xAOD::Muon::Author)m_selectedAuthors[i] == (xAOD::Muon::CaloTag)) {
                int ipar = 0;
                if (mu.parameter(ipar, xAOD::Muon::CaloMuonIDTag)) { ; }
                if (ipar < 11) continue;
            }

            // filter SiliconAssociatedForwardMuons
            if (mu.muonType() != (xAOD::Muon::MuonType)xAOD::Muon::SiliconAssociatedForwardMuon || !m_doSeparateSAFMuons)
                m_oRecoMuonPlots_perAuthor[i]->fill(mu, weight);
        }
    }
    // fill SiliconAssociatedForwardMuons
    for (unsigned int i = 0; i < m_oTruthRelatedMuonPlots_SiAssocFwrdMu.size(); i++) {
        if (mu.muonType() == (xAOD::Muon::MuonType)xAOD::Muon::SiliconAssociatedForwardMuon)
            m_oRecoMuonPlots_SiAssocFwrdMu[i]->fill(mu, weight);
    }
    // fill CaloTagLoose (one additional plot in plot list)
    unsigned int counter = m_selectedAuthors.size();
    if (counter + 1 == m_oRecoMuonPlots_perAuthor.size()) {
        if (mu.isAuthor(xAOD::Muon::CaloTag)) m_oRecoMuonPlots_perAuthor[counter]->fill(mu, weight);
    }
}

void MuonValidationPlots::fillTruthMuonPlots(const xAOD::TruthParticle& truthMu, float weight) {
    m_oTruthMuonPlots[0]->fill(truthMu, weight);  // no selections
    if (isGoodTruthTrack(truthMu)) {  // in MS acceptance (minimum precision hits)
        m_oTruthMuonPlots[1]->fill(truthMu, weight);
    }
}

void MuonValidationPlots::fill(const xAOD::Muon& mu, float weight) { fillRecoMuonPlots(mu, weight); }

void MuonValidationPlots::fill(const xAOD::TruthParticle& truthMu, float weight) { fillTruthMuonPlots(truthMu, weight); }

void MuonValidationPlots::fill(const xAOD::TruthParticle* truthMu, const xAOD::Muon* mu, const xAOD::TrackParticleContainer* MSTracks,
                               float weight) {
    if (truthMu) fillTruthMuonPlots(*truthMu, weight);
    if (mu) fillRecoMuonPlots(*mu, weight);

    if ((mu) && (truthMu)) {
        // plots for all
        m_oTruthRelatedMuonPlots->fill(*truthMu, *mu, MSTracks, weight);
        // fill SiliconAssociatedForwardMuons
        for (unsigned int i = 0; i < m_oTruthRelatedMuonPlots_SiAssocFwrdMu.size(); i++) {
            if (mu->muonType() == xAOD::Muon::SiliconAssociatedForwardMuon || !m_doSeparateSAFMuons)
                m_oTruthRelatedMuonPlots_SiAssocFwrdMu[i]->fill(*truthMu, *mu, MSTracks, weight);
        }

        // plots per quality
        xAOD::Muon::Quality muqual = mu->quality();
        for (unsigned int i = 0; i < m_selectedWPs.size(); i++) {
            if (muqual <= (xAOD::Muon::Quality)m_selectedWPs[i]) {
                m_oTruthRelatedMuonPlots_perQuality[i]->fill(*truthMu, *mu, MSTracks, weight);
            }
        }
        // plots per author
        for (unsigned int i = 0; i < m_selectedAuthors.size(); i++) {
            if (mu->isAuthor((xAOD::Muon::Author)m_selectedAuthors[i])) {
                if ((xAOD::Muon::Author)m_selectedAuthors[i] == xAOD::Muon::CaloTag) {
                    int ipar = 0;
                    if (mu->parameter(ipar, xAOD::Muon::CaloMuonIDTag)) { ; }
                    if (ipar < 11) continue;
                }
                // filter SilicionAssociatedForwardMuons
                if (mu->muonType() != xAOD::Muon::SiliconAssociatedForwardMuon || !m_doSeparateSAFMuons)
                    m_oTruthRelatedMuonPlots_perAuthor[i]->fill(*truthMu, *mu, MSTracks, weight);
            }
        }
        // fill CaloTagLoose (one additional plot in plot list)
        unsigned int counter = m_selectedAuthors.size();
        if (counter + 1 == m_oRecoMuonPlots_perAuthor.size()) {
            if (mu->isAuthor(xAOD::Muon::CaloTag)) m_oTruthRelatedMuonPlots_perAuthor[counter]->fill(*truthMu, *mu, MSTracks, weight);
        }
    }
}

Muon::MuonTree* MuonValidationPlots::getMuonTree() { return m_MuonTree.get(); }

void MuonValidationPlots::fillTreeBranches(const xAOD::Muon& mu) {
    if (m_MuonTree) m_MuonTree->fillRecoMuonBranches(mu);
}

void MuonValidationPlots::fillTreeBranches(const xAOD::TruthParticle& truthMu) {
    if (m_MuonTree) m_MuonTree->fillTruthMuonBranches(truthMu, isGoodTruthTrack(truthMu));
}

void MuonValidationPlots::fillTreeBranches(const xAOD::TruthParticle* truthMu, const xAOD::Muon* mu,
                                           const xAOD::TrackParticleContainer* MSTracks) {
    if (!m_MuonTree) return;
    if (mu) m_MuonTree->fillRecoMuonBranches(*mu);
    if (truthMu) m_MuonTree->fillTruthMuonBranches(*truthMu, isGoodTruthTrack(*truthMu));

    if (truthMu && mu) { m_MuonTree->fillTruthMuonBranches(*truthMu, *mu, MSTracks, isGoodTruthTrack(*truthMu)); }
}

void MuonValidationPlots::fillTree(const xAOD::EventInfo* eventInfo, bool isData) {
    if (!m_MuonTree) return;
    m_MuonTree->fillEventBranches(eventInfo, isData);
    m_MuonTree->getTree()->Fill();
    m_MuonTree->postFillTreeActions();
}

bool MuonValidationPlots::isGoodTruthTrack(const xAOD::TruthParticle& truthMu) {
    static const std::array<std::string,6> hitTypes{"innerSmallHits",  "innerLargeHits", "middleSmallHits",
                               "middleLargeHits", "outerSmallHits", "outerLargeHits"};  // MDT + CSC
    int minPrecHits = 5;

    int nPrecHits = 0;
    bool hasEnoughPrecHits = false;

    for (const auto& hitTypeItr : hitTypes) {
        if (truthMu.isAvailable<uint8_t>(hitTypeItr)) {
            nPrecHits += truthMu.auxdata<uint8_t>(hitTypeItr);
            if (nPrecHits >= minPrecHits) {
                hasEnoughPrecHits = true;
                break;
            }
        }
    }
    return (hasEnoughPrecHits);
}
