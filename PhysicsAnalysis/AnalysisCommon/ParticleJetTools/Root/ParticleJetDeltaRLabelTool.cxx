/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "ParticleJetTools/ParticleJetDeltaRLabelTool.h"
#include "ParticleJetTools/ParticleJetLabelCommon.h"
#include "xAODJet/JetContainer.h"
#include "AsgDataHandles/ReadHandle.h"
#include "AsgMessaging/Check.h"

using namespace std;
using namespace xAOD;

ParticleJetDeltaRLabelTool::ParticleJetDeltaRLabelTool(const std::string& name)
        : AsgTool(name) {
    declareProperty("LabelName", m_labelnames.singleint="", "Name of the jet label attribute to be added.");
    declareProperty("DoubleLabelName", m_labelnames.doubleint="", "Name of the jet label attribute to be added (with the possibility of up to 2 matched hadrons).");
    declareProperty("LabelPtName", m_labelnames.pt="HadronConeExclTruthLabelPt", "Name of attribute for labelling particle pt");
    declareProperty("LabelPtScaledName", m_labelnames.pt_scaled="HadronConeExclTruthLabelPtScaled", "Name of attribute for labelling particle pt divided by jet pt");
    declareProperty("LabelLxyName", m_labelnames.Lxy="HadronConeExclTruthLabelLxy", "Name of attribute for Lxy of labelling particle");
    declareProperty("LabelDRName", m_labelnames.dr="HadronConeExclTruthLabelDR", "Name of attribute for dR(part, jet) for labelling particle");
    declareProperty("LabelPdgIdName", m_labelnames.pdgId="HadronConeExclTruthLabelPdgId", "Name of attribute for pdgID of labelling particle");
    declareProperty("BLabelName", m_bottomlabelname="", "Name of the attribute to be added for matched B hadrons.");
    declareProperty("CLabelName", m_charmlabelname="", "Name of the attribute to be added for matched C hadrons.");
    declareProperty("TauLabelName", m_taulabelname="", "Name of the attribute to be added for matched taus.");
    declareProperty("PartPtMin", m_partptmin=5000, "Minimum pT of particles for labeling (MeV)");
    declareProperty("JetPtMin", m_jetptmin=10000, "Minimum pT of jets to be lebeled (MeV)");
    declareProperty("DRMax", m_drmax=0.3, "Maximum deltaR between a particle and jet to be labeled");
    declareProperty("MatchMode", m_matchmode="MinDR",
            "In the case that a particle matches two jets, the closest (MatchMode=MinDR) or highest-pT (MatchMode=MaxPt) jet will be labeled");
}

namespace {


}


StatusCode ParticleJetDeltaRLabelTool::initialize(){

  // initialize keys
  ATH_CHECK(m_tauPartCollectionKey.initialize());
  ATH_CHECK(m_bottomPartCollectionKey.initialize());
  ATH_CHECK(m_charmPartCollectionKey.initialize());

  // get the keys to the containers we just read, to build element
  // links later
  using Linker = ParticleJetTools::IParticleLinker;
  m_blinker.reset(new Linker(m_bottomPartCollectionKey, m_bottomlabelname));
  m_clinker.reset(new Linker(m_charmPartCollectionKey, m_charmlabelname));
  m_taulinker.reset(new Linker(m_tauPartCollectionKey, m_taulabelname));

  // build label decorators
  m_labeldecs.reset(new ParticleJetTools::LabelDecorators(m_labelnames));

  return StatusCode::SUCCESS;
}


StatusCode ParticleJetDeltaRLabelTool::decorate(const JetContainer& jets) const {

  ATH_MSG_VERBOSE("In " << name() << "::modify()");

    // Retrieve the particle and jet containers
    SG::ReadHandle<xAOD::TruthParticleContainer> truthtausReadHandle(m_tauPartCollectionKey);
    SG::ReadHandle<xAOD::TruthParticleContainer> truthbsReadHandle(m_bottomPartCollectionKey);
    SG::ReadHandle<xAOD::TruthParticleContainer> truthcsReadHandle(m_charmPartCollectionKey);

    if (!truthtausReadHandle.isValid()){
      ATH_MSG_DEBUG(" Invalid ReadHandle for xAOD::ParticleContainer with key: " << truthtausReadHandle.key());
      return StatusCode::FAILURE;
    }

    if (!truthbsReadHandle.isValid()){
      ATH_MSG_DEBUG(" Invalid ReadHandle for xAOD::ParticleContainer with key: " << truthbsReadHandle.key());
      return StatusCode::FAILURE;
    }

    if (!truthcsReadHandle.isValid()){
      ATH_MSG_DEBUG(" Invalid ReadHandle for xAOD::ParticleContainer with key: " << truthcsReadHandle.key());
      return StatusCode::FAILURE;
    }

    vector<vector<const TruthParticle*> > jetlabelpartsb = match(*truthbsReadHandle, jets);
    vector<vector<const TruthParticle*> > jetlabelpartsc = match(*truthcsReadHandle, jets);
    vector<vector<const TruthParticle*> > jetlabelpartstau = match(*truthtausReadHandle, jets);


    for (unsigned int iJet = 0; iJet < jets.size(); iJet++) {
        // remove children whose parent hadrons are also in the jet.
        // don't care about double tau jets
        // so leave them for now.

        using ParticleJetTools::childrenRemoved;
        childrenRemoved(jetlabelpartsb[iJet], jetlabelpartsb[iJet]);
        childrenRemoved(jetlabelpartsb[iJet], jetlabelpartsc[iJet]);
        childrenRemoved(jetlabelpartsc[iJet], jetlabelpartsc[iJet]);


        const Jet& jet = *jets.at(iJet);

        m_blinker->decorate(jet, jetlabelpartsb.at(iJet));
        m_clinker->decorate(jet, jetlabelpartsc.at(iJet));
        m_taulinker->decorate(jet, jetlabelpartstau.at(iJet));

        // set truth label to -99 for jets below pt threshold
        if (jet.pt() < m_jetptmin) {
            m_labeldecs->singleint(jet) = -99;
            m_labeldecs->doubleint(jet) = -99;
            continue;
        }

        // set truth label for jets above pt threshold
        // hierarchy: b > c > tau > light
        ParticleJetTools::Particles particles;
        particles.b = jetlabelpartsb.at(iJet);
        particles.c = jetlabelpartsc.at(iJet);
        particles.tau = jetlabelpartstau.at(iJet);
        ParticleJetTools::setJetLabels(jet, particles, *m_labeldecs);
    }

    return StatusCode::SUCCESS;
}


vector<vector<const TruthParticle*> >
ParticleJetDeltaRLabelTool::match(
        const TruthParticleContainer& parts,
        const JetContainer& jets) const {

  ATH_MSG_VERBOSE("In " << name() << "::match()");
  
    vector<vector<const TruthParticle*> > jetlabelparts(jets.size(), vector<const TruthParticle*>());

    // determine the match mode
    bool mindrmatch;
    if (m_matchmode == "MinDR")
        mindrmatch = true;
    else if (m_matchmode == "MaxPt")
        mindrmatch = false;
    else {
        ATH_MSG_FATAL("MatchMode must be MinDR or MaxPt");
        return jetlabelparts;
    }

    // loop over particles and find the best matched jet
    for (TruthParticleContainer::const_iterator part_itr = parts.begin();
            part_itr != parts.end(); ++part_itr) {

        const TruthParticle* part = *part_itr;

        if (part->pt() < m_partptmin)
            continue;

        double mindr = DBL_MAX;
        double maxpt = 0;
        int mindrjetidx = -1;
        int maxptjetidx = -1;
        for (unsigned int iJet = 0; iJet < jets.size(); iJet++) {

            const Jet& jet = *jets.at(iJet);

            double pt = jet.pt();
            if (pt < m_jetptmin)
                continue;

            double dr = jet.p4().DeltaR(part->p4());
            // too far for matching criterion
            if (dr > m_drmax)
                continue;

            // store the matched jet
            if (dr < mindr) {
                mindr = dr;
                mindrjetidx = iJet;
            }

            if (pt > maxpt) {
                maxpt = pt;
                maxptjetidx = iJet;
            }
        }

        // store the label particle with the jet
        if (mindrmatch && mindrjetidx >= 0)
            jetlabelparts.at(mindrjetidx).push_back(part);
        else if (maxptjetidx >= 0)
            jetlabelparts.at(maxptjetidx).push_back(part);
    }

    return jetlabelparts;
}

