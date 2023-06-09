/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// local include(s)
#include "tauRecTools/MvaTESEvaluator.h"
#include "tauRecTools/HelperFunctions.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>

// tools include(s) 

//#include "TauAnalysisTools/HelperFunctions.h"

//_____________________________________________________________________________
MvaTESEvaluator::MvaTESEvaluator(const std::string& name)
  : TauRecToolBase(name)
  , reader(0)
  , mu(0)
  , nVtxPU(0)    
  , center_lambda(0)
  , first_eng_dens(0)
  , second_lambda(0)
  , presampler_frac(0)
  , em_probability(0)    
  , ptCombined(0)
  , ptLC_D_ptCombined(0)
  , ptConstituent_D_ptCombined(0)
  , etaConstituent(0)    
  , PanTauBDT_1p0n_vs_1p1n(0)
  , PanTauBDT_1p1n_vs_1pXn(0)
  , PanTauBDT_3p0n_vs_3pXn(0)
  , nTracks(0)
  , PFOEngRelDiff(0)    
  , truthPtVis(0)
  , pt(0)
  , ptPanTauCellBased(0)
  , ptDetectorAxis(0)
  , truthDecayMode(0)
  , PanTau_DecayMode(0)
{
  declareProperty( "WeightFileName", m_sWeightFileName = "MvaTES_20170207_v2_BDTG.weights.root" );
  //R20.7 pi0-fix version: "MvaTES_20161015_pi0fix_BDTG.weights.root"
}

//_____________________________________________________________________________
MvaTESEvaluator::~MvaTESEvaluator()
{
}

//_____________________________________________________________________________
StatusCode MvaTESEvaluator::initialize(){

  // Declare input variables to the reader
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.mu", &mu) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.nVtxPU", &nVtxPU) );
  
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ClustersMeanCenterLambda", &center_lambda) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ClustersMeanFirstEngDens", &first_eng_dens) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ClustersMeanSecondLambda", &second_lambda) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ClustersMeanPresamplerFrac", &presampler_frac) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ClustersMeanEMProbability", &em_probability) );
  
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.pt_combined", &ptCombined) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ptDetectorAxis/TauJetsAuxDyn.pt_combined", &ptLC_D_ptCombined) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ptPanTauCellBased/TauJetsAuxDyn.pt_combined", &ptConstituent_D_ptCombined) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.etaPanTauCellBased", &etaConstituent) );
  
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.PanTau_BDTValue_1p0n_vs_1p1n", &PanTauBDT_1p0n_vs_1p1n) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.PanTau_BDTValue_1p1n_vs_1pXn", &PanTauBDT_1p1n_vs_1pXn) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.PanTau_BDTValue_3p0n_vs_3pXn", &PanTauBDT_3p0n_vs_3pXn) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.nTracks", &nTracks) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.PFOEngRelDiff", &PFOEngRelDiff) );

  // Spectator variables declared in the training have to be added to the reader, too
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.truthPtVis", &truthPtVis) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.pt", &pt) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ptPanTauCellBased", &ptPanTauCellBased) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.ptDetectorAxis", &ptDetectorAxis) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.truthDecayMode", &truthDecayMode) );
  m_availableVars.insert( std::make_pair("TauJetsAuxDyn.PanTau_DecayMode", &PanTau_DecayMode) );
  
  std::string weightFile = find_file(m_sWeightFileName);

  reader = tauRecTools::configureMVABDT( m_availableVars, weightFile.c_str() );
  if(reader==0) {
    ATH_MSG_FATAL("Couldn't configure MVA");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;

}

//_____________________________________________________________________________
StatusCode MvaTESEvaluator::eventInitialize()
{
  // HACK HACK HACK: Get nVtxPU, AuxElement::ConstAccessor doesn't work
  nVtxPU = 0;
  if(evtStore()->contains<xAOD::VertexContainer>("PrimaryVertices")){
    ATH_CHECK(evtStore()->retrieve(m_xVertexContainer, "PrimaryVertices"));  
    for (auto xVertex : *m_xVertexContainer)
    if (xVertex->vertexType() == xAOD::VxType::PileUp)
      nVtxPU++;
  }
  else {
    ATH_MSG_WARNING("No xAOD::VertexContainer, setting nVtxPU to 0");
    nVtxPU=0;
  }

  return StatusCode::SUCCESS;
}

//_____________________________________________________________________________
StatusCode MvaTESEvaluator::execute(xAOD::TauJet& xTau){

  // Retrieve input variables
  
  // Retrieve event info
  static SG::AuxElement::ConstAccessor<double> acc_mu("mu");
  //static SG::AuxElement::ConstAccessor<int> acc_nVtxPU("nVtxPU");
  mu = acc_mu(xTau);
  //nVtxPU = acc_nVtxPU(xTau);

  // Retrieve seed jet info
  xTau.detail(xAOD::TauJetParameters::ClustersMeanCenterLambda, center_lambda);
  xTau.detail(xAOD::TauJetParameters::ClustersMeanFirstEngDens, first_eng_dens);
  xTau.detail(xAOD::TauJetParameters::ClustersMeanEMProbability,em_probability);
  xTau.detail(xAOD::TauJetParameters::ClustersMeanSecondLambda, second_lambda);
  xTau.detail(xAOD::TauJetParameters::ClustersMeanPresamplerFrac, presampler_frac);

  // Retrieve pantau and LC-precalib TES
  etaConstituent = xTau.etaPanTauCellBased();
  float ptLC = xTau.ptDetectorAxis();
  float ptConstituent = xTau.ptPanTauCellBased();
  static SG::AuxElement::ConstAccessor<float> acc_pt_combined("pt_combined");
  ptCombined = acc_pt_combined(xTau);
  ptLC_D_ptCombined = ptLC / ptCombined;
  ptConstituent_D_ptCombined = ptConstituent / ptCombined;
  
  // Retrieve substructures info
  static SG::AuxElement::ConstAccessor<float> acc_PanTauBDT_1p0n_vs_1p1n("PanTau_BDTValue_1p0n_vs_1p1n");
  static SG::AuxElement::ConstAccessor<float> acc_PanTauBDT_1p1n_vs_1pXn("PanTau_BDTValue_1p1n_vs_1pXn");
  static SG::AuxElement::ConstAccessor<float> acc_PanTauBDT_3p0n_vs_3pXn("PanTau_BDTValue_3p0n_vs_3pXn");
  PanTauBDT_1p0n_vs_1p1n = acc_PanTauBDT_1p0n_vs_1p1n(xTau);
  PanTauBDT_1p1n_vs_1pXn = acc_PanTauBDT_1p1n_vs_1pXn(xTau);
  PanTauBDT_3p0n_vs_3pXn = acc_PanTauBDT_3p0n_vs_3pXn(xTau);
  nTracks = (float)xTau.nTracks();
  xTau.detail(xAOD::TauJetParameters::PFOEngRelDiff, PFOEngRelDiff);

  float ptMVA = float( ptCombined * reader->GetResponse() );
  if(ptMVA<1) ptMVA=1;
  xTau.setP4(xAOD::TauJetParameters::FinalCalib, ptMVA, etaConstituent, xTau.phiPanTauCellBased(), 0);

  return StatusCode::SUCCESS;

}
