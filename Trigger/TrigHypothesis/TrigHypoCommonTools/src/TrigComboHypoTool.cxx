/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigComboHypoTool.h"

#include <cmath>

using namespace TrigCompositeUtils;

TrigComboHypoTool::TrigComboHypoTool(const std::string& type, 
				     const std::string& name, 
				     const IInterface*  parent)
  : ComboHypoToolBase(type, name, parent)
{}


StatusCode TrigComboHypoTool::initialize()
{
  ATH_MSG_DEBUG("AcceptAll  = " << m_acceptAll );
  ATH_MSG_DEBUG("Variable   = " << m_varTag );
  ATH_MSG_DEBUG("UseCut min = " << m_useMin );
  ATH_MSG_DEBUG("UseCut max = " << m_useMax );
  ATH_MSG_DEBUG("varCut min = " << m_varMin );
  ATH_MSG_DEBUG("varCut max = " << m_varMax );
  ATH_MSG_DEBUG("LegA       = " << m_legA );
  ATH_MSG_DEBUG("LegB       = " << m_legB );

  if ( not m_monTool.name().empty() ) {
    ATH_CHECK( m_monTool.retrieve() );
    ATH_MSG_DEBUG("m_monTool name: " << m_monTool);
  }

  if (m_legA==""){
    ATH_MSG_ERROR("LegA not set!");
    return StatusCode::FAILURE;
  }
  if (m_legB==""){
    ATH_MSG_ERROR("LegB not set!");
    return StatusCode::FAILURE;
  }
  if ((!m_useMin) && (!m_useMax)){
    ATH_MSG_ERROR("Trying to configure the Tool without setting UseMin and UseMax!");
    return StatusCode::FAILURE;
  }
  

  ATH_MSG_DEBUG("Initialization completed successfully");

  return StatusCode::SUCCESS;
}

bool TrigComboHypoTool::executeAlg(std::vector<LegDecision> &combination) const {
  ATH_MSG_DEBUG("On combination executeAlg");
  auto varOfAccepted  = Monitored::Scalar( m_varTag+"OfAccepted"   , -1.0 );
  auto varOfProcessed = Monitored::Scalar( m_varTag+"OfProcessed"  , -1.0 );
  auto monitorIt      = Monitored::Group( m_monTool, varOfAccepted, varOfProcessed);

  //check that we found the two legs
  int nCombs(combination.size());
  if (nCombs < 2){
    ATH_MSG_ERROR("Number of legs found is less than 2! N_legs = " << combination.size() );
    return false;
  }
  int           legA_index(-1), legB_index(-1);

  ATH_MSG_DEBUG("Legs available = "<< combination);
  for (int i=0; i<nCombs; ++i){
    auto combId = HLT::Identifier(combination[i].first);
    if (!TrigCompositeUtils::isLegId(combId))
      continue;
    std::string   legName = combId.name().substr(0,6);
    if (legName == m_legA){
      legA_index = i;
    }else  if (legName == m_legB){
      legB_index = i;
    }
    ATH_MSG_DEBUG("\t Leg: "<< legName <<", full name:"<<combId.name());
  }

  if ( legA_index<0){
    ATH_MSG_ERROR("legA = "<< m_legA << " NOT FOUND!");
    return false;
  }
  if ( legB_index<0){
    ATH_MSG_ERROR("legB = "<< m_legB << " NOT FOUND!");
    return false;
  }

  auto EL= combination[legA_index].second;    
  auto legA_pLink = TrigCompositeUtils::findLink<xAOD::IParticleContainer>( *EL, featureString() ).link;
  if (!legA_pLink.isValid()){
    ATH_MSG_ERROR("link for "<<m_legA<<" not valid");
    return false;
  }
  ATH_MSG_DEBUG("link for legA: "<<m_legA<<" is valid");

  EL = combination[legB_index].second;
  auto legB_pLink = TrigCompositeUtils::findLink<xAOD::IParticleContainer>( *EL, featureString() ).link;
  if (!legB_pLink.isValid()){
    ATH_MSG_ERROR("link for "<<m_legB<<" not valid");
    return false;
  }
  ATH_MSG_DEBUG("link for legB: "<<m_legB<<" is valid");

  TLorentzVector hlv1 = (*legA_pLink)->p4();
  TLorentzVector hlv2 = (*legB_pLink)->p4();  

  // apply the cut
  bool  pass(true);
  float value(-9999.);

  //should we make a switch? (if this list of observables is used only here probably not...)
  std::array<std::string, 2> valid_varTags = {"dR","invm"};
  if(m_varTag ==  valid_varTags[0]) {
    value =  hlv1.DeltaR(hlv2);
  }else if (m_varTag == valid_varTags[1]){
    TLorentzVector hlvtot = hlv1+hlv2;
    value = hlvtot.M();
  }else {
    ATH_MSG_ERROR("m_varTag =  "<<m_varTag<<" not present in the list of valid_tags : " << valid_varTags);
    return false;
  }
  varOfProcessed = value;

  ATH_MSG_DEBUG("Found a combination with " << varOfProcessed);

  if (m_useMin && m_useMax){
    if (varOfProcessed < m_varMin || varOfProcessed > m_varMax){ 
      ATH_MSG_DEBUG("Combination failed var cut: "<< m_varTag <<"= "<< varOfProcessed << " not in [" << m_varMin << "," <<  m_varMax << "]");
      pass=false;
    }else{
      varOfAccepted = value;
      ATH_MSG_DEBUG( m_varTag <<"= "<< varOfAccepted << " is  within [" <<m_varMin<< "," << m_varMax << "] This selection passed! ");
    }
  }else if (m_useMin){
    if (varOfProcessed < m_varMin ){ 
      ATH_MSG_DEBUG("Combination failed var cut: "<< m_varTag <<"= "<< varOfProcessed << " not > " << m_varMin);
      pass=false;
    }else{
      varOfAccepted = value;
      ATH_MSG_DEBUG( m_varTag <<"= "<< varOfAccepted << " < " <<m_varMin << " This selection passed! ");
    }
  }else if (m_useMax){
    if (varOfProcessed > m_varMax){ 
      ATH_MSG_DEBUG("Combination failed var cut: "<< m_varTag <<"= "<< varOfProcessed << " not < " << m_varMax);
      pass=false;
    }else{
      varOfAccepted = value;
      ATH_MSG_DEBUG( m_varTag <<"= "<< varOfAccepted << " > " << m_varMax << " This selection passed! ");
    }
  }
  return pass;

}




