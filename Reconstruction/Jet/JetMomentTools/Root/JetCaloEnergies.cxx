/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "JetMomentTools/JetCaloEnergies.h"

#include "xAODJet/JetAccessorMap.h"

#include "JetUtils/JetCaloQualityUtils.h"
#include "PFlowUtils/FEHelpers.h"

#include "xAODCaloEvent/CaloCluster.h"
#include "xAODPFlow/PFO.h"
#include "xAODPFlow/FlowElement.h"

#include <vector>

//**********************************************************************

JetCaloEnergies::JetCaloEnergies(const std::string& name)
: AsgTool(name) { }

//**********************************************************************

StatusCode JetCaloEnergies::initialize() {
  ATH_MSG_INFO("Initializing JetCaloEnergies " << name());
  
  if(m_jetContainerName.empty()){
    ATH_MSG_ERROR("JetCaloEnergies needs to have its input jet container configured!");
    return StatusCode::FAILURE;
  }
  
  m_ePerSamplingKey = m_jetContainerName + "." + m_ePerSamplingKey.key();
  m_emFracKey = m_jetContainerName + "." + m_emFracKey.key();
  m_hecFracKey = m_jetContainerName + "." + m_hecFracKey.key();
  m_psFracKey = m_jetContainerName + "." + m_psFracKey.key();

  if(m_calcClusterBasedVars){
    m_ePerSamplingClusterKey = m_jetContainerName + "." + m_ePerSamplingClusterKey.key();
    m_emFracClusterKey = m_jetContainerName + "." + m_emFracClusterKey.key();
    m_hecFracClusterKey = m_jetContainerName + "." + m_hecFracClusterKey.key();
    m_psFracClusterKey = m_jetContainerName + "." + m_psFracClusterKey.key();
  }

  ATH_CHECK(m_ePerSamplingKey.initialize());
  ATH_CHECK(m_emFracKey.initialize());
  ATH_CHECK(m_hecFracKey.initialize());
  ATH_CHECK(m_psFracKey.initialize());

  ATH_CHECK(m_ePerSamplingClusterKey.initialize(m_calcClusterBasedVars));
  ATH_CHECK(m_emFracClusterKey.initialize(m_calcClusterBasedVars));
  ATH_CHECK(m_hecFracClusterKey.initialize(m_calcClusterBasedVars));
  ATH_CHECK(m_psFracClusterKey.initialize(m_calcClusterBasedVars));

  return StatusCode::SUCCESS;
}

//**********************************************************************

StatusCode JetCaloEnergies::decorate(const xAOD::JetContainer& jets) const {
  ATH_MSG_VERBOSE("Begin decorating jets.");
  for(const xAOD::Jet* jet : jets) {
    SG::WriteDecorHandle<xAOD::JetContainer, std::vector<float> > ePerSamplingHandle(m_ePerSamplingKey);
    size_t numConstit = jet->numConstituents();
    ePerSamplingHandle(*jet) = std::vector<float>(CaloSampling::Unknown, 0.);
    std::vector<float>& ePerSampling = ePerSamplingHandle(*jet);
    for ( float& e : ePerSampling ) e = 0.0; // re-initialize

    if ( numConstit == 0 ) {
      ATH_MSG_VERBOSE("Jet has no constituents.");
      continue;
    }
    
    // should find a more robust solution than using 1st constit type.
    xAOD::Type::ObjectType ctype = jet->rawConstituent( 0 )->type();
    if ( ctype  == xAOD::Type::CaloCluster ) {
      ATH_MSG_VERBOSE("  Constituents are calo clusters.");
      fillEperSamplingCluster(*jet, ePerSampling);

    } else if (ctype  == xAOD::Type::ParticleFlow) {
      ATH_MSG_VERBOSE("  Constituents are pflow objects.");
      fillEperSamplingPFO(*jet, ePerSampling);

    } else if (ctype  == xAOD::Type::FlowElement) {
      ATH_MSG_VERBOSE("  Constituents are FlowElements.");
      fillEperSamplingFE(*jet, ePerSampling);

      // In addition, calculate variables using the underlying cluster rather than
      // the energy-subtracted FlowElements (improved implementation)
      if(m_calcClusterBasedVars){
	ATH_MSG_VERBOSE("  Constituents are FlowElements - Additional calculation");

	SG::WriteDecorHandle<xAOD::JetContainer, std::vector<float> > ePerSamplingClusterHandle(m_ePerSamplingClusterKey);
	ePerSamplingClusterHandle(*jet) = std::vector<float>(CaloSampling::Unknown, 0.);
	std::vector<float>& ePerSamplingCluster = ePerSamplingClusterHandle(*jet);
	for ( float& e : ePerSamplingCluster ) e = 0.0; // re-initialize

	fillEperSamplingFEClusterBased(*jet, ePerSamplingCluster);
      }

    }else {
      ATH_MSG_VERBOSE("Constituents are not CaloClusters, PFOs, or FlowElements.");
    }

  }
  return StatusCode::SUCCESS;
}

void JetCaloEnergies::fillEperSamplingCluster(const xAOD::Jet& jet, std::vector<float> & ePerSampling ) const {
  // loop over raw constituents
  size_t numConstit = jet.numConstituents();    
  for ( size_t i=0; i<numConstit; i++ ) {
    if(jet.rawConstituent(i)->type()!=xAOD::Type::CaloCluster) {
      ATH_MSG_WARNING("Tried to call fillEperSamplingCluster with a jet constituent that is not a cluster!");
      continue;
    }
    const xAOD::CaloCluster* constit = static_cast<const xAOD::CaloCluster*>(jet.rawConstituent(i));      
    for ( size_t s= CaloSampling::PreSamplerB; s< CaloSampling::Unknown; s++ ) {
      ePerSampling[s] += constit->eSample( (xAOD::CaloCluster::CaloSample) s );
    }
  }
  SG::WriteDecorHandle<xAOD::JetContainer, float> emFracHandle(m_emFracKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> hecFracHandle(m_hecFracKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> psFracHandle(m_psFracKey);
  
  emFracHandle(jet) = jet::JetCaloQualityUtils::emFraction( ePerSampling );
  hecFracHandle(jet) = jet::JetCaloQualityUtils::hecF( &jet );
  psFracHandle(jet) = jet::JetCaloQualityUtils::presamplerFraction( &jet );
}

#define FillESamplingPFO( LAYERNAME )                                        \
  float E_##LAYERNAME = 0.0;                                                \
  if (constit->attribute(xAOD::PFODetails::eflowRec_LAYERENERGY_##LAYERNAME, E_##LAYERNAME)) { \
    ePerSampling[CaloSampling::LAYERNAME] += E_##LAYERNAME;                \
  }
  
void JetCaloEnergies::fillEperSamplingPFO(const xAOD::Jet & jet, std::vector<float> & ePerSampling ) const {

  float emTot=0;
  float hecTot=0;
  float eTot =0;
  size_t numConstit = jet.numConstituents();

  for ( size_t i=0; i<numConstit; i++ ) {
    if (jet.rawConstituent(i)->type()==xAOD::Type::ParticleFlow){
      const xAOD::PFO* constit = static_cast<const xAOD::PFO*>(jet.rawConstituent(i));
      if ( fabs(constit->charge())>FLT_MIN ){
        eTot += constit->track(0)->e();
      } else {
        eTot += constit->e();
        FillESamplingPFO(PreSamplerB);
        FillESamplingPFO(EMB1);
        FillESamplingPFO(EMB2);
        FillESamplingPFO(EMB3);

        FillESamplingPFO(PreSamplerE);
        FillESamplingPFO(EME1);
        FillESamplingPFO(EME2);
        FillESamplingPFO(EME3);

        FillESamplingPFO(HEC0);
        FillESamplingPFO(HEC1);
        FillESamplingPFO(HEC2);
        FillESamplingPFO(HEC3);

        FillESamplingPFO(TileBar0);
        FillESamplingPFO(TileBar1);
        FillESamplingPFO(TileBar2);

        FillESamplingPFO(TileGap1);
        FillESamplingPFO(TileGap2);
        FillESamplingPFO(TileGap3);

        FillESamplingPFO(TileExt0);
        FillESamplingPFO(TileExt1);
        FillESamplingPFO(TileExt2);

        FillESamplingPFO(FCAL0);
        FillESamplingPFO(FCAL1);
        FillESamplingPFO(FCAL2);

        FillESamplingPFO(MINIFCAL0);
        FillESamplingPFO(MINIFCAL1);
        FillESamplingPFO(MINIFCAL2);
        FillESamplingPFO(MINIFCAL3);        

        emTot += ( E_PreSamplerB+E_EMB1+E_EMB2+E_EMB3+
                  E_PreSamplerE+E_EME1+E_EME2+E_EME3+
                  E_FCAL0 );
        hecTot += ( E_HEC0+E_HEC1+E_HEC2+E_HEC3 );
        
      }//only consider neutral PFO
    } else {
      ATH_MSG_WARNING("Tried to call fillEperSamplingPFlow with a jet constituent that is not a PFO!");
    }
  }

  SG::WriteDecorHandle<xAOD::JetContainer, float> emFracHandle(m_emFracKey);
  if(eTot != 0.0){
    emFracHandle(jet) = emTot/eTot; 
    /*
     * Ratio of EM layer calorimeter energy of neutrals to sum of all constituents 
     * at EM scale (note charged PFO have an EM scale at track scale, and charged weights are ignored)
     * */
  }
  else {
    emFracHandle(jet)  = 0.;
  }

  SG::WriteDecorHandle<xAOD::JetContainer, float> hecFracHandle(m_hecFracKey);
  if (eTot != 0.0){
    hecFracHandle(jet) = hecTot/eTot;
  }
  else{
    hecFracHandle(jet) = 0.;
  }
  
}

// R21 way of calculating energy-per-layer using directly the FE decorations
// In R21, the link from PFOs to clusters was broken and thus energy after charged subtraction was used
void JetCaloEnergies::fillEperSamplingFE(const xAOD::Jet& jet, std::vector<float> & ePerSampling ) const {
  float emTot = 0.;
  float hecTot = 0.;
  float psTot = 0.;
  float eTot = 0.;
  size_t numConstit = jet.numConstituents();    

  std::vector<int> clusterIndices;

  for ( size_t i=0; i<numConstit; i++ ) {
    if(jet.rawConstituent(i)->type()!=xAOD::Type::FlowElement) {
      ATH_MSG_WARNING("Tried to call fillEperSamplingFE with a jet constituent that is not a FlowElement!");
      continue;
    }
    const xAOD::FlowElement* constit = static_cast<const xAOD::FlowElement*>(jet.rawConstituent(i));

    // Need to distinguish two cases:
    // (1) Jet is a PFlow jet (constituents are charged or neutral FEs) or
    // (2) Jet is a UFO jet (need to get the underlying charged and neutral FEs first)

    //For PFlow jets, we can directly get the information from the constituent
    if(constit->signalType() & xAOD::FlowElement::PFlow){

      //Charged FlowElements:
      if(constit->isCharged()){
	eTot += constit->chargedObject(0)->e();
      }
      //Neutral FlowElements
      else{
	eTot += constit->e();
	//Get the energy-per-layer information from the FE, not the underlying cluster (i.e. after subtraction)
	std::vector<float> constitEPerSampling = FEHelpers::getEnergiesPerSampling(*constit);
	for ( size_t s = CaloSampling::PreSamplerB; s < CaloSampling::Unknown; s++ ) {
	  ePerSampling[s] += constitEPerSampling[s];
	}
	emTot += (constitEPerSampling[CaloSampling::PreSamplerB] + constitEPerSampling[CaloSampling::EMB1]
		   + constitEPerSampling[CaloSampling::EMB2]        + constitEPerSampling[CaloSampling::EMB3]
		   + constitEPerSampling[CaloSampling::PreSamplerE] + constitEPerSampling[CaloSampling::EME1]
		   + constitEPerSampling[CaloSampling::EME2]        + constitEPerSampling[CaloSampling::EME3]
		   + constitEPerSampling[CaloSampling::FCAL0]);

	hecTot += (constitEPerSampling[CaloSampling::HEC0] + constitEPerSampling[CaloSampling::HEC1]
		   + constitEPerSampling[CaloSampling::HEC2] + constitEPerSampling[CaloSampling::HEC3]);

	psTot += (constitEPerSampling[CaloSampling::PreSamplerB] + constitEPerSampling[CaloSampling::PreSamplerE]);
      }
    }
    else{
      //For UFO jets, we first need to get the charged and neutral FE + corresponding energy from combined UFOs

      // UFO is simply a charged FlowElement
      if(constit->signalType() == xAOD::FlowElement::Charged){
	eTot += constit->chargedObject(0)->e();
      }
      //UFO is simply a neutral Flowelement
      else if(constit->signalType() == xAOD::FlowElement::Neutral){
	// For neutral UFOs, there is only one "other object" stored which is the neutral FE
	// Protection in case there is something wrong with this FE
	if(constit->otherObjects().size() != 1 || !constit->otherObject(0)){
	  continue;
	}

	// Cast other object as FlowElement
	const xAOD::FlowElement* nFE = static_cast<const xAOD::FlowElement*>(constit->otherObject(0));

	eTot += nFE->e();

	std::vector<float> neutralEPerSampling = FEHelpers::getEnergiesPerSampling(*nFE);
        for ( size_t s = CaloSampling::PreSamplerB; s < CaloSampling::Unknown; s++ ) {
          ePerSampling[s] += neutralEPerSampling[s];
	}
        emTot += (neutralEPerSampling[CaloSampling::PreSamplerB] + neutralEPerSampling[CaloSampling::EMB1]
                   + neutralEPerSampling[CaloSampling::EMB2]        + neutralEPerSampling[CaloSampling::EMB3]
                   + neutralEPerSampling[CaloSampling::PreSamplerE] + neutralEPerSampling[CaloSampling::EME1]
                   + neutralEPerSampling[CaloSampling::EME2]        + neutralEPerSampling[CaloSampling::EME3]
                   + neutralEPerSampling[CaloSampling::FCAL0]);

        hecTot += (neutralEPerSampling[CaloSampling::HEC0] + neutralEPerSampling[CaloSampling::HEC1]
                   + neutralEPerSampling[CaloSampling::HEC2] + neutralEPerSampling[CaloSampling::HEC3]);

	psTot += (neutralEPerSampling[CaloSampling::PreSamplerB] + neutralEPerSampling[CaloSampling::PreSamplerE]);
      }
      else if(constit->signalType() == xAOD::FlowElement::Combined){
	// For the combined UFOs, otherObjects are neutral FEs
	for (size_t n = 0; n < constit->otherObjects().size(); ++n) {
	  if(! constit->otherObject(n)) continue;
	  const xAOD::FlowElement* nFE_from_combined = static_cast<const xAOD::FlowElement*>(constit->otherObject(n));
	  //One neutral FE can be matched to various tracks and therefore be used for several UFOs
	  //We do not want to double count the energy and only add it once
	  if(std::find(clusterIndices.begin(), clusterIndices.end(), nFE_from_combined->index()) == clusterIndices.end()){
	    eTot += nFE_from_combined->e();
	    std::vector<float> neutralFromCombEPerSampling = FEHelpers::getEnergiesPerSampling(*nFE_from_combined);
	    for ( size_t s = CaloSampling::PreSamplerB; s < CaloSampling::Unknown; s++ ) {
	      ePerSampling[s] += neutralFromCombEPerSampling[s];
	    }
	    emTot += (neutralFromCombEPerSampling[CaloSampling::PreSamplerB] + neutralFromCombEPerSampling[CaloSampling::EMB1]
		      + neutralFromCombEPerSampling[CaloSampling::EMB2]        + neutralFromCombEPerSampling[CaloSampling::EMB3]
		      + neutralFromCombEPerSampling[CaloSampling::PreSamplerE] + neutralFromCombEPerSampling[CaloSampling::EME1]
		      + neutralFromCombEPerSampling[CaloSampling::EME2]        + neutralFromCombEPerSampling[CaloSampling::EME3]
		      + neutralFromCombEPerSampling[CaloSampling::FCAL0]);

	    hecTot += (neutralFromCombEPerSampling[CaloSampling::HEC0] + neutralFromCombEPerSampling[CaloSampling::HEC1]
		       + neutralFromCombEPerSampling[CaloSampling::HEC2] + neutralFromCombEPerSampling[CaloSampling::HEC3]);

	    psTot += (neutralFromCombEPerSampling[CaloSampling::PreSamplerB] + neutralFromCombEPerSampling[CaloSampling::PreSamplerE]);

	    clusterIndices.push_back(nFE_from_combined->index());
	  }
	}
      }
    }
  }

  SG::WriteDecorHandle<xAOD::JetContainer, float> emFracHandle(m_emFracKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> hecFracHandle(m_hecFracKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> psFracHandle(m_psFracKey);
  
  emFracHandle(jet)  = eTot != 0. ? emTot/eTot  : 0.;
  hecFracHandle(jet) = eTot != 0. ? hecTot/eTot : 0.;
  psFracHandle(jet)  = eTot != 0. ? psTot/eTot  : 0.;
}


// New way of calculating energy per layer for FE-based jets
// The underlying cluster is used to calculate the variables rather than the energy-subtracted FE

void JetCaloEnergies::fillEperSamplingFEClusterBased(const xAOD::Jet& jet, std::vector<float> & ePerSampling ) const {
  float emTot = 0.;
  float hecTot = 0.;
  float psTot = 0.;
  float eTot = 0.;
  size_t numConstit = jet.numConstituents();
  std::unique_ptr<std::vector<const xAOD::CaloCluster*> > constitV_tot = std::unique_ptr<std::vector<const xAOD::CaloCluster*>>(new std::vector<const xAOD::CaloCluster*>);

  for ( size_t i=0; i<numConstit; i++ ) {
    if(jet.rawConstituent(i)->type()!=xAOD::Type::FlowElement) {
      ATH_MSG_WARNING("Tried to call fillEperSamplingFE with a jet constituent that is not a FlowElement!");
      continue;
    }
    const xAOD::FlowElement* constit = static_cast<const xAOD::FlowElement*>(jet.rawConstituent(i));

    for (size_t n = 0; n < constit->otherObjects().size(); ++n) {
      if(! constit->otherObject(n)) continue;
      int index_pfo = constit->otherObject(n)->index();
      if(index_pfo<0) continue;

      const auto* fe = (constit->otherObject(n));
      const xAOD::CaloCluster* cluster = nullptr;

      //If we have a cluster, we can directly access the calorimeter information
      if(fe->type() == xAOD::Type::CaloCluster){
        cluster = dynamic_cast<const xAOD::CaloCluster*> (fe);
      }
      //If we have a PFO, we should still get the associated cluster first
      else {
        const xAOD::FlowElement* pfo = dynamic_cast<const xAOD::FlowElement*>(fe);
        if(pfo->otherObjects().size() > 0 && pfo->otherObject(0) && pfo->otherObject(0)->type() == xAOD::Type::CaloCluster){
          cluster = dynamic_cast<const xAOD::CaloCluster*> (pfo->otherObject(0));
        }
      }
      if(!cluster) continue;

      if(std::find(constitV_tot->begin(), constitV_tot->end(), cluster) == constitV_tot->end()){
        for ( size_t s= CaloSampling::PreSamplerB; s< CaloSampling::Unknown; s++ ) {
          ePerSampling[s] += cluster->eSample( (xAOD::CaloCluster::CaloSample) s );
        }
        eTot += cluster->rawE();

        emTot += (cluster->eSample( CaloSampling::PreSamplerB) + cluster->eSample( CaloSampling::EMB1)
		  + cluster->eSample( CaloSampling::EMB2)        + cluster->eSample( CaloSampling::EMB3)
		  + cluster->eSample( CaloSampling::PreSamplerE) + cluster->eSample( CaloSampling::EME1)
		  + cluster->eSample( CaloSampling::EME2)        + cluster->eSample( CaloSampling::EME3)
		  + cluster->eSample( CaloSampling::FCAL0));

        hecTot += (cluster->eSample( CaloSampling::HEC0) + cluster->eSample( CaloSampling::HEC1)
		   + cluster->eSample( CaloSampling::HEC2) + cluster->eSample( CaloSampling::HEC3));

        psTot += (cluster->eSample( CaloSampling::PreSamplerB) + cluster->eSample( CaloSampling::PreSamplerE));

        constitV_tot->push_back(cluster);
      }
    }
  }

  SG::WriteDecorHandle<xAOD::JetContainer, float> emFracClusterHandle(m_emFracClusterKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> hecFracClusterHandle(m_hecFracClusterKey);
  SG::WriteDecorHandle<xAOD::JetContainer, float> psFracClusterHandle(m_psFracClusterKey);

  emFracClusterHandle(jet)  = eTot != 0. ? emTot/eTot  : 0.;
  hecFracClusterHandle(jet) = eTot != 0. ? hecTot/eTot : 0.;
  psFracClusterHandle(jet)  = eTot != 0. ? psTot/eTot  : 0.;
}
