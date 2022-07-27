/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

// ********************************************************************
//
// NAME:     TrigCaloTowerMaker.cxx
// PACKAGE:  Trigger/TrigAlgorithms/TrigCaloRec
//
// AUTHOR:   P.A. Delsart
//           This is an Hlt algorithm that creates Calorimeter towers
//           It assumes CaloCellContainers are already produced and accessable
//           through the input TriggerElement. Based on TrigCaloRec
//
// ********************************************************************
//
//
#include <sstream>
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ListItem.h"

#include "AthenaMonitoringKernel/Monitored.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "CaloEvent/CaloCellContainer.h"
#include "EventKernel/INavigable4Momentum.h"
#include "NavFourMom/INavigable4MomentumCollection.h"

#include "CaloUtils/CaloTowerBuilderToolBase.h"
#include "CaloEvent/CaloTowerSeg.h"
#include "CaloEvent/CaloTowerContainer.h"

#include "TrigCaloTowerMaker.h"

class ISvcLocator;

/////////////////////////////////////////////////////////////////////
// CONSTRUCTOR:
/////////////////////////////////////////////////////////////////////
//
TrigCaloTowerMaker::TrigCaloTowerMaker(const std::string& name, ISvcLocator* pSvcLocator)
  : AthAlgorithm(name, pSvcLocator)
  , m_nEtaTowers(50)
  , m_nPhiTowers(64)
  , m_minEta(-2.5)
  , m_maxEta(2.5)
  , m_deta(0.5)
  , m_dphi(0.5)
  , m_includeFcal(false)
{


  // Name of Tower Maker Tools
  declareProperty( "TowerMakerTools",   m_towerMakerNames );
  declareProperty( "NumberOfEtaTowers", m_nEtaTowers );
  declareProperty( "NumberOfPhiTowers", m_nPhiTowers );
  declareProperty( "EtaMin", m_minEta );
  declareProperty( "EtaMax", m_maxEta );
  declareProperty( "DeltaEta", m_deta );
  declareProperty( "DeltaPhi", m_dphi );



}

/////////////////////////////////////////////////////////////////////
    // INITIALIZE:
    // The initialize method will create all the required algorithm objects
    // Note that it is NOT NECESSARY to run the initialize of individual
    // sub-algorithms.  The framework takes care of it.
    /////////////////////////////////////////////////////////////////////


    StatusCode TrigCaloTowerMaker::initialize()
{
  ATH_MSG_DEBUG ("in initialize()");

  if (!m_monTool.empty()) {
    ATH_MSG_DEBUG("Retrieving monTool");
    CHECK(m_monTool.retrieve());
  } else {
    ATH_MSG_INFO("No monTool configured => NO MONITOING");
  }

  // Cache pointer to ToolSvc
  IToolSvc* toolSvc = 0;// Pointer to Tool Service
  if (service("ToolSvc", toolSvc).isFailure()) {
    ATH_MSG_FATAL ( " Tool Service not found " );
    return StatusCode::FAILURE;
  }

  std::vector<std::string>::iterator itrName;
  std::vector<std::string>::iterator endName;

  itrName = m_towerMakerNames.begin();
  endName = m_towerMakerNames.end();

  CaloTowerSeg theTowerSeg(m_nEtaTowers,m_nPhiTowers,m_minEta,m_maxEta);
  
  for (; itrName!=endName; ++itrName) {
    
    ListItem theItem(*itrName);
    IAlgTool* algtool;
    if ( (*itrName).find("FC") != std::string::npos ) m_includeFcal=true;
    
    if( toolSvc->retrieveTool(theItem.type(), theItem.name(), algtool,this).isFailure() ) {
      ATH_MSG_FATAL ( "Unable to find tool for " << (*itrName) );
      return StatusCode::FAILURE;
    } else {
      CaloTowerBuilderToolBase* mytool = dynamic_cast<CaloTowerBuilderToolBase*>(algtool);
      if(mytool == NULL) {
	ATH_MSG_FATAL("dynamic cast to CaloTowerBuilderToolBase failed for" << theItem.name());
	return StatusCode::FAILURE;
      }
      else {
	ATH_MSG_DEBUG ( " set towerSet for this tool " << theItem.name() );
	mytool->setTowerSeg(theTowerSeg);
	m_towerMakerPointers.push_back(mytool);
      }
    }
  }

  // end of helpers...
 
  //  if (timerSvc()) {
  //m_totalTimer     = addTimer("TimeTot");
  //m_towerContTimer = addTimer("TCRRTowerContainer");
  //}
  m_caloTowerNav4LinkKey = m_outputTowerKey.key();

ATH_CHECK( m_inputRoiKey.initialize(SG::AllowEmpty) );
ATH_CHECK( m_inputCellsKey.initialize() );
ATH_CHECK( m_outputTowerKey.initialize() );
ATH_CHECK( m_caloTowerNav4LinkKey.initialize() );

  ATH_MSG_DEBUG
        ( "Initialization of TrigCaloTowerMaker completed successfully"
        );

  return StatusCode::SUCCESS;
}



StatusCode TrigCaloTowerMaker::execute()
{
  // Monitoring initialization...
  auto timer = Monitored::Timer("TIME_execute");
  auto time_tools = Monitored::Timer("TIME_tools");
  auto mon_towerContainerSize = Monitored::Scalar("towerContainerSize", 0.);

  auto monitorIt = Monitored::Group( m_monTool, timer, time_tools, mon_towerContainerSize);

    ATH_MSG_DEBUG ( "in execute()" );
    

  auto ctx = getContext();
    
  // Get RoiDescriptor
  const bool seedLess = m_inputRoiKey.empty();
  std::unique_ptr<TrigRoiDescriptor> roiDescriptor;
  if (!seedLess){
      SG::ReadHandle<TrigRoiDescriptorCollection> roiCollection = SG::makeHandle( m_inputRoiKey, ctx);
      if ( not roiCollection.isValid() ) {
        ATH_MSG_ERROR("Cell maker did not get a valid RoIs collection");
        return StatusCode::FAILURE;
      }
      if (roiCollection->size()==0) {
        ATH_MSG_DEBUG(" RoI collection size = 0");
        return StatusCode::SUCCESS;
      }
      if (roiCollection->size() > 1) ATH_MSG_WARNING("Misconfiguration - Called with " << roiCollection->size() << " ROI, but it should be called with 1 RoI - Will only process the first RoI"); 

      roiDescriptor = std::make_unique<TrigRoiDescriptor>(*roiCollection->begin());
      ATH_MSG_DEBUG("Operating on " << roiCollection->size() <<"RoI(s)");
  }
  else {
    roiDescriptor = std::make_unique<TrigRoiDescriptor>(TrigRoiDescriptor(true));
  }
  ATH_MSG_DEBUG ( " RoI id " << roiDescriptor->roiId()
	  << " located at   phi = " <<  roiDescriptor->phi()
	  << ", eta = " << roiDescriptor->eta() );
  
  auto  caloTowerContainer =   SG::makeHandle (m_outputTowerKey, ctx); 

  /// ho hum, this needs a flag for using own wdiths rather than those from the roiDescriptor  
  /// in addition, this will *not* work properly for composite RoIs
    
  if ( roiDescriptor->composite() ) { 
    ATH_MSG_WARNING ( " Attempting to use composite RoI as a normal RoI - this is probably *not* what you want to do "
	  << *roiDescriptor );
  }
  
  double   eta0 = roiDescriptor->eta();
  double   phi0 = roiDescriptor->phi();
      
  while (phi0 > 2.*M_PI) phi0 -= 2. * M_PI;
  while (phi0 < 0. )     phi0 += 2. * M_PI;
  double etamin = eta0-m_deta/2.;
  double etamax = eta0+m_deta/2.;
  // What if phimin < 0 or phimax > 2Pi ??  question
  double phimin = phi0 - m_dphi/2.;
  double phimax = phi0 + m_dphi/2.;

  ATH_MSG_DEBUG(" eta0 = "<< eta0 << " phi0 = "<< phi0 << 
		" etamin = "<< etamin << " etamax = "<< etamax << 
		" phimin = "<< phimin << " phimax = "<< phimax);
  
  
  CaloTowerSeg myseg (m_nEtaTowers,m_nPhiTowers,m_minEta,m_maxEta);
  CaloTowerSeg::SubSeg subseg = myseg.subseg ( roiDescriptor->eta(), m_deta,  roiDescriptor->phi(), m_dphi );
  if ( m_includeFcal ) {
    ATH_CHECK( caloTowerContainer .record (std::make_unique<CaloTowerContainer>(myseg)));
  } else {
    ATH_CHECK( caloTowerContainer .record (std::make_unique<CaloTowerContainer>( subseg.segmentation() )));
  }
  CaloTowerContainer* pCaloTowerContainer=caloTowerContainer .ptr();
  
  ATH_CHECK( caloTowerContainer.symLink (m_caloTowerNav4LinkKey));
  ATH_MSG_DEBUG( "CaloTowerContainer" << caloTowerContainer.name() << " symlinked to  INavigable4MomentumCollection in StoreGate");

  //  if (timerSvc()) m_towerContTimer->start();  // Measures the time to retrieve the cells in the RoI  

  auto caloCellContainer =   SG::makeHandle (m_inputCellsKey, ctx); 
 
        
  // Get the last container in the vector. Should be th one produced by the previous TrigCaloCellMaker.
  const CaloCellContainer* theCellCont = caloCellContainer.ptr();

  ATH_MSG_DEBUG(" REGTEST: Retrieved a Cell Container of Size= " << theCellCont->size());

  // Now Build the towers -----------------------------------------------------------------

  if(!m_towerMakerNames.empty()) {
    
    time_tools.start();
    std::vector<CaloTowerBuilderToolBase*>::const_iterator itrtwr=m_towerMakerPointers.begin();
    std::vector<CaloTowerBuilderToolBase*>::const_iterator endtwr=m_towerMakerPointers.end();
    int index=0;
    for (; itrtwr!=endtwr; ++itrtwr) {
      //if((*itrtwr)->setProperty( StringProperty("CellContainerName",cellCollKey) ).isFailure()) {
	//ATH_MSG_ERROR ( "ERROR setting the CaloCellContainer name in the offline tool" );
	//return HLT::TOOL_FAILURE;
      //}
      if ( m_includeFcal ) {
        if ( (*itrtwr)->execute(ctx, pCaloTowerContainer, theCellCont).isFailure() ) {
	  ATH_MSG_ERROR("Error executing tool " << m_towerMakerNames[index]);
        } else {
	  ATH_MSG_DEBUG("Executed tool " << m_towerMakerNames[index]);
        }
      } else {
        if ( (*itrtwr)->execute(ctx, pCaloTowerContainer,theCellCont, &subseg).isFailure() ) {
	  ATH_MSG_ERROR("Error executing tool " << m_towerMakerNames[index]);
        } else {	
	  ATH_MSG_DEBUG("Executed tool " << m_towerMakerNames[index]);
        }
      }
      ++index;
      time_tools.stop();

    }

    ATH_MSG_DEBUG(" REGTEST: Produced a Tower Container " << caloTowerContainer.name() << " at " << caloTowerContainer.cptr() << " of Size= " << pCaloTowerContainer->size());
    ATH_MSG_DEBUG(" handle info:  " << caloTowerContainer);
     mon_towerContainerSize = (float)pCaloTowerContainer->size() ;
  } 

  return StatusCode::SUCCESS;
}
