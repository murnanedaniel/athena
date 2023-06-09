// -*- C++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**************************************************************************
 **
 **   File: Trigger/TrigHypothesis/TrigBPhysHypo/TrigL2BMuMuXHypo.cxx
 **
 **   Description: L2 hypothesis algorithms for B+/- -> K+/- mu+ mu-
 **                Create TrigL2Bphys
 **
 **   Author: Pavel.Strachota@cern.ch
 **
 **   Created:   10.01.2010 based on C. Adorisio's code TrigL2DiMuHypo.cxx & D. Scheirich's code TrigL2LbMuMuLHypo.cxx 
 **
 **   Modified:  26.09.2012 : "B_c -> D_s(*) (->Phi(->K+K-)) mu+mu-" added (Leonid Gladilin <gladilin@mail.cern.ch>)
 **
 **************************************************************************/

#include "TrigL2BMuMuXHypo.h"

//#include "EventInfo/EventInfo.h"
//#include "EventInfo/EventID.h"

//#include "TrigMuonEvent/TrigCombDiMuonContainer.h"
//#include "TrigSteeringEvent/TrigRoiDescriptor.h"

//#include "TrigParticle/TrigL2Bphys.h"
//#include "TrigParticle/TrigL2BphysContainer.h"
//#include "TrigInDetEvent/TrigVertexCollection.h"
//#include "TrigSteeringEvent/TrigPassBits.h"
#include "xAODTrigger/TrigPassBits.h"
#include "TrigNavigation/Navigation.h"

// additions of xAOD objects
#include "TrigBphysHelperUtilsTool.h"
// additions of xAOD objects
//#include "xAODEventInfo/EventInfo.h"
#include "xAODTrigBphys/TrigBphys.h"
#include "xAODTrigBphys/TrigBphysContainer.h"

/*------------------------------------------------------------------------------------------*/
TrigL2BMuMuXHypo::TrigL2BMuMuXHypo(const std::string & name, ISvcLocator* pSvcLocator):
  HLT::HypoAlgo(name, pSvcLocator)
,m_bphysHelperTool("TrigBphysHelperUtilsTool")
{

  // Read cuts
  declareProperty("AcceptAll", m_acceptAll = false);
          // false
  declareProperty("LowerBMassCut", m_lowerBMassCut = 4500.);                    //4500.00
  declareProperty("UpperBMassCut", m_upperBMassCut = 5900.);               //5900.00
  declareProperty("LowerBdMassCut", m_lowerBdMassCut = 4600.);                  //4600.00
  declareProperty("UpperBdMassCut", m_upperBdMassCut = 6000.);             //6000.00
  declareProperty("LowerKstarMassCut", m_lowerKstarMassCut = 600.);                  //4600.00
  declareProperty("UpperKstarMassCut", m_upperKstarMassCut = 1500.);             //6000.00
  declareProperty("LowerBsMassCut", m_lowerBsMassCut = 4600.);                  //4600.00
  declareProperty("UpperBsMassCut", m_upperBsMassCut = 6000.);             //6000.00
  declareProperty("LowerPhi1020MassCut", m_lowerPhi1020MassCut = 990.);                  //4600.00
  declareProperty("UpperPhi1020MassCut", m_upperPhi1020MassCut = 1050.);             //6000.00
  declareProperty("LowerLbMassCut", m_lowerLbMassCut = 4600.);                  //4600.00
  declareProperty("UpperLbMassCut", m_upperLbMassCut = 6400.);             //6400.00
  declareProperty("LowerLambdaMassCut", m_lowerLambdaMassCut = 1040.);                  //4600.00
  declareProperty("UpperLambdaMassCut", m_upperLambdaMassCut = 1220.);             //6400.00

  declareProperty("LowerBcMassCut", m_lowerBcMassCut = 1800.);             // 1800.00
  declareProperty("UpperBcMassCut", m_upperBcMassCut = 7050.);             // 7050.00
  declareProperty("LowerDsMassCut", m_lowerDsMassCut = 1600.);             // 1600.00
  declareProperty("UpperDsMassCut", m_upperDsMassCut = 2400.);             // 2400.00

  // zero counters
  m_lastEvent = 999;
  m_lastEventPassed = 999;
  m_countTotalEvents=0;
  m_countTotalRoI=0;
  m_countPassedEvents=0;
  m_countPassedRoIs=0;
  m_countPassedBMass=0;
  m_countPassedBdMass=0;
  m_countPassedBsMass=0;
  m_countPassedLbMass=0;
  m_countPassedBcMass=0;

}
/*------------------------------------------------------------------------------------------*/
TrigL2BMuMuXHypo::~TrigL2BMuMuXHypo()
{ }
/*------------------------------------------------------------------------------------------*/
HLT::ErrorCode TrigL2BMuMuXHypo::hltInitialize()
{

  if(msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "In Initialize " << endmsg;
/*    msg() << MSG::INFO << "AcceptAll            = " << (m_acceptAll==true ? "True" : "False") << endmsg;
  msg() << MSG::INFO << "B+/- mass cuts: "  << m_lowerBMassCut  << " < Mass(K+/- MuMu) < " << m_upperBMassCut << endmsg;
  msg() << MSG::INFO << "Bd mass cuts: " << m_lowerBdMassCut << " < Mass(Kstar MuMu) < " << m_upperBdMassCut << endmsg;
  msg() << MSG::INFO << "Bs mass cuts: " << m_lowerBsMassCut << " < Mass(Phi MuMu) < " << m_upperBsMassCut << endmsg;
  msg() << MSG::INFO << "Lb mass cuts: " << m_lowerLbMassCut << " < Mass(Lambda MuMu) < " << m_upperLbMassCut << endmsg;
  msg() << MSG::INFO << "Bc mass cuts: " << m_lowerBcMassCut << " < Mass(Ds MuMu) < " << m_upperBcMassCut << endmsg;
}
*/
    if (m_bphysHelperTool.retrieve().isFailure()) {
        msg() << MSG::ERROR << "Can't find TrigBphysHelperUtilsTool" << endmsg;
        return HLT::BAD_JOB_SETUP;
    } else {
        if (msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "TrigBphysHelperUtilsTool found" << endmsg;
    }
    

  return HLT::OK;
}
/*------------------------------------------------------------------------------------------*/
HLT::ErrorCode TrigL2BMuMuXHypo::hltFinalize()
{
  msg() << MSG::INFO << "In Finalize()" << endmsg;
  msg() << MSG::INFO << "|----------------------- SUMMARY FROM TrigL2BMuMuXHypo -------------|" << endmsg;
  msg() << MSG::INFO << "Run on events/RoIs " << m_countTotalEvents << "/" << m_countTotalRoI <<  endmsg;
  msg() << MSG::INFO << "Passed events/RoIs " << m_countPassedEvents << "/" << m_countPassedRoIs <<  endmsg;
  msg() << MSG::INFO << "RoIs Passed BMass: " << m_countPassedBMass << endmsg;
  msg() << MSG::INFO << "RoIs Passed BdMass: " << m_countPassedBdMass << endmsg;
  msg() << MSG::INFO << "RoIs Passed BsMass: " << m_countPassedBsMass << endmsg;
  msg() << MSG::INFO << "RoIs Passed LbMass: " << m_countPassedLbMass << endmsg;
  msg() << MSG::INFO << "|-----------------------------------------------------------------------|" << endmsg;

  return HLT::OK;
}

//-----------------------------------------------------------------------------------------------
HLT::ErrorCode TrigL2BMuMuXHypo::hltExecute(const HLT::TriggerElement* outputTE, bool& pass)
{
  if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "In Execute()" << endmsg;

  bool PassedBplus=false;
  bool PassedBd=false;
  bool PassedBs=false;
  bool PassedLb=false;
  bool PassedBc=false;
  bool result = false;

    // Retrieve event info
    int IdRun   = 0;
    int IdEvent = 0;
    
    // event info
    uint32_t runNumber(0), evtNumber(0), lbBlock(0);
    if (m_bphysHelperTool->getRunEvtLb( runNumber, evtNumber, lbBlock).isFailure()) {
        msg() << MSG::ERROR << "Error retriving EventInfo" << endmsg;
    }
    IdRun = runNumber;
    IdEvent = evtNumber;

    if (IdEvent != (int) m_lastEvent) {
    m_countTotalEvents++;
    m_lastEvent=IdEvent;
  }
  m_countTotalRoI++;

  // Accept-All mode: temporary patch; should be done with force-accept
  if (m_acceptAll) {
    if ( msgLvl() <= MSG::DEBUG )
      msg() << MSG::DEBUG << "AcceptAll property is set: taking all events" << endmsg;
    pass = true;
    return HLT::OK;
  } else {
    if ( msgLvl() <= MSG::DEBUG )
      msg() << MSG::DEBUG << "AcceptAll property not set: applying selection" << endmsg;
  }

  if ( msgLvl() <= MSG::DEBUG ){
    msg() << MSG::DEBUG
        << "Using outputTE ( " << outputTE << " ) ->getId(): " << outputTE->getId() << endmsg;
  }


  // create vector for TrigL2Bphys particles
    // const TrigL2BphysContainer* trigBphysColl_b = 0;
    const xAOD::TrigBphysContainer* trigBphysColl_b(nullptr);

  // bplus part
  HLT::ErrorCode status = getFeature(outputTE, trigBphysColl_b, "L2BMuMuXFex");

  if ( status != HLT::OK ) {
    if ( msgLvl() <= MSG::WARNING) {
      msg() << MSG::WARNING << "Failed to get TrigBphysics collection" << endmsg;
    }
    return HLT::OK;
  }

  if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "Retrieved Bphys collection  trigBphysColl_b = " << trigBphysColl_b << endmsg;
  if ( trigBphysColl_b == 0 ) {
    if ( msgLvl() <= MSG::DEBUG )
      msg() << MSG::DEBUG << "No Bphys particles to analyse, leaving!" << endmsg;
    return HLT::OK;
  }

  if ( msgLvl() <= MSG::DEBUG ) {
    msg() << MSG::DEBUG << "Got TrigBphys collection with " << trigBphysColl_b->size() << " TrigBphys particles " << endmsg;
  }

  // if no Bphys particles were found, just leave TrigBphysColl. empty and leave
  if ( trigBphysColl_b->size() == 0 ) {
    if ( msgLvl() <= MSG::DEBUG )
      msg() << MSG::DEBUG << "No Bphys particles to analyse, leaving!" << endmsg;
    return HLT::OK;
  }

  //TrigPassBits *bits = HLT::makeTrigPassBits(trigBphysColl_b);
  std::unique_ptr<xAOD::TrigPassBits> xBits = xAOD::makeTrigPassBits<xAOD::TrigBphysContainer>(trigBphysColl_b);

  // now loop over Bphys particles to see if one passes cuts
  for ( xAOD::TrigBphysContainer::const_iterator bphysIter = trigBphysColl_b->begin(); bphysIter !=  trigBphysColl_b->end(); ++bphysIter) {
    
    //determine decay mode
    xAOD::TrigBphys::pType decayType = (*bphysIter)->particleType();
    std::string decayName("Unknown");
    if(decayType == xAOD::TrigBphys::BKMUMU)      decayName = "B+ -> mu mu K+";
    if(decayType == xAOD::TrigBphys::BDKSTMUMU)   decayName = "Bd -> mu mu K*";
    if(decayType == xAOD::TrigBphys::BSPHIMUMU)   decayName = "Bs -> mu mu Phi";
    if(decayType == xAOD::TrigBphys::LBLMUMU)     decayName = "Lambda_b -> mu mu Lambda";
    if(decayType == xAOD::TrigBphys::BCDSMUMU)    decayName = "Bc -> mu mu Ds";
    
    if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "Bphys particle type " << decayName << ", " << decayType << " with mass " << (*bphysIter)->mass() << endmsg;


    if (decayType == xAOD::TrigBphys::BKMUMU    || decayType == xAOD::TrigBphys::BDKSTMUMU ||
        decayType == xAOD::TrigBphys::BSPHIMUMU || decayType == xAOD::TrigBphys::LBLMUMU   ||
        decayType == xAOD::TrigBphys::BCDSMUMU ) {

      bool thisPassedBplus=false;
      bool thisPassedBd=false;
      bool thisPassedBs=false;
      bool thisPassedLb=false;
      bool thisPassedBc=false;

      float BMass = (*bphysIter)->mass();
      if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "BMass = " << BMass << endmsg;
      
      if(decayType == xAOD::TrigBphys::BKMUMU) {
        if ( BMass > m_lowerBMassCut && BMass < m_upperBMassCut ) {
          if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << " Mass =  " << BMass  << " -> B+ Mass passed " << endmsg;
          PassedBplus=true;
          thisPassedBplus=true;
        }
      } else { // Bd, Bs, Lb, Bc
        const xAOD::TrigBphys* trigPartX = (*bphysIter)->secondaryDecay();
        if(!trigPartX) {
          if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "No secondary decay pointer in Bphys particle of type " << (*bphysIter)->particleType() << endmsg;
        } else {
          float XMass = trigPartX->mass();
          if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "XMass = " << XMass << endmsg;
          
          if(decayType == xAOD::TrigBphys::BDKSTMUMU) {
            if ( BMass > m_lowerBdMassCut && BMass < m_upperBdMassCut ) {
              if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << " BMass =  " << BMass  << " -> Bd Mass passed " << endmsg;
              if( XMass > m_lowerKstarMassCut && XMass < m_upperKstarMassCut ) {
                if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "  XMass =  " << XMass  << " -> K* Mass passed " << endmsg;
                PassedBd=true;
                thisPassedBd=true;
              }
            }
          }
          if(decayType == xAOD::TrigBphys::BSPHIMUMU) {
            if ( BMass > m_lowerBsMassCut && BMass < m_upperBsMassCut ) {
              if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << " BMass =  " << BMass  << " -> Bs Mass passed " << endmsg;
              if( XMass > m_lowerPhi1020MassCut && XMass < m_upperPhi1020MassCut ) {
                if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "  XMass =  " << XMass  << " -> Phi Mass passed " << endmsg;
                PassedBs=true;
                thisPassedBs=true;
              }
            }
          }
          if(decayType == xAOD::TrigBphys::LBLMUMU) {
            if ( BMass > m_lowerLbMassCut && BMass < m_upperLbMassCut ) {
              if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << " BMass =  " << BMass  << " -> Lb Mass passed " << endmsg;
              if( XMass > m_lowerLambdaMassCut && XMass < m_upperLambdaMassCut ) {
                if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "  XMass =  " << XMass  << " -> Lambda Mass passed " << endmsg;
                PassedLb=true;
                thisPassedLb=true;
              }
            }
          }
          if(decayType == xAOD::TrigBphys::BCDSMUMU) {
            if ( BMass > m_lowerBcMassCut && BMass < m_upperBcMassCut ) {
              if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << " BMass =  " << BMass  << " -> Bc Mass passed " << endmsg;
              if( XMass > m_lowerDsMassCut && XMass < m_upperDsMassCut ) {
                if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "  XMass =  " << XMass  << " -> Ds Mass passed " << endmsg;
                PassedBc=true;
                thisPassedBc=true;
              }
            }
          }
        } // end if trigPartX
      } // end Bd, Bs, Lb, Bc
      
      if ( thisPassedBplus || thisPassedBd || thisPassedBs || thisPassedLb || thisPassedBc ) {
        result=true;
        //HLT::markPassing(bits, *bphysIter, trigBphysColl_b);
        xBits->markPassing((*bphysIter),trigBphysColl_b,true);
      }
    }
  }

  if (PassedBplus) m_countPassedBMass++;
  if (PassedBd) m_countPassedBdMass++;
  if (PassedBs) m_countPassedBsMass++;
  if (PassedLb) m_countPassedLbMass++;
  if (PassedBc) m_countPassedBcMass++;
  if (result) {
    m_countPassedRoIs++;
    if (IdEvent!= (int) m_lastEventPassed) {
      m_countPassedEvents++;
      m_lastEventPassed=IdEvent;
    }
  }

  // Reach this point successfully
  if ( msgLvl() <= MSG::DEBUG) msg() << MSG::DEBUG << "Run: " << IdRun << " Event: " << IdEvent << "  result is " << result << endmsg;


  if (!m_acceptAll) pass = result;

  // store result
  //if ( attachBits(outputTE, bits) != HLT::OK ) {
  //  msg() << MSG::ERROR << "Problem attaching TrigPassBits! " << endmsg;
  //}
  if(attachFeature(outputTE, xBits.release(),"passbits") != HLT::OK)
      ATH_MSG_ERROR("Could not store TrigPassBits! ");

  return HLT::OK;
}
