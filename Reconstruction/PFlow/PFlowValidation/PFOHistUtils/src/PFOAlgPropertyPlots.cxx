/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "PFOHistUtils/PFOAlgPropertyPlots.h"

namespace PFO {

  PFOAlgPropertyPlots::PFOAlgPropertyPlots(PlotBase* pParent, std::string sDir, std::string sFEContainerName) : PlotBase(pParent, sDir), m_sFEContainerName(sFEContainerName){    
    m_FE_isInDenseEnvironment = nullptr;
    m_FE_tracksExpectedEnergyDeposit = nullptr;

    m_FE_isInDenseEnvironment_etaBinA = nullptr;
    m_FE_tracksExpectedEnergyDeposit_etaBinA = nullptr;

    m_FE_isInDenseEnvironment_etaBinB = nullptr;
    m_FE_tracksExpectedEnergyDeposit_etaBinB = nullptr;

    m_FE_isInDenseEnvironment_etaBinC = nullptr;
    m_FE_tracksExpectedEnergyDeposit_etaBinC = nullptr;
  }

  void PFOAlgPropertyPlots::initializePlots(){    
    // book FlowElement histograms
    if(!m_sFEContainerName.empty()){
      m_FE_isInDenseEnvironment = Book1D("_isInDenseEnvironment",m_sFEContainerName+"_isInDenseEnvironment",3,-1,2);
      m_FE_tracksExpectedEnergyDeposit = Book1D("_tracksExpectedEnergyDeposit",m_sFEContainerName+"_tracksExpectedEnergyDeposit",11,-1,10);

      m_FE_isInDenseEnvironment_etaBinA = Book1D("_isInDenseEnvironment_binA",m_sFEContainerName+"_isInDenseEnvironment (|eta| < 1)",3,-1,2);
      m_FE_tracksExpectedEnergyDeposit_etaBinA = Book1D("_tracksExpectedEnergyDeposit_binA)",m_sFEContainerName+"_tracksExpectedEnergyDeposit (|eta| < 1)",11,-1,10);
      
      m_FE_isInDenseEnvironment_etaBinB = Book1D("_isInDenseEnvironment_binB",m_sFEContainerName+"_isInDenseEnvironment (1 <= |eta| < 2)",3,-1,2);
      m_FE_tracksExpectedEnergyDeposit_etaBinB = Book1D("_tracksExpectedEnergyDeposit_binB",m_sFEContainerName+"_tracksExpectedEnergyDeposit (1 <= |eta| < 2)",11,-1,10);

      m_FE_isInDenseEnvironment_etaBinC = Book1D("_isInDenseEnvironment_binC",m_sFEContainerName+"_isInDenseEnvironment (|eta| >= 2)",3,-1,2);
      m_FE_tracksExpectedEnergyDeposit_etaBinC = Book1D("_tracksExpectedEnergyDeposit_binC",m_sFEContainerName+"_tracksExpectedEnergyDeposit (|eta| >= 2)",11,-1,10);      
    }
  }

 void PFOAlgPropertyPlots::fill(const xAOD::FlowElement& FE, const xAOD::EventInfo& eventInfo){

   static const SG::AuxElement::ConstAccessor<int> acc_IsInDenseEnvironment("IsInDenseEnvironment");
   // dump the "isInDenseEnvironment
   if(acc_IsInDenseEnvironment.isAvailable(FE)){
     int isInDenseEnvironment=acc_IsInDenseEnvironment(FE);
     m_FE_isInDenseEnvironment->Fill(isInDenseEnvironment,eventInfo.beamSpotWeight());
     if (fabs(FE.eta()) < 1) m_FE_isInDenseEnvironment_etaBinA->Fill(isInDenseEnvironment,eventInfo.beamSpotWeight());
     else if (fabs(FE.eta()) < 2) m_FE_isInDenseEnvironment_etaBinB->Fill(isInDenseEnvironment,eventInfo.beamSpotWeight());
     else m_FE_isInDenseEnvironment_etaBinC->Fill(isInDenseEnvironment,eventInfo.beamSpotWeight());
   }     
   else{ 
     m_FE_isInDenseEnvironment->Fill(-1.0,eventInfo.beamSpotWeight());
     if (fabs(FE.eta()) < 1) m_FE_isInDenseEnvironment_etaBinA->Fill(-1.0,eventInfo.beamSpotWeight());
     else if (fabs(FE.eta()) < 2) m_FE_isInDenseEnvironment_etaBinB->Fill(-1.0,eventInfo.beamSpotWeight());
     else m_FE_isInDenseEnvironment_etaBinC->Fill(-1.0,eventInfo.beamSpotWeight());     
   }
   static const SG::AuxElement::ConstAccessor<float> acc_FE_tracksExpectedEnergyDeposit("TracksExpectedEnergyDeposit");
   
   if(acc_FE_tracksExpectedEnergyDeposit.isAvailable(FE)){
     float expectedEnergy=acc_FE_tracksExpectedEnergyDeposit(FE);
     m_FE_tracksExpectedEnergyDeposit->Fill(expectedEnergy/1000.0,eventInfo.beamSpotWeight());
     if(fabs(FE.eta())<1) 
       m_FE_tracksExpectedEnergyDeposit_etaBinA->Fill(expectedEnergy/1000.0,eventInfo.beamSpotWeight());
     else if(fabs(FE.eta())<2)
       m_FE_tracksExpectedEnergyDeposit_etaBinB->Fill(expectedEnergy/1000.0,eventInfo.beamSpotWeight());
     else
       m_FE_tracksExpectedEnergyDeposit_etaBinC->Fill(expectedEnergy/1000.0,eventInfo.beamSpotWeight());
   }// end of accessor block on tracks expected energy deposit
   else{
     m_FE_tracksExpectedEnergyDeposit->Fill(-1.0,eventInfo.beamSpotWeight());
     if( fabs(FE.eta())<1) m_FE_tracksExpectedEnergyDeposit_etaBinA->Fill(-1.0,eventInfo.beamSpotWeight());
     else if ((fabs(FE.eta())<2)) m_FE_tracksExpectedEnergyDeposit_etaBinB->Fill(-1.0,eventInfo.beamSpotWeight());
     else
       m_FE_tracksExpectedEnergyDeposit_etaBinC->Fill(-1.0,eventInfo.beamSpotWeight());
   }

 } // end of PFOAlgPropertyPlots::fill(const xAOD::FlowElement& FE) 
} // end of namespace PFO
