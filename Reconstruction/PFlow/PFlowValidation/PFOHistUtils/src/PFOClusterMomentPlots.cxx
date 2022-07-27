/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "PFOHistUtils/PFOClusterMomentPlots.h"

namespace PFO {

  PFOClusterMomentPlots::PFOClusterMomentPlots(PlotBase* pParent, std::string sDir, std::string sFEContainerName) : PlotBase(pParent, sDir), m_sFEContainerName(sFEContainerName){
        
    m_FE_SECOND_R = nullptr;
    m_FE_CENTER_LAMBDA = nullptr;
    m_FE_ISOLATION = nullptr;
    m_FE_ENG_BAD_CELLS = nullptr;
    m_FE_N_BAD_CELLS = nullptr;
    m_FE_BADLARQ_FRAC = nullptr;
    m_FE_ENG_POS = nullptr;
    m_FE_AVG_LAR_Q = nullptr;
    m_FE_AVG_TILE_Q = nullptr;
    m_FE_EM_PROBABILITY = nullptr;
    m_FE_SECOND_LAMBDA = nullptr;

    m_FE_SECOND_R_etaBinA = nullptr;
    m_FE_CENTER_LAMBDA_etaBinA = nullptr;
    m_FE_ISOLATION_etaBinA = nullptr;
    m_FE_ENG_BAD_CELLS_etaBinA = nullptr;
    m_FE_N_BAD_CELLS_etaBinA = nullptr;
    m_FE_BADLARQ_FRAC_etaBinA = nullptr;
    m_FE_ENG_POS_etaBinA = nullptr;
    m_FE_AVG_LAR_Q_etaBinA = nullptr;
    m_FE_AVG_TILE_Q_etaBinA = nullptr;
    m_FE_EM_PROBABILITY_etaBinA = nullptr;
    m_FE_SECOND_LAMBDA_etaBinA = nullptr;

    m_FE_SECOND_R_etaBinB = nullptr;
    m_FE_CENTER_LAMBDA_etaBinB = nullptr;
    m_FE_ISOLATION_etaBinB = nullptr;
    m_FE_ENG_BAD_CELLS_etaBinB = nullptr;
    m_FE_N_BAD_CELLS_etaBinB = nullptr;
    m_FE_BADLARQ_FRAC_etaBinB = nullptr;
    m_FE_ENG_POS_etaBinB = nullptr;
    m_FE_AVG_LAR_Q_etaBinB = nullptr;
    m_FE_AVG_TILE_Q_etaBinB = nullptr;
    m_FE_EM_PROBABILITY_etaBinB = nullptr;
    m_FE_SECOND_LAMBDA_etaBinB = nullptr;

    m_FE_SECOND_R_etaBinC = nullptr;
    m_FE_CENTER_LAMBDA_etaBinC = nullptr;
    m_FE_ISOLATION_etaBinC = nullptr;
    m_FE_ENG_BAD_CELLS_etaBinC = nullptr;
    m_FE_N_BAD_CELLS_etaBinC = nullptr;
    m_FE_BADLARQ_FRAC_etaBinC = nullptr;
    m_FE_ENG_POS_etaBinC = nullptr;
    m_FE_AVG_LAR_Q_etaBinC = nullptr;
    m_FE_AVG_TILE_Q_etaBinC = nullptr;
    m_FE_EM_PROBABILITY_etaBinC = nullptr;
    m_FE_SECOND_LAMBDA_etaBinC = nullptr;

    m_FE_SECOND_R_etaBinD = nullptr;
    m_FE_CENTER_LAMBDA_etaBinD = nullptr;
    m_FE_ISOLATION_etaBinD = nullptr;
    m_FE_ENG_BAD_CELLS_etaBinD = nullptr;
    m_FE_N_BAD_CELLS_etaBinD = nullptr;
    m_FE_BADLARQ_FRAC_etaBinD = nullptr;
    m_FE_ENG_POS_etaBinD = nullptr;
    m_FE_AVG_LAR_Q_etaBinD = nullptr;
    m_FE_AVG_TILE_Q_etaBinD = nullptr;
    m_FE_EM_PROBABILITY_etaBinD = nullptr;
    m_FE_SECOND_LAMBDA_etaBinD = nullptr;
  }

  void PFOClusterMomentPlots::initializePlots(){
    
    // FlowElement
    if(!m_sFEContainerName.empty()){
      m_FE_SECOND_R = Book1D("_SECOND_R",m_sFEContainerName + "_SECOND_R",60,-1.0,50.0); 
      m_FE_CENTER_LAMBDA = Book1D("_CENTER_LAMBDA",m_sFEContainerName + "_CENTER_LAMBDA",60,-50.0,3000.0);
      m_FE_ISOLATION = Book1D("_ISOLATION",m_sFEContainerName + "_ISOLATION",60,-1.0,2.0);
      m_FE_ENG_BAD_CELLS = Book1D("_ENG_BAD_CELLS",m_sFEContainerName + "_ENG_BAD_CELLS",60,-1.0,5);
      m_FE_N_BAD_CELLS = Book1D("_N_BAD_CELLS",m_sFEContainerName + "_N_BAD_CELLS",30,-1.0,2.0);
      m_FE_BADLARQ_FRAC = Book1D("_BADLARQ_FRAC",m_sFEContainerName + "_BADLARQ_FRAC",25,-1.0,1.5);
      m_FE_ENG_POS = Book1D("_ENG_POS",m_sFEContainerName + "_ENG_POS",60,-100.0,10000.0);
      m_FE_AVG_LAR_Q = Book1D("_AVG_LAR_Q",m_sFEContainerName + "_AVG_LAR_Q",31,-1000.0,30000.0);
      m_FE_AVG_TILE_Q = Book1D("_AVG_TILE_Q",m_sFEContainerName + "_AVG_TILE_Q",21,-10.0,200.0);
      m_FE_EM_PROBABILITY = Book1D("_EM_PROBABILITY",m_sFEContainerName + "_EM_PROBABILITY",21,-1.0,1.0);
      m_FE_SECOND_LAMBDA = Book1D("_SECOND_LAMBDA",m_sFEContainerName + "_SECOND_LAMBDA",60,-1.0,3000.0);
      
      m_FE_SECOND_R_etaBinA = Book1D("_SECOND_R_A",m_sFEContainerName + "_SECOND_R (|eta| < 1.5)",60,-1.0,50.0); 
      m_FE_CENTER_LAMBDA_etaBinA = Book1D("_CENTER_LAMBDA_A",m_sFEContainerName + "_CENTER_LAMBDA (|eta| < 1.5)",60,-50.0,3000.0);
      m_FE_ISOLATION_etaBinA = Book1D("_ISOLATION_A",m_sFEContainerName + "_ISOLATION (|eta| < 1.5)",60,-1.0,2.0);
      m_FE_ENG_BAD_CELLS_etaBinA = Book1D("_ENG_BAD_CELLS_A",m_sFEContainerName + "_ENG_BAD_CELLS (|eta| < 1.5)",60,-1.0,5);
      m_FE_N_BAD_CELLS_etaBinA = Book1D("_N_BAD_CELLS_A",m_sFEContainerName + "_N_BAD_CELLS (|eta| < 1.5)",30,-1.0,2.0);
      m_FE_BADLARQ_FRAC_etaBinA = Book1D("_BADLARQ_FRAC_A",m_sFEContainerName + "_BADLARQ_FRAC (|eta| < 1.5)",25,-1.0,1.5);
      m_FE_ENG_POS_etaBinA = Book1D("_ENG_POS_A",m_sFEContainerName + "_ENG_POS (|eta| < 1.5)",60,-100.0,10000.0);
      m_FE_AVG_LAR_Q_etaBinA = Book1D("_AVG_LAR_Q_A",m_sFEContainerName + "_AVG_LAR_Q (|eta| < 1.5)",31,-1000.0,30000.0);
      m_FE_AVG_TILE_Q_etaBinA = Book1D("_AVG_TILE_Q_A",m_sFEContainerName + "_AVG_TILE_Q (|eta| < 1.5)",21,-10.0,200.0);
      m_FE_EM_PROBABILITY_etaBinA = Book1D("_EM_PROBABILITY_A",m_sFEContainerName + "_EM_PROBABILITY (|eta| < 1.5)",21,-1.0,1.0);
      m_FE_SECOND_LAMBDA_etaBinA = Book1D("_SECOND_LAMBDA_A",m_sFEContainerName + "_SECOND_LAMBDA (|eta| < 1.5)",60,-1.0,3000.0);
      
      m_FE_SECOND_R_etaBinB = Book1D("_SECOND_R_B",m_sFEContainerName + "_SECOND_R (1.5 <= |eta| < 2.5)",60,-1.0,50.0); 
      m_FE_CENTER_LAMBDA_etaBinB = Book1D("_CENTER_LAMBDA_B",m_sFEContainerName + "_CENTER_LAMBDA (1.5 <= |eta| < 2.5)",60,-50.0,3000.0);
      m_FE_ISOLATION_etaBinB = Book1D("_ISOLATION_B",m_sFEContainerName + "_ISOLATION (1.5 <= |eta| < 2.5)",60,-1.0,2.0);
      m_FE_ENG_BAD_CELLS_etaBinB = Book1D("_ENG_BAD_CELLS_B",m_sFEContainerName + "_ENG_BAD_CELLS (1.5 <= |eta| < 2.5)",60,-1.0,5);
      m_FE_N_BAD_CELLS_etaBinB = Book1D("_N_BAD_CELLS_B",m_sFEContainerName + "_N_BAD_CELLS (1.5 <= |eta| < 2.5)",30,-1.0,2.0);
      m_FE_BADLARQ_FRAC_etaBinB = Book1D("_BADLARQ_FRAC_B",m_sFEContainerName + "_BADLARQ_FRAC (1.5 <= |eta| < 2.5)",25,-1.0,1.5);
      m_FE_ENG_POS_etaBinB = Book1D("_ENG_POS_B",m_sFEContainerName + "_ENG_POS (1.5 <= |eta| < 2.5)",60,-100.0,10000.0);
      m_FE_AVG_LAR_Q_etaBinB = Book1D("_AVG_LAR_Q_B",m_sFEContainerName + "_AVG_LAR_Q (1.5 <= |eta| < 2.5)",31,-1000.0,30000.0);
      m_FE_AVG_TILE_Q_etaBinB = Book1D("_AVG_TILE_Q_B",m_sFEContainerName + "_AVG_TILE_Q (1.5 <= |eta| < 2.5)",21,-10.0,200.0);
      m_FE_EM_PROBABILITY_etaBinB = Book1D("_EM_PROBABILITY_B",m_sFEContainerName + "_EM_PROBABILITY (1.5 <= |eta| < 2.5)",21,-1.0,1.0);
      m_FE_SECOND_LAMBDA_etaBinB = Book1D("_SECOND_LAMBDA_B",m_sFEContainerName + "_SECOND_LAMBDA (1.5 <= |eta| < 2.5)",60,-1.0,3000.0);

      m_FE_SECOND_R_etaBinC = Book1D("_SECOND_R_C",m_sFEContainerName + "_SECOND_R (2.5 <= |eta| < 3.2)",60,-1.0,50.0); 
      m_FE_CENTER_LAMBDA_etaBinC = Book1D("_CENTER_LAMBDA_C",m_sFEContainerName + "_CENTER_LAMBDA (2.5 <= |eta| < 3.2)",60,-50.0,3000.0);
      m_FE_ISOLATION_etaBinC = Book1D("_ISOLATION_C",m_sFEContainerName + "_ISOLATION (2.5 <= |eta| < 3.2)",60,-1.0,2.0);
      m_FE_ENG_BAD_CELLS_etaBinC = Book1D("_ENG_BAD_CELLS_C",m_sFEContainerName + "_ENG_BAD_CELLS (2.5 <= |eta| < 3.2)",60,-1.0,5);
      m_FE_N_BAD_CELLS_etaBinC = Book1D("_N_BAD_CELLS_C",m_sFEContainerName + "_N_BAD_CELLS (2.5 <= |eta| < 3.2)",30,-1.0,2.0);
      m_FE_BADLARQ_FRAC_etaBinC = Book1D("_BADLARQ_FRAC_C",m_sFEContainerName + "_BADLARQ_FRAC (2.5 <= |eta| < 3.2)",25,-1.0,1.5);
      m_FE_ENG_POS_etaBinC = Book1D("_ENG_POS_C",m_sFEContainerName + "_ENG_POS (2.5 <= |eta| < 3.2)",60,-100.0,10000.0);
      m_FE_AVG_LAR_Q_etaBinC = Book1D("_AVG_LAR_Q_C",m_sFEContainerName + "_AVG_LAR_Q (2.5 <= |eta| < 3.2)",31,-1000.0,30000.0);
      m_FE_AVG_TILE_Q_etaBinC = Book1D("_AVG_TILE_Q_C",m_sFEContainerName + "_AVG_TILE_Q (2.5 <= |eta| < 3.2)",21,-10.0,200.0);
      m_FE_EM_PROBABILITY_etaBinC = Book1D("_EM_PROBABILITY_C",m_sFEContainerName + "_EM_PROBABILITY (2.5 <= |eta| < 3.2)",21,-1.0,1.0);
      m_FE_SECOND_LAMBDA_etaBinC = Book1D("_SECOND_LAMBDA_C",m_sFEContainerName + "_SECOND_LAMBDA (2.5 <= |eta| < 3.2)",60,-1.0,3000.0);
      
      m_FE_SECOND_R_etaBinD = Book1D("_SECOND_R_D",m_sFEContainerName + "_SECOND_R (|eta| >= 3.2)",60,-1.0,50.0); 
      m_FE_CENTER_LAMBDA_etaBinD = Book1D("_CENTER_LAMBDA_D",m_sFEContainerName + "_CENTER_LAMBDA (|eta| >= 3.2)",60,-50.0,3000.0);
      m_FE_ISOLATION_etaBinD = Book1D("_ISOLATION_D",m_sFEContainerName + "_ISOLATION (|eta| >= 3.2)",60,-1.0,2.0);
      m_FE_ENG_BAD_CELLS_etaBinD = Book1D("_ENG_BAD_CELLS_D",m_sFEContainerName + "_ENG_BAD_CELLS (|eta| >= 3.2)",60,-1.0,5);
      m_FE_N_BAD_CELLS_etaBinD = Book1D("_N_BAD_CELLS_D",m_sFEContainerName + "_N_BAD_CELLS (|eta| >= 3.2)",30,-1.0,2.0);
      m_FE_BADLARQ_FRAC_etaBinD = Book1D("_BADLARQ_FRAC_D",m_sFEContainerName + "_BADLARQ_FRAC (|eta| >= 3.2)",25,-1.0,1.5);
      m_FE_ENG_POS_etaBinD = Book1D("_ENG_POS_D",m_sFEContainerName + "_ENG_POS (|eta| >= 3.2)",60,-100.0,10000.0);
      m_FE_AVG_LAR_Q_etaBinD = Book1D("_AVG_LAR_Q_D",m_sFEContainerName + "_AVG_LAR_Q (|eta| >= 3.2)",31,-1000.0,30000.0);
      m_FE_AVG_TILE_Q_etaBinD = Book1D("_AVG_TILE_Q_D",m_sFEContainerName + "_AVG_TILE_Q (|eta| >= 3.2)",21,-10.0,200.0);
      m_FE_EM_PROBABILITY_etaBinD = Book1D("_EM_PROBABILITY_D",m_sFEContainerName + "_EM_PROBABILITY (|eta| >= 3.2)",21,-1.0,1.0);
      m_FE_SECOND_LAMBDA_etaBinD = Book1D("_SECOND_LAMBDA_D",m_sFEContainerName + "_SECOND_LAMBDA (|eta| >= 3.2)",60,-1.0,3000.0);
    } 
  }

  
  void PFOClusterMomentPlots::fill(const xAOD::FlowElement& FE, const xAOD::EventInfo& eventInfo){
     float moment_SECOND_R = -1.0;
     float moment_CENTER_LAMBDA = -1.0;
     float moment_ISOLATION = -1.0;
     float moment_ENG_BAD_CELLS=-1.0;
     float moment_N_BAD_CELLS = -1.0;     
     float moment_BADLARQ_FRAC = -1.0;
     float moment_ENG_POS = -1.0;
     float moment_AVG_LAR_Q = -1.0;
     float moment_AVG_TILE_Q = -1.0; 
     float moment_EM_PROBABILITY = -1.0;
     float moment_SECOND_LAMBDA = -1.0;

     //as opposed to PFO which uses specific functions to grab the cluster moments, the auxdata is used
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_SECOND_R("SECOND_R");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_CENTER_LAMBDA("CENTER_LAMBDA");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_ISOLATION("ISOLATION");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_ENG_BAD_CELLS("ENG_BAD_CELLS");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_N_BAD_CELLS("N_BAD_CELLS");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_BADLARQ_FRAC("BADLARQ_FRAC");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_ENG_POS("ENG_POS");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_AVG_LAR_Q("AVG_LAR_Q");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_AVG_TILE_Q("AVG_TILE_Q");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_EM_PROBABILITY("EM_PROBABILITY");
     static const SG::AuxElement::ConstAccessor<float>acc_FE_moment_SECOND_LAMBDA("SECOND_LAMBDA");
     
     //use accessors to retrieve the auxvars
     if(acc_FE_moment_SECOND_R.isAvailable(FE))
       moment_SECOND_R=acc_FE_moment_SECOND_R(FE);
     
     if(acc_FE_moment_CENTER_LAMBDA.isAvailable(FE))
       moment_CENTER_LAMBDA=acc_FE_moment_CENTER_LAMBDA(FE);
     
     if(acc_FE_moment_ISOLATION.isAvailable(FE))
       moment_ISOLATION=acc_FE_moment_ISOLATION(FE);
     
     if(acc_FE_moment_ENG_BAD_CELLS.isAvailable(FE))
       moment_ENG_BAD_CELLS=acc_FE_moment_ENG_BAD_CELLS(FE);
     
     if(acc_FE_moment_N_BAD_CELLS.isAvailable(FE))
       moment_N_BAD_CELLS=acc_FE_moment_N_BAD_CELLS(FE);
     
     if(acc_FE_moment_BADLARQ_FRAC.isAvailable(FE))
       moment_BADLARQ_FRAC=acc_FE_moment_BADLARQ_FRAC(FE);
     
     if(acc_FE_moment_ENG_POS.isAvailable(FE))
       moment_ENG_POS=acc_FE_moment_ENG_POS(FE);
     
     if(acc_FE_moment_AVG_LAR_Q.isAvailable(FE))
       moment_AVG_LAR_Q=acc_FE_moment_AVG_LAR_Q(FE);
     
     if(acc_FE_moment_AVG_TILE_Q.isAvailable(FE))
       moment_AVG_TILE_Q=acc_FE_moment_AVG_TILE_Q(FE);
     
     if(acc_FE_moment_EM_PROBABILITY.isAvailable(FE))
       moment_EM_PROBABILITY=acc_FE_moment_EM_PROBABILITY(FE);
     
     if(acc_FE_moment_SECOND_LAMBDA.isAvailable(FE))
       moment_SECOND_LAMBDA=acc_FE_moment_SECOND_LAMBDA(FE);
     
     float FE_eta=FE.eta();
     m_FE_SECOND_R->Fill(moment_SECOND_R,eventInfo.beamSpotWeight());     
     m_FE_CENTER_LAMBDA->Fill(moment_CENTER_LAMBDA,eventInfo.beamSpotWeight());
     m_FE_ISOLATION->Fill(moment_ISOLATION,eventInfo.beamSpotWeight());     
     m_FE_ENG_BAD_CELLS->Fill(moment_ENG_BAD_CELLS,eventInfo.beamSpotWeight());
     m_FE_N_BAD_CELLS->Fill(moment_N_BAD_CELLS,eventInfo.beamSpotWeight());
     m_FE_BADLARQ_FRAC->Fill(moment_BADLARQ_FRAC,eventInfo.beamSpotWeight());
     m_FE_ENG_POS->Fill(moment_ENG_POS,eventInfo.beamSpotWeight());
     m_FE_AVG_LAR_Q->Fill(moment_AVG_LAR_Q,eventInfo.beamSpotWeight());
     m_FE_AVG_TILE_Q->Fill(moment_AVG_TILE_Q,eventInfo.beamSpotWeight());
     m_FE_EM_PROBABILITY->Fill(moment_EM_PROBABILITY,eventInfo.beamSpotWeight());
     m_FE_SECOND_LAMBDA->Fill(moment_SECOND_LAMBDA,eventInfo.beamSpotWeight());

     if (fabs(FE_eta) < 1.5){
       m_FE_SECOND_R_etaBinA->Fill(moment_SECOND_R,eventInfo.beamSpotWeight());
       m_FE_CENTER_LAMBDA_etaBinA->Fill(moment_CENTER_LAMBDA,eventInfo.beamSpotWeight());
       m_FE_ISOLATION_etaBinA->Fill(moment_ISOLATION,eventInfo.beamSpotWeight());
       m_FE_ENG_BAD_CELLS_etaBinA->Fill(moment_ENG_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_N_BAD_CELLS_etaBinA->Fill(moment_N_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_BADLARQ_FRAC_etaBinA->Fill(moment_BADLARQ_FRAC,eventInfo.beamSpotWeight());
       m_FE_ENG_POS_etaBinA->Fill(moment_ENG_POS,eventInfo.beamSpotWeight());
       m_FE_AVG_LAR_Q_etaBinA->Fill(moment_AVG_LAR_Q,eventInfo.beamSpotWeight());
       m_FE_AVG_TILE_Q_etaBinA->Fill(moment_AVG_TILE_Q,eventInfo.beamSpotWeight());
       m_FE_EM_PROBABILITY_etaBinA->Fill(moment_EM_PROBABILITY,eventInfo.beamSpotWeight());
       m_FE_SECOND_LAMBDA_etaBinA->Fill(moment_SECOND_LAMBDA,eventInfo.beamSpotWeight());
     }//|eta| < 1.5
     else if (fabs(FE_eta) < 2.5){
       m_FE_SECOND_R_etaBinB->Fill(moment_SECOND_R,eventInfo.beamSpotWeight());
       m_FE_CENTER_LAMBDA_etaBinB->Fill(moment_CENTER_LAMBDA,eventInfo.beamSpotWeight());
       m_FE_ISOLATION_etaBinB->Fill(moment_ISOLATION,eventInfo.beamSpotWeight());
       m_FE_ENG_BAD_CELLS_etaBinB->Fill(moment_ENG_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_N_BAD_CELLS_etaBinB->Fill(moment_N_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_BADLARQ_FRAC_etaBinB->Fill(moment_BADLARQ_FRAC,eventInfo.beamSpotWeight());
       m_FE_ENG_POS_etaBinB->Fill(moment_ENG_POS,eventInfo.beamSpotWeight());
       m_FE_AVG_LAR_Q_etaBinB->Fill(moment_AVG_LAR_Q,eventInfo.beamSpotWeight());
       m_FE_AVG_TILE_Q_etaBinB->Fill(moment_AVG_TILE_Q,eventInfo.beamSpotWeight());
       m_FE_EM_PROBABILITY_etaBinB->Fill(moment_EM_PROBABILITY,eventInfo.beamSpotWeight());
       m_FE_SECOND_LAMBDA_etaBinB->Fill(moment_SECOND_LAMBDA,eventInfo.beamSpotWeight());
     }
     else if (fabs(FE_eta) < 3.2){
       m_FE_SECOND_R_etaBinC->Fill(moment_SECOND_R,eventInfo.beamSpotWeight());
       m_FE_CENTER_LAMBDA_etaBinC->Fill(moment_CENTER_LAMBDA,eventInfo.beamSpotWeight());
       m_FE_ISOLATION_etaBinC->Fill(moment_ISOLATION,eventInfo.beamSpotWeight());
       m_FE_ENG_BAD_CELLS_etaBinC->Fill(moment_ENG_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_N_BAD_CELLS_etaBinC->Fill(moment_N_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_BADLARQ_FRAC_etaBinC->Fill(moment_BADLARQ_FRAC,eventInfo.beamSpotWeight());
       m_FE_ENG_POS_etaBinC->Fill(moment_ENG_POS,eventInfo.beamSpotWeight());
       m_FE_AVG_LAR_Q_etaBinC->Fill(moment_AVG_LAR_Q,eventInfo.beamSpotWeight());
       m_FE_AVG_TILE_Q_etaBinC->Fill(moment_AVG_TILE_Q,eventInfo.beamSpotWeight());
       m_FE_EM_PROBABILITY_etaBinC->Fill(moment_EM_PROBABILITY,eventInfo.beamSpotWeight());
       m_FE_SECOND_LAMBDA_etaBinC->Fill(moment_SECOND_LAMBDA,eventInfo.beamSpotWeight());
     }
     else{
       m_FE_SECOND_R_etaBinD->Fill(moment_SECOND_R,eventInfo.beamSpotWeight());
       m_FE_CENTER_LAMBDA_etaBinD->Fill(moment_CENTER_LAMBDA,eventInfo.beamSpotWeight());
       m_FE_ISOLATION_etaBinD->Fill(moment_ISOLATION,eventInfo.beamSpotWeight());
       m_FE_ENG_BAD_CELLS_etaBinD->Fill(moment_ENG_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_N_BAD_CELLS_etaBinD->Fill(moment_N_BAD_CELLS,eventInfo.beamSpotWeight());
       m_FE_BADLARQ_FRAC_etaBinD->Fill(moment_BADLARQ_FRAC,eventInfo.beamSpotWeight());
       m_FE_ENG_POS_etaBinD->Fill(moment_ENG_POS,eventInfo.beamSpotWeight());
       m_FE_AVG_LAR_Q_etaBinD->Fill(moment_AVG_LAR_Q,eventInfo.beamSpotWeight());
       m_FE_AVG_TILE_Q_etaBinD->Fill(moment_AVG_TILE_Q,eventInfo.beamSpotWeight());
       m_FE_EM_PROBABILITY_etaBinD->Fill(moment_EM_PROBABILITY,eventInfo.beamSpotWeight());
       m_FE_SECOND_LAMBDA_etaBinD->Fill(moment_SECOND_LAMBDA,eventInfo.beamSpotWeight());
     }     
     
     
  }
}
