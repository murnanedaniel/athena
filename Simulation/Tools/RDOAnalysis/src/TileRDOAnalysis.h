/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#ifndef TILE_RDO_ANALYSIS_H
#define TILE_RDO_ANALYSIS_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"

#include "TileEvent/TileRawChannelContainer.h"
#include "TileEvent/TileRawChannelCollection.h"
#include "TileEvent/TileContainer.h"
#include "TileEvent/TileDigitsContainer.h"
#include "TileEvent/TileDigitsCollection.h"

#include <string>
#include <vector>
#include "TH1.h"

class TTree;
class TH1;

class TileRDOAnalysis : public AthAlgorithm {

 public:
  TileRDOAnalysis(const std::string& name, ISvcLocator* pSvcLocator);
  ~TileRDOAnalysis(){}

  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();

 private:
  // TileRawChannel
  // AMP, TIME, QUAL REALLY VECTORS - CHECK SIZE/OUTPUT
  std::vector<unsigned long long>* m_adcID;
  std::vector<unsigned long long>* m_pmtID;
  std::vector<unsigned long long>* m_cellID;
  std::vector<unsigned long long>* m_ttID;
  std::vector<unsigned long long>* m_mtID;
  std::vector<int>* m_fragID;
  std::vector<float>* m_rawAmp;
  std::vector<float>* m_rawTime;
  std::vector<float>* m_rawQual;
  std::vector<float>* m_rawPed;
  std::vector<unsigned long long>* m_adcID_mu;
  std::vector<unsigned long long>* m_pmtID_mu;
  std::vector<unsigned long long>* m_cellID_mu;
  std::vector<unsigned long long>* m_ttID_mu;
  std::vector<unsigned long long>* m_mtID_mu;
  std::vector<int>* m_fragID_mu;
  std::vector<float>* m_rawAmp_mu;
  std::vector<float>* m_rawTime_mu;
  std::vector<float>* m_rawQual_mu;
  std::vector<float>* m_rawPed_mu;
  // TileMuonReceiver
  std::vector<int>* m_muRcvID;
  std::vector<bool>* m_muRcv_dec;
  std::vector<float>* m_muRcv_thresh;
  std::vector<float>* m_muRcv_energy;
  std::vector<float>* m_muRcv_time;
  // TileTTL1
  std::vector<unsigned long long>* m_ttl1MBTS_ID;
  std::vector< std::vector<double> >* m_ttl1MBTS_digits;
  std::vector<unsigned long long>* m_ttl1_ID;
  std::vector< std::vector<double> >* m_ttl1_digits;
  // TileL2
  std::vector<int>* m_L2ID;
  std::vector< std::vector<unsigned int> >* m_L2val;
  std::vector< std::vector<float> >* m_L2eta;
  std::vector<float>* m_L2phi;
  std::vector< std::vector<float> >* m_L2energyA;
  std::vector< std::vector<float> >* m_L2energyBC;
  std::vector< std::vector<float> >* m_L2energyD;
  std::vector< std::vector<unsigned int> >* m_L2qual;
  std::vector< std::vector<float> >* m_L2sumE;
  // Tile Digits
  std::vector<uint32_t>* m_fragSize;
  std::vector<uint32_t>* m_fragBCID;
  std::vector< std::vector<double> >* m_digits;
  std::vector<uint32_t>* m_muFragSize;
  std::vector<uint32_t>* m_muFragBCID;
  std::vector< std::vector<double> >* m_muDigits;

  // HISTOGRAMS
  TH1* h_adcID;
  TH1* h_rawAmp;
  TH1* h_rawTime;
  TH1* h_rawQual;
  TH1* h_rawPed;
  TH1* h_adcID_mu;
  TH1* h_rawAmp_mu;
  TH1* h_rawTime_mu;
  TH1* h_rawQual_mu;
  TH1* h_rawPed_mu;
  TH1* h_muRcvID;
  TH1* h_muRcv_dec;
  TH1* h_muRcv_thresh;
  TH1* h_muRcv_energy;
  TH1* h_muRcv_time;
  TH1* h_ttl1MBTS_ID;
  TH1* h_ttl1MBTS_digits;
  TH1* h_ttl1_ID;
  TH1* h_ttl1_digits;
  TH1* h_L2ID;
  TH1* h_L2val;
  TH1* h_L2eta;
  TH1* h_L2phi;
  TH1* h_L2energyA;
  TH1* h_L2energyBC;
  TH1* h_L2energyD;
  TH1* h_L2qual;
  TH1* h_L2sumE;
  TH1* h_digits;
  TH1* h_muDigits;

  
  TTree *m_tree;
  std::string m_ntupleFileName;
  std::string m_ntupleDirName;
  std::string m_ntupleTreeName;
  std::string m_path;
  ServiceHandle<ITHistSvc> m_thistSvc;

};

#endif // TILE_RDO_ANALYSIS_H
