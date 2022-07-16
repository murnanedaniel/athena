/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
 
//////////////////////////////////////////////////////////////////
// TFCGEtaSlice.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef ISF_TFCGETASLICE_H
#define ISF_TFCGETASLICE_H 1

// ISF includes
#include <map>
#include <vector>

#include "ISF_FastCaloSimEvent/TFCSTruthState.h"
#include "ISF_FastCaloSimEvent/TFCSSimulationState.h"
#include "ISF_FastCaloSimEvent/TFCSExtrapolationState.h"

#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/parse_json.hh"

#include <fstream>

class TFCGEtaSlice {
  public:
  TFCGEtaSlice(int pid, int etaMin, int etaMax);
    
  struct FitResult{
    double constant;
    double exp;
  };

  typedef std::map<int, std::vector<FitResult>> FitResultsPerLayer ;
  typedef std::map<std::string, std::map<std::string, double> >  NetworkInputs ;
  typedef std::map<std::string, double> NetworkOutputs;

  bool LoadGAN(std::string inputFolderName, bool runSingleGAN);
  void CalculateMeanPointFromDistributionOfR();
  NetworkOutputs GetNetworkOutputs(const TFCSTruthState* truth, TFCSExtrapolationState extrapol, TFCSSimulationState simulstate);
  
  bool IsGanCorrectlyLoaded();
  FitResultsPerLayer GetFitResults();
  void SetLatentSpaceSize(int latentDim);
  void SetGANVersion(int version);
  int GetGANVersion();

  private:
  int m_pid;
  int m_etaMin;
  int m_etaMax;
  int m_ganVersion;
  int m_latentDim;
  
  std::string m_inputFolderName;
  
  std::unique_ptr<lwt::LightweightGraph> m_graph;
  std::unique_ptr<lwt::LightweightGraph> m_graph_high;
  std::unique_ptr<lwt::LightweightGraph> m_graph_low;
  FitResultsPerLayer m_allFitResults;

  bool LoadGANFromFileV1(std::string inputFileName);
  bool LoadGANFromFileV2(std::string inputFileName, std::string energyRange);
};

#endif //> !ISF_TFCGETASLICE_H
