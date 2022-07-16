/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// TFCGEtaSlice.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include "TFCGEtaSlice.h"

#include "CLHEP/Random/RandFlat.h"

#include "TFitResult.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"

TFCGEtaSlice::TFCGEtaSlice(int pid, int etaMin, int etaMax) {
  m_pid = pid;
  m_etaMin = etaMin;
  m_etaMax = etaMax;
}

void TFCGEtaSlice::SetGANVersion(int version){
  m_ganVersion = version;
}

void TFCGEtaSlice::SetLatentSpaceSize(int latentDim){
  m_latentDim = latentDim;
}

bool TFCGEtaSlice::IsGanCorrectlyLoaded(){
  if (m_ganVersion == 1){
    if(m_graph==nullptr){
      return false;
    } 
  }
  else{
    if(m_graph_high==nullptr || m_graph_low==nullptr){
      return false;
    }     
  }
  return true;
}

TFCGEtaSlice::FitResultsPerLayer TFCGEtaSlice::GetFitResults(){
  return m_allFitResults;
}

bool TFCGEtaSlice::LoadGAN(std::string inputFolderName, std::vector<int> relevantLayers){
  m_inputFolderName = inputFolderName;
  
  std::string inputFileName;
  if (m_ganVersion == 1 || m_pid == 211){
    inputFileName = m_inputFolderName + "/neural_net_" + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_" + std::to_string(m_etaMax) +".json";
    return LoadGANFromFileV1(inputFileName, relevantLayers);    
  }
  else{
    bool returnValue;
    inputFileName = m_inputFolderName + "/neural_net_" + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_" + std::to_string(m_etaMax) +"_High12.json";
    returnValue = LoadGANFromFileV2(inputFileName, "high", relevantLayers);    
    if (!returnValue) {
      std::cout<<"ERROR loading GAN"<< std::endl;
      return returnValue;
    }

    inputFileName = m_inputFolderName + "/neural_net_" + std::to_string(m_pid) + "_eta_" + std::to_string(m_etaMin) + "_" + std::to_string(m_etaMax) + "_UltraLow12.json";
    return LoadGANFromFileV2(inputFileName, "low", relevantLayers);    
  }
}

bool TFCGEtaSlice::LoadGANFromFileV1(std::string inputFileName, std::vector<int> relevantLayers){
  std::string inputFile=PathResolverFindCalibFile(inputFileName);
  if (inputFile.empty()){
    //std::cout<<"Could not find json file " << inputFile << std::endl;
    return false;
  } 
  else {
    //std::cout<<"For pid: " << m_pid <<" and eta " << m_etaMin <<"-" << m_etaMax <<", loading json file " << inputFile << std::endl;
    std::ifstream input(inputFile);
    // build the graph
    m_graph=std::make_unique<lwt::LightweightGraph>(lwt::parse_json_graph(input) );
    if (m_graph==nullptr){
      //std::cout<<"Could not create LightWeightGraph from  " << inputFile << std::endl;
      return false;
    }

    if (m_region->GetGANVersion() > 1){
      //std::cout<<"Pions in v2"<<std::endl;
      CalculateMeanPointFromDistributionOfR(relevantLayers);
    }
  }
  
  return true;
}

bool TFCGEtaSlice::LoadGANFromFileV2(std::string inputFileName, std::string energyRange, std::vector<int> relevantLayers){
  std::string inputFile=PathResolverFindCalibFile(inputFileName);
  if (inputFile.empty()){
    std::cout<<"Could not find json file " << inputFile << std::endl;
    return false;
  } 
  else {
    //std::cout<<"For pid: " << m_pid <<" and eta " << m_etaMin <<"-" << m_etaMax <<", loading json file " << inputFile << std::endl;
    std::ifstream input(inputFile);
    // build the graph
    if (energyRange == "high"){
      m_graph_high=std::make_unique<lwt::LightweightGraph>(lwt::parse_json_graph(input) );
      if (m_graph_high==nullptr){
        //std::cout<<"Could not create LightWeightGraph from  " << inputFile << std::endl;
        return false;
      }
    }
    else if (energyRange == "low"){
      m_graph_low=std::make_unique<lwt::LightweightGraph>(lwt::parse_json_graph(input) );
      if (m_graph_low==nullptr){
        std::cout<<"Could not create LightWeightGraph from  " << inputFile << std::endl;
        return false;
      }
    }
    else{
      return false;
    }
    
    CalculateMeanPointFromDistributionOfR();
  }
  
  return true;
}

void TFCGEtaSlice::CalculateMeanPointFromDistributionOfR(, std::vector<int> relevantLayers){
  std::string rootFileName = m_inputFolderName + "/../../../rootFiles/pid"+std::to_string(m_pid)+"_E1048576_eta_"+std::to_string(m_etaMin)+"_"+std::to_string(m_etaMin+5)+".root";
  std::cout<<"Opening file "<<rootFileName<<std::endl;
  TFile* file = TFile::Open(rootFileName.c_str(), "read");
  for (int layer : relevantLayers){
    std::string histoName="r"+std::to_string(layer)+"w";
    TH1D* h1 = (TH1D*)file->Get(histoName.c_str());
    if (TMath::IsNaN(h1->Integral())){
        histoName="r"+std::to_string(layer);
        h1 = (TH1D*)file->Get(histoName.c_str());
    }
    
    TAxis* x = (TAxis*)h2->GetXaxis();
    for (int ix = 1; ix <= xBinNum; ++ix){
      h1->GetXaxis()->SetRangeUser(x->GetBinLowEdge(ix), x->GetBinUpEdge(ix));
      //ATH_MSG_VERBOSE("Layer "<<layer<<" Bin "<<ix -1);
      //ATH_MSG_VERBOSE(x->GetBinLowEdge(ix)<<" "<< x->GetBinUpEdge(ix)<< " " <<h1->Integral());
      if (h1->Integral() > 0){
        TFitResultPtr res(0); 
        FitResult result;

        res = h1->Fit("expo","SQ");
        if (res >=0 && !isnan(res->Parameter(0))){
          result.constant = res->Parameter(0);
          result.exp = res->Parameter(1);   
        }
        else{
          result.constant = 1;
          result.exp = 0;   
        }
        //ATH_MSG_VERBOSE(result.constant<<" "<<result.exp);
        m_allFitResults[layer].push_back(result);
      }
    }
  }
}

TFCGEtaSlice::NetworkOutputs TFCGEtaSlice::GetNetworkOutputs(const TFCSTruthState* truth, TFCSExtrapolationState extrapol, TFCSSimulationState simulstate) {
  double randUniformZ = 0.;
  NetworkInputs inputs;

  int maxExp;
  if (truth->P() > 4096){ //This is the momentum, not the energy, because the split is based on the samples which are produced with the momentum
    maxExp = 22;
  }
  else{
    maxExp = 12;
  }

  for (int i = 0; i< m_latentDim; i ++)
  {
    randUniformZ = CLHEP::RandFlat::shoot(simulstate.randomEngine(), -1., 1.);
    inputs["node_0"].insert ( std::pair<std::string,double>(std::to_string(i), randUniformZ) );
  }
  
  //std::cout << "Check label: " <<trueEnergy <<" "<< std::pow(2,22)<<" "<<trueEnergy/std::pow(2,22)<< std::endl;
  inputs["node_1"].insert ( std::pair<std::string,double>("0", truth->Ekin()/(std::pow(2,maxExp))) );

  if(m_ganVersion >= 2){
    if (false){ //conditioning on eta, should only be needed in transition regions and added only to the GANs that use it, for now all GANs have 3 conditioning inputs so filling zeros
      inputs["node_1"].insert ( std::pair<std::string,double>("2", fabs(extrapol.CaloSurface_eta())) );
    }
    else {
      inputs["node_1"].insert ( std::pair<std::string,double>("2", 0) );
    }
  }
  
  if (m_ganVersion == 1 || m_pid == 211){
    return m_graph->compute(inputs);
  }
  else{
    if (isfp.momentum().mag() > 4096){ //This is the momentum, not the energy, because the split is based on the samples which are produced with the momentum
      return m_graph_high->compute(inputs);
    }
    else{
      return m_graph_low->compute(inputs);
    }
  }    
}
