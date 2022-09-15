/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// **********************************************************************
// $Id: MDTMLOverview.cxx,v 1.0 2008/10/08 Valerio Consorti
// **********************************************************************

#include "dqm_algorithms/MDTMLOverview.h"



#include <TClass.h>
#include <TH1.h>
#include <TAxis.h>
#include <TF1.h>
#include <TProfile.h>

#include "dqm_core/exceptions.h"
#include "dqm_core/AlgorithmConfig.h"
#include "dqm_core/AlgorithmManager.h"
#include "dqm_core/Result.h"
#include "dqm_algorithms/tools/AlgorithmHelper.h"
#include "ers/ers.h"

#include <cmath>
#include <iostream>
#include <map>
#include <list>
#include <string>
#include <vector>

static dqm_algorithms::MDTMLOverview staticInstance;


namespace dqm_algorithms {

// *********************************************************************
// Public Methods
// *********************************************************************

MDTMLOverview::MDTMLOverview()
  : m_name("MDTMLoverview")
{
  dqm_core::AlgorithmManager::instance().registerAlgorithm( m_name, this );
}


MDTMLOverview::~MDTMLOverview()
{
}


dqm_core::Algorithm*
MDTMLOverview::clone()
{
  return new MDTMLOverview(*this);
}


dqm_core::Result*
MDTMLOverview::execute( const std::string& name, const TObject& object, const dqm_core::AlgorithmConfig& config)
{
  const TH1 * hist;
  const TH1 * ref;
  
  if( object.IsA()->InheritsFrom( "TH1" ) ) {
    hist = static_cast<const TH1*>(&object);
    if (hist->GetDimension() >= 2 ){
      throw dqm_core::BadConfig( ERS_HERE, name, "dimension >= 2 " );
    }
  } else {
    throw dqm_core::BadConfig( ERS_HERE, name, "does not inherit from TH1" );
  }   
  
  //Get Parameters and Thresholds
  bool ref_y_n=0;
  double thresh; 
  double minstat; 
  double greenTh;
  double redTh;
  try {
    ref_y_n      = dqm_algorithms::tools::GetFirstFromMap("ref", config.getParameters(), 0);
    thresh       = dqm_algorithms::tools::GetFirstFromMap("thresh", config.getParameters());
    minstat      = dqm_algorithms::tools::GetFirstFromMap("MinStat", config.getParameters(), 50);
    redTh        = dqm_algorithms::tools::GetFromMap( "Limits", config.getRedThresholds());
    greenTh      = dqm_algorithms::tools::GetFromMap( "Limits", config.getGreenThresholds() );
  }
  catch ( dqm_core::Exception & ex ) {
    throw dqm_core::BadConfig( ERS_HERE, name, ex.what(), ex );
  }

  //Get Reference Histo
 
  try {
    ref = static_cast<const TH1*>( config.getReference() );
  }
  catch ( dqm_core::Exception & ex ) {
    throw dqm_core::BadRefHist(ERS_HERE,name," Could not retrieve reference");
  }

  if (hist->GetDimension() != ref->GetDimension() ) {
    throw dqm_core::BadRefHist( ERS_HERE, name, "Reference VS histo: Different dimension!" );
  } 
  if (hist->GetNbinsX() != ref->GetNbinsX() ) {
    throw dqm_core::BadRefHist( ERS_HERE, name, "Reference VS histo: Different bin number in X axis!" );
  } 
    
 //Check of statistics
  if (hist->GetEntries() < minstat ) {
    ERS_INFO("Histogram does not satisfy MinStat requirement " <<hist->GetName());
    dqm_core::Result *result = new dqm_core::Result(dqm_core::Result::Undefined);
    result->tags_["InsufficientEntries"] = hist->GetEntries();
    return result;
  }
  ERS_DEBUG(1,"Statistics: "<< hist->GetEntries()<< " entries ");

 
 //Algo

  int binX = hist->GetNbinsX();
  //double x_center=0;

  std::vector<int> hist_buffer;
  std::vector<int> ref_buffer;
  
  std::vector< double > new_empty_bins;
  int count=0;

  for(int x_index=1; x_index<=binX; x_index++){
    if( hist->GetBinContent(x_index) != 0 ) hist_buffer.push_back((int)hist->GetBinContent(x_index));
    if(ref_y_n==1) {
      if( ref->GetBinContent(x_index) != 0 ) ref_buffer.push_back((int)ref->GetBinContent(x_index));
    }
  };
  
  
 
  auto median=[](std::vector<int> & v)->int {
    const auto midPoint = v.begin()+v.size()/2;
    std::nth_element(v.begin(), midPoint, v.end());
    return v[v.size()/2];
  };
  
  
  double percentThreshold = thresh*0.01;
  const int mediana_hist = median(hist_buffer);
  
  if(ref_y_n==1) {
    const int mediana_ref = median(ref_buffer);
    for(int xi=1; xi<=binX; xi++){
      if( hist->GetBinContent(xi) < mediana_hist*percentThreshold && ref->GetBinContent(xi) >= mediana_ref*percentThreshold){
        count++;
        new_empty_bins.push_back(hist->GetBinCenter(xi));
      };
    };
  } else if(ref_y_n==0){ 
    for(int xi=1; xi<=binX; xi++){
      if( hist->GetBinContent(xi) < mediana_hist*percentThreshold){
        count++;
        new_empty_bins.push_back(hist->GetBinCenter(xi));
      };
    };
  };

  
  //double tot=binX;
  //if(ref_y_n==0) tot=binX;

 //Tag generator

  dqm_core::Result* result = new dqm_core::Result(); 

  std::string message;
  std::string eta="ETA_";
  std::string ml="__ML";
  char eta_num[3];
  std::string ml_num;
  double num_ml_off=1;
  int counter=0;
  
  int empty_bin_number=new_empty_bins.size();

  for(int i=0; i<empty_bin_number; i++){
    counter=-1;
    do{
      counter++;
      if(counter>=20) break;
    }while(std::abs(counter-new_empty_bins[i])>0.3);
      
    snprintf(eta_num,sizeof(eta_num),"%d",counter);
    if( (counter-new_empty_bins[i])>=0)ml_num="1";
    if( (counter-new_empty_bins[i])<0) ml_num="2";
    message=eta+(std::string)eta_num+ml+ml_num;
    result->tags_[message]=num_ml_off;
  };

 //Result

  result->tags_["00-number_of_off_ML"] = count;

  if (count<greenTh) {
    result->status_ = dqm_core::Result::Green; 
    ERS_DEBUG(1,"Green");
  } else if ( count>=greenTh && count<redTh ) {
    result->status_ = dqm_core::Result::Yellow;
    ERS_DEBUG(1,"Yellow");
  } else if (count>=redTh) {
    result->status_ = dqm_core::Result::Red;
    ERS_DEBUG(1,"Red");
  }
  return result;
}


void 
MDTMLOverview::printDescription(std::ostream& out)
{
  std::string message;
  message += "\n";
  message += "Algorithm: \"" + m_name + "\"\n";
  message += "Description: compare the histo with reference and check if there are any ML switched off\n";
  message += "Mandatory Green/Red Threshold: Limits: number of new ML switched off while are on in reference\n";
  message += "                     thresh  = % of the mean entries per ML under which the algo will set ML off\n";
  message += "Optional Parameters: MinStat = Minimum histogram statistics needed to perform Algorithm\n";
  message += "			   ref     = Bool value: 1 if you want perform the check comparing with ref\n";
  message += "			                         0 if you want perform the check without the comparison with ref\n";
  message += "\n";
  
  out << message;
}

} // namespace dqm_algorithms
