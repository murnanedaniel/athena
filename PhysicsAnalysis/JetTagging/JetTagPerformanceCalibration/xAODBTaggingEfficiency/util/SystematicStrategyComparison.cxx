/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


//#ifdef XAOD_STANDALONE
#include "xAODRootAccess/TStore.h"
//#endif // XAOD_STANDALONE

#include <AsgTools/StandaloneToolHandle.h>
#include "FTagAnalysisInterfaces/IBTaggingEfficiencyTool.h"

#include <string>
#include <iomanip>

using CP::CorrectionCode;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// This executable cycles through the different calibrated flavour taggers                             //
//// in a CDI file, and tests out both the SFEigen and SFGlobalEigen systematic strategies              //
//// to compare their output and performance.                                                         //
//// Use this in unison with the 'validate_reduction' function in CalibrationDataEigenVariations.cxx//
//// to plot the bin-to-bin correlations - both for the "true" covariance matrix, and the         //
//// approximate covariance matrix constructed **after** the PCA reduction takes place          //
////////////////////////////////////////////////////////////////////////////////////////////////



int main() {
  bool retval = true;
  xAOD::TStore store;
  

  //std::string taggerName = "DL1";
  //std::string taggerName = "MV2c10";
  //std::string workingPointName = "FixedCutBEff_77";
  //std::string jetauthorName = "AntiKt4EMPFlowJets_BTagging201903";

  std::vector<std::string> taggers;
  taggers.push_back("MV2c10");
  taggers.push_back("DL1");
  taggers.push_back("DL1r");
  //taggers.push_back("DL1rmu");
  //taggers.push_back("DL1dv00");
  
  std::vector<std::string> wps;
  wps.push_back("FixedCutBEff_60");
  wps.push_back("FixedCutBEff_70");
  wps.push_back("FixedCutBEff_77");
  wps.push_back("FixedCutBEff_85");
  wps.push_back("Continuous");

  std::vector<std::string> strats;
  strats.push_back("SFGlobalEigen");
  strats.push_back("SFEigen");
  
  //std::vector<std::string> jets;
  //jets.push_back("AntiKt4EMPFlowJets_BTagging201903");
  //jets.push_back("AntiKt4EMPFlowJets_BTagging201810");
  //jets.push_back("AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903");
  //jets.push_back("AntiKtVR30RMax4RMin02TrackJets_BTagging201810");
  //jets.push_back("AntiKt4EMTopoJets_BTagging201810");
  
  //std::string workingPointName = "Continuous";

  // let's create matrix comparison plots for every common combo
  //*
  for(std::string tag : taggers){
  std::vector<std::string> jets;
    // let's deal with the different jet collections semi-manually here...
    if (tag == "DL1"){
      //jets.push_back("AntiKt4EMPFlowJets_BTagging201810");
      //jets.push_back("AntiKt4EMPFlowJets_BTagging201903");
      jets.push_back("AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903");
    } else if (tag == "DL1r" || tag == "DL1rmu"){
      jets.push_back("AntiKt4EMPFlowJets_BTagging201903");
      jets.push_back("AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903");
      //jets.push_back("AntiKt4EMPFlowJets");
      //jets.push_back("AntiKtVR30Rmax4Rmin02TrackJets");
    } else if (tag == "DL1dv00") {
      jets.push_back("AntiKt4EMPFlowJets");
    } else if (tag == "MV2c10") {
      jets.push_back("AntiKt4EMPFlowJets_BTagging201810");
    }
    for(std::string jeta : jets){
      for(std::string wp : wps){
        for(std::string strat : strats ){
          
          asg::StandaloneToolHandle<IBTaggingEfficiencyTool> tool("BTaggingEfficiencyTool/BTagEffTest");
          StatusCode code1 = tool.setProperty("ScaleFactorFileName","xAODBTaggingEfficiency/13TeV/2020-21-13TeV-MC16-CDI-2021-04-16_v1.root");//2019-21-13TeV-MC16-CDI-2019-10-07_v1.root" );
          //StatusCode code1 = tool.setProperty("ScaleFactorFileName","/eos/user/w/wleinone/plotting/2020-21-13TeV-MC16-CDI.root");
          StatusCode code2 = tool.setProperty("TaggerName",    tag  );
          StatusCode code3 = tool.setProperty("OperatingPoint", wp);
          StatusCode code4 = tool.setProperty("JetAuthor",     jeta);// "AntiKt4EMPFlowJets_BTagging201810" );
          StatusCode code5 = tool.setProperty("MinPt", 1);
          StatusCode code7 = tool.setProperty("SystematicsStrategy", strat);//"SFGlobalEigen");
          //tool.setProperty("EigenvectorReductionB", "Tight");
          //tool.setProperty("EigenvectorReductionC", "Tight");
          //tool.setProperty("EigenvectorReductionLight", "Tight"); // <---- does this leave 9? Or does it leave it untouched? A: Leaves it untouched.
          StatusCode code6 = tool.initialize();   
          if (code1 != StatusCode::SUCCESS || code7 != StatusCode::SUCCESS || code2 != StatusCode::SUCCESS || code3 != StatusCode::SUCCESS || code4 != StatusCode::SUCCESS || code5 != StatusCode::SUCCESS  || code6 != StatusCode::SUCCESS) 
            {
            std::cout << "Initialization of tool " << tool->name() << " failed! " << std::endl;
            return -1;
          }
          else {
            std::cout << "Initialization of tool " << tool->name() << " finished." << std::endl;
          }
  
          // select your efficiency map based on the DSID of your sample:
          //unsigned int sample_dsid = 410250;//470;

          //tool->setMapIndex(sample_dsid);

          
        } // end strats loop
      } // end wps loop
    } // end jets loop
  } // end taggers loop
  //*/
  
  

  std::cout << " Great, now the BTaggingEfficiencyToolTester is finished! " << std::endl;

  return retval;
}
