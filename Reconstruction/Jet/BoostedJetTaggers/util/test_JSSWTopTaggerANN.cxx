/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

// System include(s):
#include <string>

// ROOT include(s):
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

// EDM include(s):
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/tools/IOStats.h"

// Tool testing include(s):
#include "BoostedJetTaggers/JSSWTopTaggerANN.h"
#include "JetUncertainties/JetUncertaintiesTool.h"

int main( int argc, char* argv[] ) {

  // The application's name:
  char* APP_NAME = argv[ 0 ];

  // arguments
  TString fileName = "/eos/user/d/dmelini/MassDecorrelatedTagger/test_samples/mc16_13TeV.364707.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ7WithSW.deriv.DAOD_JETM8.e7142_s3126_r10201_p4049/DAOD_JETM8.20146864._000047.pool.root.1";
  TString calibArea="/afs/cern.ch/work/d/dmelini/MassDecorrelatedTaggers/toBoostedJetTagger/UFO_Jan22/ANN/";
  TString configFileName = "JSSANN50Tagger_AntiKt10UFOCSSKSoftDrop_Jan22.dat";
  int  ievent=-1;
  int  nevents=-1;
  bool m_isMC=true;
  bool verbose=false;
  bool m_calcSF=false;
  Info( APP_NAME, "==============================================" );
  Info( APP_NAME, "Usage: $> %s [xAOD file name]", APP_NAME );
  Info( APP_NAME, " $> %s       | Run on default file", APP_NAME );
  Info( APP_NAME, " $> %s -f X  | Run on xAOD file X", APP_NAME );
  Info( APP_NAME, " $> %s -a X  | Run in calibArea X", APP_NAME );
  Info( APP_NAME, " $> %s -c X  | Run with config file X", APP_NAME );
  Info( APP_NAME, " $> %s -n X  | X = number of events you want to run on", APP_NAME );
  Info( APP_NAME, " $> %s -e X  | X = specific number of the event to run on - for debugging", APP_NAME );
  Info( APP_NAME, " $> %s -d X  | X = dataset ID", APP_NAME );
  Info( APP_NAME, " $> %s -m X  | X = isMC", APP_NAME );
  Info( APP_NAME, " $> %s -v    | run in verbose mode   ", APP_NAME );
  Info( APP_NAME, " $> %s -sf   | run also ScaleFactor  ", APP_NAME );
  Info( APP_NAME, "==============================================" );

  // Check if we received a file name:
  if( argc < 2 ) {
    Info( APP_NAME, "No arguments - using default file" );
    Info( APP_NAME, "Executing on : %s", fileName.Data() );
  }

  ////////////////////////////////////////////////////
  //:::  parse the options
  ////////////////////////////////////////////////////
  std::string options;
  for( int i=0; i<argc; i++){
    options+=(argv[i]);
  }

  if(options.find("-f")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-f")==0){
        fileName = argv[ipos+1];
        Info( APP_NAME, "Argument (-f) : Running on file # %s", fileName.Data() );
        break;
      }
    }
  }

  if(options.find("-a")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-a")==0){
        calibArea= argv[ipos+1];
        Info( APP_NAME, "Argument (-a) : Running in calibArea # %s", calibArea.Data() );
        break;
      }
    }
  }

  if(options.find("-c")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-c")==0){
        configFileName = argv[ipos+1];
        Info( APP_NAME, "Argument (-c) : Running with config # %s", configFileName.Data() );
        break;
      }
    }
  }

  if(options.find("-event")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-event")==0){
        ievent = atoi(argv[ipos+1]);
        Info( APP_NAME, "Argument (-event) : Running only on event # %i", ievent );
        break;
      }
    }
  }

  if(options.find("-m")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-m")==0){
        m_isMC = atoi(argv[ipos+1]);
        Info( APP_NAME, "Argument (-m) : IsMC = %i", m_isMC );
        break;
      }
    }
  }

  if(options.find("-n")!=std::string::npos){
    for( int ipos=0; ipos<argc ; ipos++ ) {
      if(std::string(argv[ipos]).compare("-n")==0){
        nevents = atoi(argv[ipos+1]);
        Info( APP_NAME, "Argument (-n) : Running on NEvents = %i", nevents );
        break;
      }
    }
  }

  if(options.find("-v")!=std::string::npos){
    verbose=true;
    Info( APP_NAME, "Argument (-v) : Setting verbose");
  }

  if(options.find("-sf")!=std::string::npos){
    m_calcSF=true;
    Info( APP_NAME, "Argument (-sf) : Also calculating scale factors");
  }


  ////////////////////////////////////////////////////
  //:::  initialize the application and get the event
  ////////////////////////////////////////////////////
  xAOD::Init( APP_NAME );
  xAOD::TReturnCode::enableFailure();

  // Open the input file:
  TFile* ifile( TFile::Open( fileName, "READ" ) );
  if( !ifile ) Error( APP_NAME, "Cannot find file %s",fileName.Data() );

  TChain *chain = new TChain ("CollectionTree","CollectionTree");
  chain->Add(fileName);

  // Create a TEvent object:
  xAOD::TEvent event( (TTree*)chain, xAOD::TEvent::kAthenaAccess );
  Info( APP_NAME, "Number of events in the file: %i", static_cast< int >( event.getEntries() ) );

  // Create a transient object store. Needed for the tools.
  xAOD::TStore store;

  // Decide how many events to run over:
  Long64_t entries = event.getEntries();

  // Fill a validation true with the tag return value
  std::unique_ptr<TFile> outputFile(TFile::Open("output_JSSWTopTaggerANN.root", "recreate"));
  int pass,truthLabel,idx;
  float sf,pt,eta,m;
  TTree* Tree = new TTree( "tree", "test_tree" );
  Tree->Branch( "pass", &pass, "pass/I" );
  Tree->Branch( "sf", &sf, "sf/F" );
  Tree->Branch( "pt", &pt, "pt/F" );
  Tree->Branch( "m", &m, "m/F" );
  Tree->Branch( "eta", &eta, "eta/F" );
  Tree->Branch( "idx", &idx, "idx/I" );
  Tree->Branch( "truthLabel", &truthLabel, "truthLabel/I" );

  
  std::unique_ptr<JetUncertaintiesTool> m_jetUncToolSF(new JetUncertaintiesTool(("JetUncProvider_SF")));
  std::vector<CP::SystematicSet> m_jetUnc_sysSets;

  if(m_calcSF){
    m_jetUncToolSF->setProperty("JetDefinition", "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets");
    m_jetUncToolSF->setProperty("Path", "/eos/user/g/gang/public/BoostedJetTaggers/JSSWTopTaggerANN/");
    m_jetUncToolSF->setProperty("ConfigFile", "TagSFUncert_JSSANNTagger_AntiKt10LCTopoTrimmed.config");
    m_jetUncToolSF->setProperty("MCType", "MC16a");
    m_jetUncToolSF->initialize();
    std::vector<std::string> pulls = {"__1down", "__1up"};
    CP::SystematicSet jetUnc_sysSet = m_jetUncToolSF->recommendedSystematics();
    const std::set<std::string> sysNames = jetUnc_sysSet.getBaseNames();
    for (std::string sysName: sysNames) {
      for (std::string pull : pulls) {
	std::string sysPulled = sysName + pull;
	m_jetUnc_sysSets.push_back(CP::SystematicSet(sysPulled));
      }
    }
  }
  

  ////////////////////////////////////////////
  /////////// START TOOL SPECIFIC ////////////
  ////////////////////////////////////////////

  ////////////////////////////////////////////////////
  //::: Tool setup
  // setup the tool handle as per the
  // recommendation by ASG - https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AthAnalysisBase#How_to_use_AnaToolHandle
  ////////////////////////////////////////////////////
  std::cout<<"Initializing JSSWTopTaggerANN Tagger"<<std::endl;
  asg::AnaToolHandle<JSSWTopTaggerANN> m_Tagger; //!
  ASG_SET_ANA_TOOL_TYPE( m_Tagger, JSSWTopTaggerANN);
  m_Tagger.setName("myANN50Tagger");
  if(verbose) m_Tagger.setProperty("OutputLevel", MSG::DEBUG);
  m_Tagger.setProperty( "CalibArea",    calibArea);
  m_Tagger.setProperty( "ConfigFile",   configFileName);
  if(m_calcSF)
    m_Tagger.setProperty("TruthJetContainerName", "AntiKt10TruthSoftDropBeta100Zcut10Jets");
  m_Tagger.setProperty("IsMC", m_isMC);
  m_Tagger.retrieve();


  std::cout << "Total Events in File : " << entries << std::endl;

  ////////////////////////////////////////////////////
  // Loop over the events
  ////////////////////////////////////////////////////
  for( Long64_t entry = 0; entry < entries; ++entry ) {

    if( nevents!=-1 && entry > nevents ) break;
    // Tell the object which entry to look at:
    event.getEntry( entry );

    // Print some event information
    const xAOD::EventInfo* evtInfo = 0;
    if(event.retrieve( evtInfo, "EventInfo" ) != StatusCode::SUCCESS){
      continue;
    }
    if(ievent!=-1 && static_cast <int> (evtInfo->eventNumber())!=ievent) {
      continue;
    }

    // Get the jets
    const xAOD::JetContainer* myJets = 0;
    if( event.retrieve( myJets, "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets" ) != StatusCode::SUCCESS)
      continue ;

    // Loop over jet container
    std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > jets_shallowCopy = xAOD::shallowCopyContainer( *myJets );
    std::unique_ptr<xAOD::JetContainer> shallowJets(jets_shallowCopy.first);
    std::unique_ptr<xAOD::ShallowAuxContainer> shallowAux(jets_shallowCopy.second);
    idx=0;
    for( xAOD::Jet* jetSC : *shallowJets ){      
      if(verbose) std::cout<<"Testing ANN Tagger "<<std::endl;
      const Root::TAccept& res = m_Tagger->tag( *jetSC );
      if(verbose) std::cout<<"jet pt              = "<<jetSC->pt()<<std::endl;
      if(verbose) std::cout<<"RunningTag : "<<res<<std::endl;
      if(verbose) std::cout<<"Printing jet score : " << jetSC->auxdata<float>("ANNTagger_Score") << std::endl;
      if(verbose) std::cout<<"result masspasslow  = "<<res.getCutResult("PassMassLow")<<std::endl;
      if(verbose) std::cout<<"result masspasshigh = "<<res.getCutResult("PassMassHigh")<<std::endl;
      truthLabel = jetSC->auxdata<int>("R10TruthLabel_R21Consolidated");

      pass = res;
      sf = jetSC->auxdata<float>("ANNTagger_SF");
      pt = jetSC->pt();
      m  = jetSC->m();
      eta = jetSC->eta();
      
      Tree->Fill();
      idx++;

      if ( m_isMC && m_calcSF){
	if ( pt/1.e3 > 350 && std::abs(eta) < 2.0 && pass ) {
	  bool validForUncTool = ( pt/1.e3 >= 150 && pt/1.e3 < 2500 );
	  validForUncTool &= ( m/pt >= 0 && m/pt <= 1 );
	  validForUncTool &= ( std::abs(eta) < 2 );
	  std::cout << "Nominal SF=" << sf << " truthLabel=" << truthLabel << " (1: t->qqb)" << std::endl;
	  if( validForUncTool ){
	    for ( CP::SystematicSet sysSet : m_jetUnc_sysSets ){
	      m_Tagger->tag( *jetSC );
	      m_jetUncToolSF->applySystematicVariation(sysSet);
	      m_jetUncToolSF->applyCorrection(*jetSC);
	      std::cout << sysSet.name() << " " << jetSC->auxdata<float>("ANNTagger_SF") << std::endl;
	    }
	  }
	  
	}
      }
    }

  Info( APP_NAME, "===>>>  done processing event #%i, run #%i %i events processed so far  <<<===", static_cast< int >( evtInfo->eventNumber() ), static_cast< int >( evtInfo->runNumber() ), static_cast< int >( entry + 1 ) );
  }

  ////////////////////////////////////////////
  /////////// END TOOL SPECIFIC //////////////
  ////////////////////////////////////////////

  // write the tree to the output file
  outputFile->cd();
  Tree->Write();
  outputFile->Close();

  // cleanup
  delete chain;

  // print the branches that were used for help with smart slimming
  std::cout<<std::endl<<std::endl;
  std::cout<<"Smart Slimming Checker :"<<std::endl;
  xAOD::IOStats::instance().stats().printSmartSlimmingBranchList();
  std::cout<<std::endl<<std::endl;

  return 0;

}

