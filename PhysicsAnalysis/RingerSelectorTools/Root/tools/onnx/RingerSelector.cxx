/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

// Local includes:
#include "RingerSelectorTools/tools/onnx/RingerSelector.h"
#include "PathResolver/PathResolver.h"
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include "TEnv.h"



namespace Ringer{

  namespace onnx{


    //==============================================================================
    RingerSelector::RingerSelector(std::string name):
      asg::AsgMessaging(name)
    {;}
    

    //==============================================================================
    StatusCode RingerSelector::read_from( std::string path , AthONNX::IONNXRuntimeSvc *svc )
    {

      std::string configFile = PathResolverFindCalibFile( path );


      std::string basepath = configFile.substr(0, configFile.find_last_of("/"));

      ATH_MSG_INFO( "Loading all tunings configs from " << configFile );

      if (configFile.empty()){
         ATH_MSG_ERROR( "Could not locate the config file: " << configFile);
         return StatusCode::FAILURE; 
      }

      // Read configuration file
      TEnv env(configFile.c_str());
   

      // Retrieve all models
      {
        // Retrieve the size
        unsigned size = env.GetValue( "Model__size" , 0 );

        unsigned barcode = env.GetValue( "Model__barcode", 0); // select the input preprocessing mode
          
        // Retrieve Et bins
        auto etmin = GetValues<float>( "Model__etmin", env );
        auto etmax  = GetValues<float>( "Model__etmax", env );
        
        // Retrieve Eta bins
        auto etamin = GetValues<float>( "Model__etamin", env );
        auto etamax  = GetValues<float>( "Model__etamax", env );
        
        // Retreive all ONNX model file paths
        auto model_paths = GetPaths( "Model__path", env );
        ATH_MSG_INFO( "Basepath is "<< basepath);

        // Loop over all models
        for ( unsigned idx = 0; idx < size; ++idx )
        {
          std::string modelPath = PathResolverFindCalibFile( basepath+"/"+model_paths[idx] );
          ATH_MSG_DEBUG( "Reading Onnx model from: " << modelPath );
          auto model = Ringer::onnx::Model( modelPath, svc, etmin[idx], etmax[idx], etamin[idx], etamax[idx], barcode) ;
          
          // Compile the model
          model.compile();
          m_models.push_back(model);
        }
      }

      // Retrieve all thresholds
      {
        // Retrieve the size
         unsigned size = env.GetValue( "Threshold__size" , 0);
  
        
           
         auto max_avgmu = GetValues<float>( "Threshold__MaxAverageMu", env );
         
         // Retrieve Et bins
         auto etmin = GetValues<float>( "Threshold__etmin", env );
         auto etmax = GetValues<float>( "Threshold__etmax", env );
         // Retrieve Eta bins
         auto etamin  = GetValues<float>( "Threshold__etamin", env );
         auto etamax  = GetValues<float>( "Threshold__etamax", env );
         // Retrieve slopes
         auto slopes = GetValues<float>( "Threshold__slope", env );
         // Retrieve offsets
         auto offsets = GetValues<float>("Threshold__offset", env );

         for ( unsigned idx = 0; idx <size; ++idx ){
           m_thresholds.push_back( Ringer::onnx::Threshold(
                                     etmin[idx],
                                     etmax[idx],
                                     etamin[idx],
                                     etamax[idx],
                                     slopes[idx],
                                     offsets[idx],
                                     max_avgmu[idx])
                                 );
          }
      }

      return StatusCode::SUCCESS;
    }


   
    //==============================================================================
    bool RingerSelector::accept( const xAOD::TrigRingerRings *ringsCluster, float discr, float avgmu ) const 
    {
      float et = ringsCluster->emCluster()->et()/Gaudi::Units::GeV;
      float eta = std::abs(ringsCluster->emCluster()->eta());
    
      ATH_MSG_DEBUG( "Event et = "<< et << ", eta = " << eta );
      for( auto& cutDef : m_thresholds ){
        if ( et < cutDef.etMin() || et >= cutDef.etMax() ) continue;
        if ( eta < cutDef.etaMin() || eta >= cutDef.etaMax() ) continue;
        return cutDef.accept( discr, avgmu );
      }// loop over all thresholds
    
      return false;
    }


    //==============================================================================
    // barcode = 0: use only rings normalized by total energy as input
    // barcode = 1: use normalized rings and shower shapes (6 from cluster) as input
    // barcode = 2: use normalized rings, shower shapes and track variables as input 
    // barcode = 3: use only a half of total rings and then normalize to use as input
    std::vector<std::vector<float>> RingerSelector::prepare_inputs(  unsigned barcode,
                                                                     const xAOD::TrigRingerRings *ringsCluster, 
                                                                     const xAOD::TrigElectron *el ) const
    {
      std::vector< std::vector< float > > inputs;
      // check if the barcode = 3 because this is the unique case that we do not use 100 rings
      if ( barcode != 3 ){ 
        const std::vector<float> rings = ringsCluster->rings();
        std::vector<float> refRings(rings.size());
        refRings.assign(rings.begin(), rings.end());
        
        float energy=0.0;
        for(auto &ring : refRings ) energy+=ring;
        for(auto &ring : refRings ) ring/=std::abs(energy);
        
        inputs.push_back( refRings );
      }

      if ( barcode == 1 || barcode == 2 ){ // barcode 1
        std::vector<float> refShowers;
        const xAOD::TrigEMCluster *pClus = ringsCluster->emCluster();
        
        float e0 = pClus->energy( CaloSampling::PreSamplerB ) + pClus->energy( CaloSampling::PreSamplerE );
        float e1 = pClus->energy( CaloSampling::EMB1 ) + pClus->energy( CaloSampling::EME1 );
        float e2 = pClus->energy( CaloSampling::EMB2 ) + pClus->energy( CaloSampling::EME2 );
        float e3 = pClus->energy( CaloSampling::EMB3 ) + pClus->energy( CaloSampling::EME3 );
        float eallsamples = e0+e1+e2+e3;
        float f3 = std::abs( eallsamples )>0. ? e3/eallsamples : 0.;

        float eratio = pClus->emaxs1() - pClus->e2tsts1();
        eratio/= (pClus->emaxs1() + pClus->e2tsts1());
        float weta2 = pClus->weta2();
        float wstot = pClus->wstot();
        float reta = pClus->e237() / pClus->e277();

        // fix eratio and wstot for some cases
        if(eratio>10) eratio=0;
        if(eratio>1) eratio=1;
        if(wstot<-99) wstot=0;

        float etot=0.0;
        for ( float e : pClus->energySample() )  etot+=e;
        
        float f1 = e1/etot;

        refShowers.push_back( reta/1.);
        refShowers.push_back( eratio/1. );
        refShowers.push_back( f1/0.6 );
        refShowers.push_back( f3/0.04 );
        refShowers.push_back( weta2/0.02 );
        refShowers.push_back( wstot/1. );

        inputs.push_back( refShowers );
      }

      if ( barcode == 2 && el){ // barcode 2
        std::vector<float> refTrack;
        refTrack.push_back( el->etOverPt() );
        refTrack.push_back( el->trkClusDeta() );
        refTrack.push_back( el->trkClusDphi() );
        inputs.push_back( refTrack );
      }


      if ( barcode == 3){ // barcode 3 (half of rings)

        std::vector< std::vector< float > > inputs;
        const std::vector<float> rings = ringsCluster->rings();
        std::vector<float> refRings(rings.size());
        refRings.assign(rings.begin(), rings.end());

        std::vector<float> halfRings; 
        
        int PS   = 8;
        int EM1  = 64;
        int EM2  = 8;
        int EM3  = 8;
        int HAD1 = 4;
        int HAD2 = 4;
        int HAD3 = 4;

        halfRings.insert(halfRings.end(), refRings.begin(), refRings.begin() + PS/2);
        halfRings.insert(halfRings.end(), refRings.begin() + PS, refRings.begin() + PS + (EM1/2));
        halfRings.insert(halfRings.end(), refRings.begin() + PS + EM1, refRings.begin() + PS + EM1 + (EM2/2));
        halfRings.insert(halfRings.end(), refRings.begin() + PS + EM1 + EM2, refRings.begin() + PS + EM1 + EM2 + (EM3/2));
        halfRings.insert(halfRings.end(), refRings.begin() + PS + EM1 + EM2 + EM3, refRings.begin() + PS + EM1 + EM2 + EM3 + (HAD1/2));
        halfRings.insert(halfRings.end(), refRings.begin() + PS + EM1 + EM2 + EM3 + HAD1, refRings.begin() + PS + EM1 + EM2 + EM3 + HAD1 + (HAD2/2));
        halfRings.insert(halfRings.end(), refRings.begin() + PS + EM1 + EM2 + EM3 + HAD1 + HAD2, refRings.begin() + PS + EM1 + EM2 + EM3 + HAD1 + HAD2 + (HAD3/2));

        // concatenate
        
        
        float energy=0.0;
        for(auto &ring : halfRings ) energy+=ring;
        for(auto &ring : halfRings ) ring/=std::abs(energy);

        inputs.push_back( halfRings );
      }
      return inputs;
    }


    //==============================================================================
    float RingerSelector::predict(const xAOD::TrigRingerRings *ringsCluster , const xAOD::TrigElectron *el ) const
    {
      float et = ringsCluster->emCluster()->et()/Gaudi::Units::GeV;
      float eta = std::abs(ringsCluster->emCluster()->eta());

      // Find the correct model and predict
      for( auto& model : m_models ){
        
        if(et<model.etMin() || et >= model.etMax()) continue;
        if(eta<model.etaMin() || eta >= model.etaMax()) continue;

        auto inputs = prepare_inputs( model.barcode(), ringsCluster, el );
        auto output = model.predict( inputs ); // propagate the input throut the model
        ATH_MSG_DEBUG( "The current model predict with output: " << output );
        return output;
      }

      return -999;
    }




    //==============================================================================
    template <typename T> bool RingerSelector::strtof(const std::string& input, T& f)
    {  
      int diff = 0 ;
      std::string tmp = input;
      std::string::size_type first(0);
      std::string::size_type last(0);
      first = ( input.find('#') ) ;
 
      //if we do not find a comment character "#" we are fine
      if (first == std::string::npos) {
        std::istringstream buffer (tmp);
        buffer>>f;
        return true;
      } 
      
      //if we have found comment character check if it is inlined between two "#"
      last = (input.find('#',first+1) );
      //if nor error
      if (last == std::string::npos) {
        return false;
      }
      //else if between two "#" remove this part 
      diff = last - first ;
      tmp= tmp.erase(first,diff+1);
      std::istringstream buffer (tmp);
      buffer>>f;
      return true;
    }


    //==============================================================================
    template <typename T>  std::vector<T> RingerSelector::GetValues (const std::string& input,  TEnv& env)
    {    
      std::vector<T> CutVector;    
      std::string env_input(env.GetValue(input.c_str(), ""));
      if (!env_input.empty()) {
        std::string::size_type end;
        do {
          end = env_input.find(';');
          T  myValue(0);
          if(strtof(env_input.substr(0,end),myValue)){
            CutVector.push_back(myValue);
          }
          if (end != std::string::npos) {
            env_input= env_input.substr(end+1);
          }
        } while (end != std::string::npos);     
      }
      return CutVector;
    }


    //==============================================================================
    std::vector<std::string> RingerSelector::GetPaths(const std::string& input, TEnv& env)
    {
      std::vector<std::string> CutVector;    
      std::string env_input(env.GetValue(input.c_str(), ""));
      if (!env_input.empty()) {
        std::string::size_type end;
        do {
          end = env_input.find(';');
          CutVector.push_back( env_input.substr(0,end) );
          if (end != std::string::npos) {
            env_input= env_input.substr(end+2);
          }
        } while (end != std::string::npos);     
      }
      return CutVector;
    }




  } // namespace onnx
} // namespace Ringer


