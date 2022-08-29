/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "xAODBTaggingEfficiency/BTaggingSelectionTool.h"
#include "xAODBTagging/BTagging.h"
#include "xAODBTagging/BTaggingUtilities.h"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include "CalibrationDataInterface/CalibrationDataVariables.h"
#include "CalibrationDataInterface/CalibrationDataContainer.h"

#include "PATInterfaces/SystematicRegistry.h"
#include "PathResolver/PathResolver.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <algorithm>
#include <string>

using std::string;

using CP::CorrectionCode;
using CP::SystematicSet;
using CP::SystematicVariation;
using CP::SystematicRegistry;

using Analysis::Uncertainty;
using Analysis::CalibrationDataVariables;
using Analysis::CalibrationDataContainer;
using Analysis::CalibResult;
using Analysis::CalibrationStatus;
using Analysis::Total;
using Analysis::SFEigen;
using Analysis::SFNamed;
using Analysis::None;

using xAOD::IParticle;

BTaggingSelectionTool::BTaggingSelectionTool( const std::string & name)
  : asg::AsgTool( name ), m_acceptinfo( "JetSelection" )
{
  m_initialised = false;
  declareProperty( "MaxEta", m_maxEta = 2.5 );
  declareProperty( "MinPt", m_minPt = -1 /*MeV*/);
  declareProperty( "MaxRangePt", m_maxRangePt = 3000000 /*MeV*/);
  declareProperty( "FlvTagCutDefinitionsFileName", m_CutFileName = "", "name of the files containing official cut definitions (uses PathResolver)");
  declareProperty( "TaggerName",                    m_taggerName="",    "tagging algorithm name");
  declareProperty( "OperatingPoint",                m_OP="",            "operating point");
  declareProperty( "JetAuthor",                     m_jetAuthor="",     "jet collection");
  declareProperty( "ErrorOnTagWeightFailure",     m_ErrorOnTagWeightFailure=true, "optionally ignore cases where the tagweight cannot be retrived. default behaviour is to give an error, switching to false will turn it into a warning");
}

StatusCode BTaggingSelectionTool::initialize() {
  m_initialised = true;

  if (""==m_OP){
    ATH_MSG_ERROR( "BTaggingSelectionTool wasn't given a working point name" );
    return StatusCode::FAILURE;
  }

  TString pathtofile =  PathResolverFindCalibFile(m_CutFileName);
  m_inf = TFile::Open(pathtofile, "read");
  if (0==m_inf) {
    ATH_MSG_ERROR( "BTaggingSelectionTool couldn't access tagging cut definitions" );
    return StatusCode::FAILURE;
  }

  m_useVeto = false;
  if( m_OP.find("Veto") != string::npos ){
    m_useVeto = true;
    //op string should follow the format nominal_WP_Veto_vetotagger_vetoWP
    //for example, FixedCutBEff_77_Veto_DL1_CTag_Loose
    TString vetotaggerstring = m_OP;
    TObjArray* wptokens = vetotaggerstring.Tokenize("_");
    if( wptokens->GetEntries() != 6 ){
        ATH_MSG_ERROR( "BTaggingSelectionTool improperly formatted WP defintion: "+m_OP );
        return StatusCode::FAILURE;
    }
    m_OP = (std::string)(((TObjString *)(wptokens->At(0)))->String());
    m_OP = m_OP+"_"+(std::string)(((TObjString *)(wptokens->At(1)))->String());

    m_taggerName_Veto = (std::string)(((TObjString *)(wptokens->At(3)))->String());
    m_OP_Veto = (std::string)(((TObjString *)(wptokens->At(4)))->String());
    m_OP_Veto = m_OP_Veto+"_"+(std::string)(((TObjString *)(wptokens->At(5)))->String());
  }

  // check the CDI file for the selected tagger and jet collection
  TString check_CDI = m_taggerName;
  if(!m_inf->Get(check_CDI)){
     ATH_MSG_ERROR( "Tagger: "+m_taggerName+" not found in this CDI file: "+m_CutFileName);
    return StatusCode::FAILURE;
  }
  check_CDI = m_taggerName+"/"+m_jetAuthor;
  if(!m_inf->Get(check_CDI)){
     ATH_MSG_ERROR( "Tagger: "+m_taggerName+" and Jet Collection : "+m_jetAuthor+"  not found in this CDI file: "+m_CutFileName);
    return StatusCode::FAILURE;
  }


  // Change the minPt cut if the user didn't touch it
  if (m_minPt < 0) {
    ATH_MSG_ERROR( "Tagger: "+m_taggerName+" and Jet Collection : "+m_jetAuthor+" do not have a minimum jet pT cut set.");
    return StatusCode::FAILURE;
  }

  // Operating point reading
    TString cutname = m_OP;

  if ("Continuous"==cutname(0,10)){  // For continuous tagging load all flat-cut WPs
    //100% efficiency => MVXWP=-infinity
    m_continuouscuts[0] = -1.e4;

    cutname = m_taggerName+"/"+m_jetAuthor+"/FixedCutBEff_85/cutvalue";
    m_tagger.constcut = (TVector*) m_inf->Get(cutname);
    if (m_tagger.constcut!=nullptr) m_continuouscuts[1] = m_tagger.constcut[0](0);
    else ATH_MSG_ERROR( "Continuous tagging is trying to use an unvalid operating point: FixedCutBEff_85" );

    cutname = m_taggerName+"/"+m_jetAuthor+"/FixedCutBEff_77/cutvalue";
    m_tagger.constcut = (TVector*) m_inf->Get(cutname);
    if (m_tagger.constcut!=nullptr) m_continuouscuts[2] = m_tagger.constcut[0](0);
    else ATH_MSG_ERROR( "Continuous tagging is trying to use an unvalid operating point: FixedCutBEff_77" );

    cutname = m_taggerName+"/"+m_jetAuthor+"/FixedCutBEff_70/cutvalue";
    m_tagger.constcut = (TVector*) m_inf->Get(cutname);
    if (m_tagger.constcut!=nullptr) m_continuouscuts[3] = m_tagger.constcut[0](0);
    else ATH_MSG_ERROR( "Continuous tagging is trying to use an unvalid operating point: FixedCutBEff_70" );

    cutname = m_taggerName+"/"+m_jetAuthor+"/FixedCutBEff_60/cutvalue";
    m_tagger.constcut = (TVector*) m_inf->Get(cutname);
    if (m_tagger.constcut!=nullptr) m_continuouscuts[4] = m_tagger.constcut[0](0);
    else ATH_MSG_ERROR( "Continuous tagging is trying to use an unvalid operating point: FixedCutBEff_60" );

    //0% efficiency => MVXWP=+infinity
    m_continuouscuts[5]= +1.e4;

    if(m_taggerName.find("DL1") != string::npos || m_taggerName.find("GN") != string::npos){
      //this call will extract the c-fraction value and set it in m_tagger
      //which is needed in case the user calls getTaggerWeight to compute the DL1 score or the GNx scores
      ExtractTaggerProperties(m_tagger,m_taggerName , "FixedCutBEff_60");
    }
  }
  else {  // Else load only one WP
    ExtractTaggerProperties(m_tagger,m_taggerName , m_OP);
    if(m_useVeto){ExtractTaggerProperties(m_vetoTagger,m_taggerName_Veto , m_OP_Veto);}
  }
  m_inf->Close();

  m_acceptinfo.addCut( "Eta", "Selection of jets according to their pseudorapidity" );
  m_acceptinfo.addCut( "Pt",  "Selection of jets according to their transverse momentum" );
  m_acceptinfo.addCut( "WorkingPoint",  "Working point for flavour-tagging of jets according to their b-tagging weight" );

  return StatusCode::SUCCESS;
}

void BTaggingSelectionTool::ExtractTaggerProperties(taggerproperties &tagger, std::string taggerName, std::string OP){

  TString cutname = OP;

  if ("FlatBEff"==cutname(0,8) || "HybBEff"==cutname(0,7) ){
    cutname = taggerName+"/"+m_jetAuthor+"/"+OP+"/cutprofile";
    tagger.spline = (TSpline3*) m_inf->Get(cutname);
    if (tagger.spline == nullptr) ATH_MSG_ERROR( "Invalid operating point" );
    tagger.constcut = nullptr;
  }
  else {
    cutname = taggerName+"/"+m_jetAuthor+"/"+OP+"/cutvalue";
    tagger.constcut = (TVector*) m_inf->Get(cutname);
    if (tagger.constcut == nullptr) ATH_MSG_ERROR( "Invalid operating point" );
    tagger.spline = nullptr;
  }

  //retrive the "fraction" used in the DL1/GNx log likelihood from the CDI, if its not there, use the hard coded values
  // (backwards compatibility)
  if(taggerName.find("DL1") != string::npos || taggerName.find("GN") != string::npos){

    TString fraction_data_name = taggerName+"/"+m_jetAuthor+"/"+OP+"/fraction";
    TVector *fraction_data = (TVector*) m_inf->Get(fraction_data_name);

    if(fraction_data!=nullptr){
      tagger.fraction = fraction_data[0](0);
    }else{
      if("DL1"    ==taggerName){ tagger.fraction = 0.08; }
      if("DL1mu"  ==taggerName){ tagger.fraction = 0.08; }
      if("DL1rnn" ==taggerName){ tagger.fraction = 0.03; }

      if(OP == "CTag_Loose" && "DL1"    == taggerName ){ tagger.fraction = 0.25; }
      if(OP == "CTag_Tight" && "DL1"    == taggerName ){ tagger.fraction = 0.18; }
      if(OP == "CTag_Loose" && "DL1rnn" == taggerName ){ tagger.fraction = 0.08; }
      if(OP == "CTag_Tight" && "DL1rnn" == taggerName ){ tagger.fraction = 0.02; }
    }
    delete fraction_data;
  }
}


CorrectionCode BTaggingSelectionTool::getTaggerWeight( const xAOD::Jet& jet, double & tagweight, bool useVetoWP) const{

  taggerproperties localtagger;
  std::string taggerName;
  std::string OP;

  if(useVetoWP && m_useVeto){
    localtagger = m_vetoTagger;
    taggerName = m_taggerName_Veto;
    OP = m_OP_Veto;
  }else{
    localtagger = m_tagger;
    taggerName = m_taggerName;
    OP = m_OP;
  }

 tagweight = -100.;

 if ( "MV2cl100_MV2c100"==taggerName ){

  ATH_MSG_ERROR("Cannot retrive weight for 2D taggers!");

  return CorrectionCode::Error;
 }

 else if ( taggerName.find("MV2") != string::npos ){


    const xAOD::BTagging* btag = xAOD::BTaggingUtilities::getBTagging( jet );

    if ((!btag) || (!btag->MVx_discriminant(taggerName, tagweight))){
      if(m_ErrorOnTagWeightFailure){
        ATH_MSG_ERROR("Failed to retrieve "+taggerName+" weight!");
        return CorrectionCode::Error;
      }else{
        ATH_MSG_WARNING("Failed to retrieve "+taggerName+" weight!");
        return CorrectionCode::Ok;
      }
    }
    ATH_MSG_VERBOSE( taggerName << " " <<  tagweight );
    return  CorrectionCode::Ok;
  }

 else if(taggerName.find("DL1") != string::npos || taggerName.find("GN") != string::npos){

  double dl1_pb(-10.);
  double dl1_pc(-10.);
  double dl1_pu(-10.);

  const xAOD::BTagging* btag = xAOD::BTaggingUtilities::getBTagging( jet );

  if ((!btag)){
   ATH_MSG_ERROR("Failed to retrieve the BTagging information");
  }

  if ( (!btag->pb(taggerName, dl1_pb ))
   || (!btag->pc(taggerName, dl1_pc ))
   || (!btag->pu(taggerName, dl1_pu )) ){

     if(m_ErrorOnTagWeightFailure){
       ATH_MSG_ERROR("Failed to retrieve "+taggerName+" weight!");
       return CorrectionCode::Error;
     }else{
       ATH_MSG_WARNING("Failed to retrieve "+taggerName+" weight!");
       return CorrectionCode::Ok;
     }
  }
  ATH_MSG_VERBOSE( "pb " <<  dl1_pb );
  ATH_MSG_VERBOSE( "pc " <<  dl1_pc );
  ATH_MSG_VERBOSE( "pu " <<  dl1_pu );

  return getTaggerWeight(dl1_pb, dl1_pc, dl1_pu, tagweight, useVetoWP);

  }

  //if we got here the tagger name is not configured properly
  ATH_MSG_ERROR("BTaggingSelectionTool doesn't support tagger: "+m_taggerName);
  return CorrectionCode::Error;

}

CorrectionCode BTaggingSelectionTool::getTaggerWeight( double pb, double pc, double pu , double & tagweight,  bool useVetoWP) const {

  taggerproperties localtagger;
  std::string taggerName;
  std::string OP;

  if(useVetoWP && m_useVeto){
    localtagger = m_vetoTagger;
    taggerName = m_taggerName_Veto;
    OP = m_OP_Veto;
  }else{
    localtagger = m_tagger;
    taggerName = m_taggerName;
    OP = m_OP;
  }


  tagweight = -100.;
  if( taggerName.find("DL1") != string::npos || m_taggerName.find("GN") != string::npos){

    bool valid_input = (!std::isnan(pu) && pb>=0 && pc>=0 && pu>=0);

    if (!valid_input){
      if(m_ErrorOnTagWeightFailure){
        ATH_MSG_ERROR("Invalid inputs for "+taggerName+" pb " << pb << " pc " << pc << " pu " << pu << " ");
        return CorrectionCode::Error;
      }else{
        ATH_MSG_WARNING("Invalid inputs for "+taggerName+" pb " << pb << " pc " << pc << " pu " << pu << " ");
        return CorrectionCode::Ok;
      }
    }

    if(OP.find("CTag") != string::npos){
     tagweight = log(pc / (localtagger.fraction * pb + (1. - localtagger.fraction) * pu));
    }
    else{
     tagweight = log(pb / (localtagger.fraction * pc + (1. - localtagger.fraction) * pu) );
    }

    return CorrectionCode::Ok;
  }

  //if we got here the tagger name is not configured properly
  ATH_MSG_ERROR("this call to getTaggerWeight only works for DL1/GNx taggers");
  return CorrectionCode::Error;

}



asg::AcceptData BTaggingSelectionTool::accept( const xAOD::IParticle* p ) const {

  // Check if this is a jet:
  if( p->type() != xAOD::Type::Jet ) {
    ATH_MSG_ERROR( "accept(...) Function received a non-jet" );
    return asg::AcceptData (&m_acceptinfo);
  }

  // Cast it to a jet:
  const xAOD::Jet* jet = dynamic_cast< const xAOD::Jet* >( p );
  if( ! jet ) {
    ATH_MSG_FATAL( "accept(...) Failed to cast particle to jet" );
    return asg::AcceptData (&m_acceptinfo);
  }

  // Let the specific function do the work:
  return accept( *jet );
}

asg::AcceptData BTaggingSelectionTool::accept( const xAOD::Jet& jet ) const {
  asg::AcceptData acceptData (&m_acceptinfo);

  if (! m_initialised) {
    ATH_MSG_ERROR("BTaggingSelectionTool has not been initialised");
    return acceptData;
  }

  if  ( "AntiKt2PV0TrackJets"== m_jetAuthor ||
        "AntiKt4PV0TrackJets"== m_jetAuthor ||
        m_jetAuthor.find("AntiKtVR30Rmax4Rmin02TrackJets") != string::npos
      ) {
    // We want at least 2 tracks in a track jet
    acceptData.setCutResult( "NConstituents", jet.numConstituents() >= 2 );
  }

  double pT = jet.pt();
  double eta = jet.eta();

  // Retrieve the tagger weight which was assigned to the jet
  TString cutname = m_OP;

  //the MV2 c-tagger is a special case, its a 2D cut on two MV2 based taggers
  if(m_taggerName=="MV2cl100_MV2c100"){

    double weight_mv2c100(-10.);
    double weight_mv2cl100(-10.);
    const xAOD::BTagging* btag = xAOD::BTaggingUtilities::getBTagging( jet );

    if ((!btag)){
      ATH_MSG_ERROR("Failed to retrieve the BTagging information");
    }

    if ( (!btag->MVx_discriminant("MV2cl100", weight_mv2cl100 ))|| (!btag->MVx_discriminant("MV2c100", weight_mv2c100 )) ){
      if(m_ErrorOnTagWeightFailure){
        ATH_MSG_ERROR("Failed to retrieve the BTagging "+m_taggerName+" weight!");
      }else{
        ATH_MSG_WARNING("Failed to retrieve the BTagging "+m_taggerName+" weight!");
      }
      return acceptData;
    }
    ATH_MSG_VERBOSE( "MV2cl100 "  <<  weight_mv2cl100 );
    ATH_MSG_VERBOSE( "MV2c100 " <<  weight_mv2c100 );
    return accept(pT, eta, weight_mv2cl100, weight_mv2c100 );

  }else if(m_taggerName.find("DL1") != string::npos || m_taggerName.find("MV2") != string::npos || m_taggerName.find("GN") != string::npos){


    //for all other taggers, use the same method
    double tagger_weight(-100);

    if( getTaggerWeight( jet ,tagger_weight)!=CorrectionCode::Ok){
      return acceptData;
    };

      if(m_useVeto){
        double taggerweightVeto(-100);
        if( getTaggerWeight( jet , taggerweightVeto, true)!=CorrectionCode::Ok){
          return acceptData;
        };
        return accept(pT, eta, tagger_weight,taggerweightVeto);
      }

      return accept(pT, eta, tagger_weight);
  }

    //if we got here the tagger name is not configured properly
    ATH_MSG_ERROR("BTaggingSelectionTool doesn't support tagger: "+m_taggerName);
    return acceptData;

}

asg::AcceptData BTaggingSelectionTool::accept(double pT, double eta, double tag_weight) const
{
  
  asg::AcceptData acceptData (&m_acceptinfo);

  if (! m_initialised) {
    ATH_MSG_ERROR("BTaggingSelectionTool has not been initialised");
    return acceptData;
  }

  eta = std::abs(eta);

  if (! checkRange(pT, eta,acceptData))
    return acceptData;

  // After initialization, either m_tagger.spline or m_tagger.constcut should be non-zero
  // Else, the initialization was incorrect and should be revisited
  double cutvalue(DBL_MAX);
  if( getCutValue(pT, cutvalue )!=CorrectionCode::Ok ){
    return acceptData;
  };

  if (  tag_weight < cutvalue ){
    return acceptData;
  }
  
  acceptData.setCutResult( "WorkingPoint", true );
  // Return the result:
  return acceptData;
}



asg::AcceptData BTaggingSelectionTool::accept(double pT, double eta, double weight_mv2cl100, double weight_mv2c100) const
{
   asg::AcceptData acceptData (&m_acceptinfo);

  if (! m_initialised) {
    ATH_MSG_ERROR("BTaggingSelectionTool has not been initialised");
    return acceptData;
  }

  eta = std::abs(eta);

  if (! checkRange(pT, eta,acceptData))
    return acceptData;

  // After initialization, either m_tagger.spline or m_tagger.constcut should be non-zero
  // Else, the initialization was incorrect and should be revisited
  double cutvalueA(DBL_MAX);
  double cutvalueB(DBL_MAX);


  if(m_useVeto){

    double nominal_taggerWeight = weight_mv2cl100;
    double veto_taggerWeight = weight_mv2c100;

   if( getCutValue(pT, cutvalueA)!=CorrectionCode::Ok){
    return acceptData;
   };
   if( getCutValue(pT, cutvalueB, true)!=CorrectionCode::Ok){
    return acceptData;
   };

      if( ! ( nominal_taggerWeight > cutvalueA && veto_taggerWeight <= cutvalueB ) ){
          return acceptData;
      }

  }else if(m_tagger.constcut != nullptr && m_tagger.spline == nullptr) {
    cutvalueA = std::as_const(m_tagger.constcut)[0](0);
    cutvalueB = std::as_const(m_tagger.constcut)[0](1);

    ATH_MSG_VERBOSE( "Cut values " << cutvalueA << " " << cutvalueB );

    if ( !(  weight_mv2cl100 > cutvalueA && weight_mv2c100 < cutvalueB )  ){
      return acceptData;
    }
  }
  else{
    ATH_MSG_ERROR( "Bad cut configuration!" );
    return acceptData;
  }

  acceptData.setCutResult( "WorkingPoint", true ); // IF we arrived here, the jet is tagged

  // Return the result:
  return acceptData;
}


asg::AcceptData BTaggingSelectionTool::accept(double pT, double eta, double pb, double pc, double pu) const
 {
   asg::AcceptData acceptData (&m_acceptinfo);

   if (! m_initialised) {
     ATH_MSG_ERROR("BTaggingSelectionTool has not been initialised");
     return acceptData;
   }

   eta = std::abs(eta);

   if (! checkRange(pT, eta,acceptData))
     return acceptData;

   // After initialization, either m_tagger.spline or m_tagger.constcut should be non-zero
   // Else, the initialization was incorrect and should be revisited
   double cutvalue(DBL_MAX);

   if( getCutValue(pT, cutvalue)!=CorrectionCode::Ok){
    return acceptData;
   };


   double tagger_weight(-100);

   if( getTaggerWeight(pb, pc, pu, tagger_weight)!=CorrectionCode::Ok){
      return acceptData;
   }

   if ( tagger_weight < cutvalue ){
     return acceptData;
   }

   //if you made it here, the jet is tagged
   acceptData.setCutResult( "WorkingPoint", true );

   return acceptData;
 }

int BTaggingSelectionTool::getQuantile( const xAOD::IParticle* p ) const {
  // Check if this is a jet:


  if( p->type() != xAOD::Type::Jet ) {
    ATH_MSG_ERROR( "accept(...) Function received a non-jet" );
    return -1;
  }

  // Cast it to a jet:
  const xAOD::Jet* jet = dynamic_cast< const xAOD::Jet* >( p );
  if( ! jet ) {
    ATH_MSG_FATAL( "accept(...) Failed to cast particle to jet" );
    return -1;
  }

  // Let the specific function do the work:
  return getQuantile( *jet );
}

int BTaggingSelectionTool::getQuantile( const xAOD::Jet& jet ) const{
  double pT = jet.pt();
  double eta = std::abs( jet.eta() );
  // Retrieve the tagger weight which was assigned to the jet
  double weight_mv2(-10.);
  if (getTaggerWeight(jet, weight_mv2)==CorrectionCode::Error){
    ATH_MSG_WARNING("getQuantile: Failed to retrieve "+m_taggerName+" weight!");
    return -1;
  }
  ATH_MSG_VERBOSE( m_taggerName << " " <<  weight_mv2 );

  return getQuantile(pT, eta, weight_mv2 );
}


int BTaggingSelectionTool::getQuantile(double pT, double eta, double weight_mv2 ) const
{
  if (! m_initialised) {
    ATH_MSG_ERROR("BTaggingSelectionTool has not been initialised");
  }
  //////////////////////
  // Cheatsheet:
  // returns 5 if between 60% and 0%
  // returns 4 if between 70% and 60%
  // returns 3 if between 77% and 70%
  // returns 2 if between 85% and 77%
  // returns 1 if between 100% and 85%
  // return -1 not in b-tagging acceptance
  //////////////////////

  int bin_index(-1);
  asg::AcceptData acceptData (&m_acceptinfo);
  if (! checkRange(pT, eta,acceptData)) return bin_index;

  // If in b-tagging acceptance, cont.tagging
  for (int i=1; i<=5; ++i) {
    if (weight_mv2 < m_continuouscuts[i]) {
      bin_index = i;
      break;
    }
    else if (weight_mv2 >= m_continuouscuts[5]){
      bin_index = 5;
      break;
    }
  }
  return bin_index;
}

bool BTaggingSelectionTool::checkRange(double pT, double eta,asg::AcceptData& acceptData) const
{
  // Do the |eta| cut:
  if( std::abs(eta) > m_maxEta ) {
    return false;
  }
  acceptData.setCutResult( "Eta", true );

  // Do the pT cut:
  ATH_MSG_VERBOSE( "Jet pT: " << pT );
  if( pT < m_minPt ) {
    return false;
  }
  acceptData.setCutResult( "Pt", true );

  return true;
}

CorrectionCode BTaggingSelectionTool::getCutValue(double pT, double & cutval,bool useVetoWP) const
{
   cutval = DBL_MAX;

  taggerproperties localtagger;

  if(useVetoWP && m_useVeto){
    localtagger = m_vetoTagger;
  }else{
    localtagger = m_tagger;
  }

   // flat cut for out of range pTs
   if (pT>m_maxRangePt)
     pT = m_maxRangePt;

   if (localtagger.spline != nullptr && localtagger.constcut == nullptr) {
     pT = pT/1000.0;
     double maxsplinept = localtagger.spline->GetXmax();
     if (pT>maxsplinept){ pT = maxsplinept; }
     cutval = localtagger.spline->Eval(pT);
   }

   else if (localtagger.constcut != nullptr && localtagger.spline == nullptr) {
     cutval = localtagger.constcut[0](0);
   }
   else{
    ATH_MSG_ERROR( "Bad cut configuration!" );
    return CorrectionCode::Error;
   }


   ATH_MSG_VERBOSE( "Cut value " << cutval );

   return CorrectionCode::Ok;
}
