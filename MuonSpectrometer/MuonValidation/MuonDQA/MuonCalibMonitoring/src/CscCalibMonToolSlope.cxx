/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "CscCalibMonToolSlope.h"



#include "CscCalibData/CscCalibReportContainer.h"
#include "CscCalibData/CscCalibReportSlope.h"

#include "TProfile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TF1.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>

CscCalibMonToolSlope::CscCalibMonToolSlope(const std::string & type, const std::string & name, 
    const IInterface* parent) : 
  CscCalibMonToolBase(type, name, parent){
  declareProperty("MaxSlopeDiff",m_slopeMaxDiff=0.5);
  declareProperty("MaxIntercept",m_interceptMax = 5.0);
  declareProperty("MaxChi2_NDF",m_chi2Max = 100);
  declareProperty("MaxFracDev",m_fracDevMax = .5);
  declareProperty("MaxPeaktDiff",m_peaktMaxDiff= 10);
  m_calibResultKey = "CscCalibResultSlope";
  declareProperty("calGraphKey", m_histKey = "cscSlopeCalibReport");
  declareProperty("subDir", m_subDir ="");
  declareProperty("DeadPulserLevelCutoff", m_deadPulserLevelCutoff = 50);   
  declareProperty("DeadADCCutoff", m_deadADCCutoff =100); //Cutoff after ped subtraction
  declareProperty("DoNeighborRatio", m_doNeighborRatios = true);
  declareProperty("HistAttenLevels", m_histAttenLevels = false , "Make histograms for all attenuation levels. Potentially requires a LOT of RAM.");


  /*##From CscCalibMonToolBase.cxx, for reference. Do not uncomment!
    declareProperty("MakeHashValueHist",m_makeHashHist=true);
    declareProperty("MakeLayerValueHists",m_makeLayerHists=false);
    declareProperty("DetailedResultHashIds",m_detailedHashIds);
   */

}

StatusCode CscCalibMonToolSlope::finalize()
{
  ATH_MSG_DEBUG( "Slope Finalizing " );
  delete m_slopeNewColl;
  delete m_slopeOldColl;
  delete m_slopeDiffColl;
  delete m_interceptColl;
  delete m_chi2Coll;
  delete m_deadNewColl;
  delete m_deadDiffColl;
  delete m_peaktNewColl;
  delete m_peaktOldColl;
  delete m_peaktDiffColl;

  return CscCalibMonToolBase::finalize();
}

StatusCode CscCalibMonToolSlope::initialize() 
{
  m_onlyExpectPrecisionHashIds = true; 
  ATH_CHECK(CscCalibMonToolBase::initialize());

  ATH_MSG_DEBUG( "Expected chamber layer is " << m_expectedChamberLayer );

  for(unsigned int hash = 0 ; hash <= m_maxHashId; hash++)
  {
    m_fracProfs.push_back(nullptr);
  }
  m_generic_path_csccalibmonitoring = "MUON_CSC_PULSER";

  return StatusCode::SUCCESS;
}


StatusCode CscCalibMonToolSlope::bookHistograms()
{
  ATH_CHECK( CscCalibMonToolBase::bookHistograms() );
  ATH_MSG_DEBUG( "CscCalibMonToolSlope : in bookHistograms()"  );

  if (newRunFlag())
  {
    std::string name,title,xaxis,yaxis;
    int highbound,lowbound,nbins;

    std::string geoPath = getGeoPath();
    std::string path = getFullPath(geoPath, "Misc", "");
    MonGroup monGroup(this, path, run, ATTRIB_MANAGED );

    //num bad histograms
    name = "h_csc_calib_numSignificant";
    title = "Number of significant results.";
    xaxis = "Catagory";
    yaxis = "Num channels with bad value.";
    lowbound = 1;
    highbound = 8;
    m_h_numBad = new TH1I(name.c_str(),title.c_str(),highbound-lowbound,lowbound,highbound);
    m_h_numBad->GetYaxis()->SetTitle(yaxis.c_str());     
    m_h_numBad->GetXaxis()->SetTitle(xaxis.c_str());
    m_h_numBad->GetXaxis()->SetBinLabel(m_slopeBadBin,"slope diff");
    m_h_numBad->GetXaxis()->SetBinLabel(m_interceptBadBin,"intercept diff");
    m_h_numBad->GetXaxis()->SetBinLabel(m_chi2BadBin,"chi^2/ndf ");
    m_h_numBad->GetXaxis()->SetBinLabel(m_peaktBadBin,"Peaking Time Diff");
    m_h_numBad->GetXaxis()->SetBinLabel(m_fracBadBin,"Fractional Deviation");
    m_h_numBad->GetXaxis()->SetBinLabel(m_deadBadBin,"New (or fixed) Dead Channel");
    m_h_numBad->GetXaxis()->SetBinLabel(m_missingBadBin,"Num chans w/ no data");
    m_h_numBad->SetFillColor(m_histCol);
    monGroup.regHist(m_h_numBad).ignore();

    //--overview histograms-----------------------------------------------------------------
    /*
       name = "h_csc_calib_peaktCompareOverview";


       title = "Differences between measured peaking time and values stored in condtions database";
       xaxis = "Difference (ns)";
       yaxis = "";
       lowbound = -10;
       highbound = 10;
       nbins = 100;
       m_h_peaktCompareOverview = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
       m_h_peaktCompareOverview->GetXaxis()->SetTitle(xaxis.c_str());
       m_h_peaktCompareOverview->GetYaxis()->SetTitle(yaxis.c_str());     
       m_h_peaktCompareOverview->SetFillColor(m_histCol);
       sc = monGroup.regHist(m_h_peaktCompareOverview);

       name = "h_csc_calib_slopeCompareOverview";
       title = "Differences between measured slope  and values stored in condtions database";
       xaxis = "Difference (fC/ADC)";
       yaxis = "";
       lowbound = -10;
       highbound = 10;
       nbins = 100;
       m_h_slopeCompareOverview = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
       m_h_slopeCompareOverview->GetXaxis()->SetTitle(xaxis.c_str());
       m_h_slopeCompareOverview->GetYaxis()->SetTitle(yaxis.c_str());     
       m_h_slopeCompareOverview->SetFillColor(m_histCol);
       sc = monGroup.regHist(m_h_slopeCompareOverview);

       name = "h_csc_calib_interceptOverview";
       title = "Intercepts";
       xaxis = "Intercept (ADC count)";
       yaxis = "";
       lowbound = 0;
       highbound = 10;
       nbins = 100;
       m_h_interceptOverview = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
       m_h_interceptOverview->GetXaxis()->SetTitle(xaxis.c_str());
       m_h_interceptOverview->GetYaxis()->SetTitle(yaxis.c_str());     
       m_h_interceptOverview->SetFillColor(m_histCol);
       sc = monGroup.regHist(m_h_interceptOverview);

       name = "h_csc_calib_chi2Overview";
       title = "chi^2/ndf for slope fits";
       xaxis = "chi^2/ndf";
       yaxis = "";
       lowbound = 0;
       highbound = 10;
       nbins = 100;
       m_h_chi2Overview = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
       m_h_chi2Overview->GetXaxis()->SetTitle(xaxis.c_str());
       m_h_chi2Overview->GetYaxis()->SetTitle(yaxis.c_str());     
       m_h_chi2Overview->SetFillColor(m_histCol);
       sc = monGroup.regHist(m_h_chi2Overview);
     */
    name = "h_csc_calib_deadOverview";
    title = "Number of dead channels";
    lowbound = 1;
    highbound = 5;
    nbins = 4;
    m_h_deadOverview = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
    m_h_deadOverview->GetXaxis()->SetBinLabel(m_totalLiveBin,"Total Live channels");
    m_h_deadOverview->GetXaxis()->SetBinLabel(m_totalDeadBin,"Total Dead channels");
    m_h_deadOverview->GetXaxis()->SetBinLabel(m_newLiveBin,"New Live Channels");
    m_h_deadOverview->GetXaxis()->SetBinLabel(m_newDeadBin,"New Dead Channels");
    m_h_deadOverview->SetFillColor(m_histColAlert);
    monGroup.regHist(m_h_deadOverview).ignore();

    name = "h_csc_calib_slopeMissingChans";
    title = "Number of dead channels";
    xaxis = "Channel (Hash ID)";
    yaxis = "";
    lowbound = 0;
    highbound = m_maxHashId+1;
    nbins = m_maxHashId+1;
    m_h_slopeMissingChans = new TH1F(name.c_str(),title.c_str(),nbins,lowbound,highbound);
    m_h_slopeMissingChans->GetXaxis()->SetTitle(xaxis.c_str());
    m_h_slopeMissingChans->GetYaxis()->SetTitle(yaxis.c_str());     
    m_h_slopeMissingChans->SetFillColor(m_histColAlert);
    monGroup.regHist(m_h_slopeMissingChans).ignore();


    std::string peaktDataName        = "peakt";
    std::string peaktDataTitle       = "Peaking Time";
    std::string peaktSubDir          = "Peakt";

    std::string slopeDataName        = "slope";
    std::string slopeDataTitle       = "Pulser Gain Slope";
    std::string slopeSubDir          = "Slope";
    
    std::string slopeRatioDataName        = "slopeRatio";
    std::string slopeRatioDataTitle       = "Ratio of N : N+1 Channel Slopes";
    std::string slopeRatioSubDir          = "SlopeRatio";

    std::string interceptDataName    = "intercept";  
    std::string interceptDataTitle    = "Intercept";
    std::string interceptSubDir      = "Intercept";

    std::string chi2DataName         = "chi2";
    std::string chi2DataTitle        = "Chi^2/ndf for gain slope fit";
    std::string chi2SubDir           = "Chi2";

    std::string deadDataName         = "dead";
    std::string deadDataTitle        = "Dead";
    std::string deadSubDir           = "Dead";

    std::string fitResDataName         = "fitRes";
    std::string fitResDataTitle        = "Fit Return Value";
    std::string fitResSubDir           = "FitResult";

    //Set naming parameters for histogram category names
    std::string newCatName       = "new";
    std::string newCatTitle      = "New";

    std::string oldCatName       = "old";
    std::string oldCatTitle      = "COOL";

    std::string diffCatName      = "diff";
    std::string diffCatTitle     = "Difference Between New and COOL";

    //axis info
    std::string peaktAxisLabel = "Peaking Time (ns)";
    std::string peaktDiffAxisLabel = "Peaking Time Difference (ns)";
    int peaktNumBins = 100;
    float peaktLowBound = 0;
    float peaktHighBound = 100;

    std::string slopeAxisLabel = "Gain (fC/ADC)";
    std::string slopeDiffAxisLabel = "Gain Difference (fC/ADC)";
    int slopeNumBins =300;
    float slopeLowBound = 0;
    float slopeHighBound = 5;
    
    std::string slopeRatioAxisLabel = "Ratio of N/(N+1) channel";
    int slopeRatioNumBins = 500;
    float slopeRatioLowBound = 0;
    float slopeRatioHighBound = 5;

    std::string interceptAxisLabel = "Intercept (ADC counts)";
    int interceptNumBins = 200;
    float interceptLowBound = -100;
    float interceptHighBound = 100;

    std::string chi2AxisLabel = "Chi^2/ndf";
    int chi2NumBins = 1000;
    float chi2LowBound = 0;
    float chi2HighBound = 3000; 

    std::string deadAxisLabel = "Is Dead";
    int deadNumBins = 3;
    float deadLowBound = -1.5;
    float deadHighBound = 1.5;

    std::string fitResAxisLabel = "Fit Return Value";
    int fitResNumBins =10000;
    float fitResLowBound = 0;
    float fitResHighBound = 10000;

    //Initialize histogram collections
    m_peaktNewColl = new HistCollection(m_maxHashId +1, m_maxHashId +1);
    m_peaktNewColl->ignoreY = false;
    m_peaktOldColl = new HistCollection(m_maxHashId +1);
    m_peaktOldColl->ignoreY = false;
    m_peaktDiffColl = new HistCollection(m_maxHashId +1);
    m_peaktDiffColl->ignoreY = false;


    if(m_doNeighborRatios){
      m_slopeRatioColl = new HistCollection(m_maxHashId +1);
      m_slopeRatioColl->ignoreY = true;
    }

    m_slopeNewColl = new HistCollection(m_maxHashId +1, m_maxHashId +1);
    m_slopeNewColl->ignoreY = true;
    m_slopeOldColl = new HistCollection(m_maxHashId +1);
    m_slopeOldColl->ignoreY = true;
    m_slopeDiffColl = new HistCollection(m_maxHashId +1);
    m_slopeDiffColl->ignoreY = true;

    m_interceptColl = new HistCollection(m_maxHashId +1, m_maxHashId +1);
    m_interceptColl->ignoreY = true;

    m_chi2Coll = new HistCollection(m_maxHashId +1);
    m_chi2Coll->ignoreY = true;


    m_deadNewColl = new HistCollection(m_maxHashId +1, m_maxHashId +1);
    m_deadNewColl->ignoreY = false;
    m_deadDiffColl = new HistCollection(m_maxHashId +1,m_maxHashId +1);
    m_deadDiffColl->ignoreY = false;

    m_fitResColl = new HistCollection(m_maxHashId +1);
    m_fitResColl->ignoreY = true;

    //Set subdirectory in monitoring root file

    //initialize, name, and book histograms in histogram collections using the bookHistCollection()
    //function

    ATH_CHECK(  bookHistCollection(m_peaktNewColl, peaktDataName, peaktDataTitle, newCatName, 
                                   newCatTitle, peaktAxisLabel, peaktNumBins, peaktLowBound, peaktHighBound, peaktSubDir) );

    ATH_MSG_DEBUG( "Registering peaktOldColl"  );
    ATH_CHECK( bookHistCollection(m_peaktOldColl, peaktDataName, peaktDataTitle, oldCatName, oldCatTitle,
                                  peaktAxisLabel, peaktNumBins, peaktLowBound, peaktHighBound, peaktSubDir) );

    ATH_MSG_DEBUG( "Registering peaktDiffColl"  );
    ATH_CHECK( bookHistCollection(m_peaktDiffColl, peaktDataName, peaktDataTitle, diffCatName, diffCatTitle,
                                  peaktDiffAxisLabel, peaktNumBins, peaktLowBound, peaktHighBound, peaktSubDir) );

    ATH_MSG_DEBUG( "Registering slopeRatioColl"  );
    ATH_CHECK( bookHistCollection(m_slopeNewColl, slopeDataName, slopeDataTitle, newCatName, 
                                  newCatTitle, slopeAxisLabel, slopeNumBins, slopeLowBound, slopeHighBound, slopeSubDir) );

    ATH_MSG_DEBUG( "Registering slopeNewColl"  );
    ATH_CHECK( bookHistCollection(m_slopeRatioColl, slopeRatioDataName, slopeRatioDataTitle, "", 
                                  "", slopeRatioAxisLabel, slopeRatioNumBins, slopeRatioLowBound, slopeRatioHighBound, slopeRatioSubDir) );

    ATH_MSG_DEBUG( "Registering slopeOldColl"  );
    ATH_CHECK( bookHistCollection(m_slopeOldColl, slopeDataName, slopeDataTitle, oldCatName, oldCatTitle,
                                  slopeAxisLabel, slopeNumBins, slopeLowBound, slopeHighBound, slopeSubDir) );

    ATH_MSG_DEBUG( "Registering slopeDiffColl"  );
    ATH_CHECK( bookHistCollection(m_slopeDiffColl, slopeDataName, slopeDataTitle, diffCatName, diffCatTitle,
                                  slopeDiffAxisLabel, slopeNumBins, slopeLowBound, slopeHighBound, slopeSubDir) );

    ATH_MSG_DEBUG( "Registering " << interceptDataTitle  );
    ATH_CHECK( bookHistCollection(m_interceptColl, interceptDataName, interceptDataTitle, "",
                                  "", interceptAxisLabel, interceptNumBins, interceptLowBound, interceptHighBound, 
                                  interceptSubDir) );

    ATH_CHECK( bookHistCollection(m_chi2Coll, chi2DataName, chi2DataTitle, "", "",
                                  chi2AxisLabel, chi2NumBins, chi2LowBound, chi2HighBound, chi2SubDir) );

    ATH_CHECK(  bookHistCollection(m_deadNewColl, deadDataName, deadDataTitle, newCatName, newCatTitle,
                                   deadAxisLabel, deadNumBins, deadLowBound, deadHighBound, deadSubDir) );
    
    ATH_CHECK( bookHistCollection(m_fitResColl, fitResDataName, fitResDataTitle, "", "",
                                  fitResAxisLabel, fitResNumBins, fitResLowBound, fitResHighBound, fitResSubDir) );

  }
  return StatusCode::SUCCESS;

}

/*StatusCode CscCalibMonToolSlope::fillHistograms()
  {
  StatusCode sc = StatusCode::SUCCESS;

  if (m_debuglevel) 
  m_log << MSG::DEBUG << "CscCalibMonToolSlope :: in fillHistograms()" << endmsg;

  return sc;
  }*/

//--handleParameter: Processes a vector of parameter values by filling the appropriate histograms
StatusCode CscCalibMonToolSlope::handleParameter(const CscCalibResultCollection* parVals)
{
  ATH_MSG_DEBUG( "CscCalibMonToolSlope : in handleParameter()"  );

  //The whole point of this funciton is to pass the correct histograms and setup info 
  //to CsccalibMonToolBase::procParameter. To organize this, we store the setup info into
  //these structs: 
  ProcSetupInfo simpleSet;

  //--setup for this parameter
  std::string parName = parVals->parName();
  /*if(parName == "peakt")
    {
    m_log << MSG::INFO << "Evaluating peaking times" << endmsg;
    simpleSet.dbName = parVals->parName();
    simpleSet.badHist = m_h_numBad;
    simpleSet.badBin = m_peaktBadBin;
    simpleSet.maxDiff = m_peaktMaxDiff;
    simpleSet.overDiff = m_h_peaktCompareOverview;
    simpleSet.overChi2 = m_h_chi2Overview;
    simpleSet.chi2BadBin = m_chi2BadBin;
    simpleSet.doChi2 = false;
    simpleSet.vals = &m_peaktNewColl->data;
    simpleSet.errors = &m_peaktNewColl->errors;
    simpleSet.oldVals = &m_peaktOldColl->data;
    simpleSet.diffs = &m_peaktDiffColl->data;
    simpleSet.missingBadBin = m_missingBadBin;
    simpleSet.missingChans = m_h_slopeMissingChans;
    simpleSet.expectedChannels = m_expectedHashIdsPrec;
    }
   */
  if(parName == "pslope")
  {
    ATH_MSG_INFO( "Evaluating slopes"  );
    simpleSet.dbName = parVals->parName();
    simpleSet.badHist = m_h_numBad;
    simpleSet.badBin = m_slopeBadBin;
    simpleSet.maxDiff = m_slopeMaxDiff;
    simpleSet.overDiff = m_h_slopeCompareOverview;
    simpleSet.overChi2 = m_h_chi2Overview;
    simpleSet.chi2BadBin = m_chi2BadBin;
    simpleSet.chi2Max = m_chi2Max;
    simpleSet.doChi2 = true;
    simpleSet.vals = &m_slopeNewColl->data;
    simpleSet.errors = &m_slopeNewColl->errors;
    simpleSet.oldVals = &m_slopeOldColl->data;
    simpleSet.diffs = &m_slopeDiffColl->data;
    simpleSet.chi2s = &m_chi2Coll->data;
    simpleSet.missingBadBin = m_missingBadBin;
    simpleSet.missingChans = m_h_slopeMissingChans;
    simpleSet.expectedChannels = m_expectedHashIdsPrec;
  }
  else if (parName == "pinter")
  {
    ATH_MSG_INFO( "Evaluating intercepts"  );
    simpleSet.expectedVal = 0;
    simpleSet.badHist = m_h_numBad;
    simpleSet.badBin = m_interceptBadBin;
    simpleSet.maxDiff = m_interceptMax;
    simpleSet.overDiff = m_h_interceptOverview;
    simpleSet.overChi2 = m_h_chi2Overview;
    simpleSet.doChi2 = false;
    simpleSet.vals = &m_interceptColl->data;
    simpleSet.errors = &m_interceptColl->errors;
    simpleSet.oldVals = nullptr ;
    simpleSet.diffs = nullptr ;
    simpleSet.missingBadBin = m_missingBadBin;
    simpleSet.missingChans = m_h_slopeMissingChans;
    simpleSet.expectedChannels = m_expectedHashIdsPrec;
  }
  else
  {
    ATH_MSG_INFO( "CscCalibMonToolSlope : Did not recognize parameter name " 
                  << parName << ". This is usually ok."  );
    return StatusCode::FAILURE;
  }

  //Process parameter by filling histograms in simpleSet and allIdsSet structures
  ATH_CHECK( procParameter(parVals,&simpleSet) );


  if(parName == "peakt")
  {
    ATH_MSG_INFO( "Generating peaking time histograms"  );

    ATH_CHECK( copyDataToHists(m_peaktNewColl) &
               copyDataToHists(m_peaktOldColl) &
               copyDataToHists(m_peaktDiffColl) );
  }
  if(parName == "pslope") 
  {

    ATH_MSG_INFO( "Generating slope histograms"  );
    
    if(m_doNeighborRatios){
      genNeighborRatios(m_slopeNewColl->data, m_slopeRatioColl->data);
      ATH_CHECK( copyDataToHists(m_slopeRatioColl) );
    }

    ATH_CHECK( copyDataToHists(m_slopeNewColl) &
               copyDataToHists(m_slopeOldColl) &
               copyDataToHists(m_slopeDiffColl) &
               copyDataToHists(m_chi2Coll) );
  }
  if(parName == "pinter")
  {
    ATH_MSG_INFO( "Generating intercept histograms"  );
    ATH_CHECK( copyDataToHists(m_interceptColl) );
  }
  return StatusCode::SUCCESS;
}

//--retrieveHistos() will retrieve the slope amplitude histograms for the channels
//requested by the user in m_detailedHashIds. 
StatusCode CscCalibMonToolSlope::postProc()
{
  ATH_MSG_DEBUG( "CscCalibMonToolSlope : in retrieveHistos()"  );

  ATH_MSG_DEBUG( "Do all detailed is: " << m_doAllDetailed  );



  IdContext chanContext = m_idHelperSvc->cscIdHelper().channel_context();

  //Get the slopeReport, checking for pointer errors along the way
  const DataHandle<CscCalibReportContainer> repCont;
  if (!evtStore()->retrieve(repCont, m_histKey).isSuccess())
  {
    ATH_MSG_ERROR( " Cannot retrieve object from storegate with key "
                   << m_histKey <<  " aborting retrieving hists "  );
    return StatusCode::RECOVERABLE;
  }
  if(repCont->size() != 1)
  {
    ATH_MSG_ERROR( "Container with key " << m_histKey
                   << " does not have a size of one. Do not know how to proceed, so aborting"
                   << " retrieving calibration histograms."  );
    return StatusCode::RECOVERABLE;
  }

  const CscCalibReportSlope * slopeReport =
    dynamic_cast<const CscCalibReportSlope *>(repCont->front());
  if(!slopeReport)
  {
    ATH_MSG_ERROR( "No report stored in the container with key " << m_histKey
                   << ". Aborting retrieving histograms."  );
    return StatusCode::RECOVERABLE;
  }
  if(slopeReport->getLabel() != "calGraphs")
  {
    ATH_MSG_ERROR( "Incorect object retrieved from container with key " << m_histKey
                   << ". Aborting hist retrieval"  );
    return StatusCode::RECOVERABLE;
  }
  
 //ampProfs 
  const std::map <int,TProfile*> * ampProfs = slopeReport->getAmpProfs();
  if(!ampProfs)
  {
    ATH_MSG_ERROR( "There are no amplitude profiles in the slope report! Can't find dead chans." );
    return StatusCode::RECOVERABLE;
  }

  if(m_histAttenLevels){
    for(const auto & [attenuationVal, pProfile] : *ampProfs){

      const float atten = attenuationVal * 0.5f;
      const std::string attenStr = std::to_string(atten);

      HistCollection * ampColl = new HistCollection(m_maxHashId +1, m_maxHashId +1);
      ampColl->ignoreY = true;

      //put into ampColls to make for easy deleting later
      m_ampColls.push_back(ampColl);


      std::string dataName = "amp_atten" + attenStr;
      std::string dataTitle = "ADC response at attenuation " + attenStr + " db" ; 
      std::string axisLabel = "ADC";
      unsigned int numBins = 300; 
      float lowBound = 0;
      float highBound = 3000;
      std::string subDir = "AmpAtten" + attenStr;

      bookHistCollection(ampColl, dataName, dataTitle, "", "", axisLabel, numBins, lowBound, highBound, subDir).ignore();
      for(unsigned int stripHash = 0; stripHash < m_maxHashId; stripHash++){
        ampColl->data[stripHash] =  pProfile->GetBinContent(stripHash +1);
      }

      ATH_CHECK( copyDataToHists(ampColl) );
    }
  }


  const std::vector<float> * fitResVec = slopeReport->getFitResults();
  m_fitResColl->data = *fitResVec; 
  ATH_CHECK( copyDataToHists(m_fitResColl) );

  //Generate fractional deviation histograms
  ATH_MSG_DEBUG( "About to generate fractional deviation graphs" );
  if(!makeFracGraphs(*slopeReport).isSuccess())
    ATH_MSG_WARNING( "Failed to generate fractional deviation graphs. Continuing anyway.." );


  ATH_MSG_DEBUG( "About to find dead channels"  );
  //Determine dead channels
  //sc = findDeadChannels(*slopeReport);
  //if(!sc.isSuccess())
  //  m_log << MSG::WARNING << "Failure while finding dead channels" << endmsg;


  //Put extra info for those channels indicated in m_detailedHashIds
  ATH_MSG_DEBUG( "Picking detailed graphs to output to root file" );
  if(m_numBad >0 || m_maxDetailedChannels < 0 || m_doAllDetailed)
  { 
    const DataVector<TGraphErrors> * calGraphs 
      = slopeReport->getCalGraphs();
    if(!calGraphs)
    {
      ATH_MSG_ERROR( "No calGraph stored inside object with key " << m_histKey
                     << ". Aborting hist retrieval."  );
      return StatusCode::RECOVERABLE;
    }
    else
      ATH_MSG_INFO( "Got calGraphs"  );

    const DataVector<TH1I> * bitHists = slopeReport->getBitHists();
    if(!bitHists)
      ATH_MSG_INFO( "No bit histogram vector found from calibration. "
                    << " Won't be in monitoring output file. "  );
    else
      ATH_MSG_INFO( "Got bitHists"  );


    //These are the channels we will get detailed forr.
    for(unsigned int idItr = 0; idItr < m_maxHashId; idItr++)
    {
      ATH_MSG_ERROR("Calgraph Address: " <<  (*calGraphs)[idItr]);
      if(m_expectedHashIdsPrec.count(idItr) && (m_detailedHashIds[idItr] || (m_doAllDetailed) ) )
      {
        ATH_MSG_ERROR("Doing detailed plots of hash " << idItr);

        Identifier chanId;
        m_idHelperSvc->cscIdHelper().get_id(IdentifierHash(idItr), chanId, &chanContext);
        int stationSize = m_idHelperSvc->cscIdHelper().stationName(chanId);
        int stationEta = m_idHelperSvc->cscIdHelper().stationEta(chanId);
        int stationPhi = m_idHelperSvc->cscIdHelper().stationPhi(chanId);
        int wireLayer = m_idHelperSvc->cscIdHelper().wireLayer(chanId);
        int measuresPhi = m_idHelperSvc->cscIdHelper().measuresPhi(chanId);
        int strip = m_idHelperSvc->cscIdHelper().strip(chanId);
        int sector = getSector(stationPhi, stationSize);

        std::string geoPath = getGeoPath(stationEta, sector, wireLayer, measuresPhi);
        std::string calGraphPath = getFullPath(geoPath, "CalGraphs","");
        std::string fracPath = getFullPath(geoPath,"FracProf","");
        std::string bitHistPath = getFullPath(geoPath, "BitHists", "");

        //Record calgraph
        TGraphErrors  * sourceGraph = 
          const_cast<TGraphErrors*>((*calGraphs)[idItr]);
        if(!sourceGraph)
        {
          ATH_MSG_ERROR( "The requested calgraph for hash "
                         << idItr << " doesn't exist in CscCalibReport object!"  );
        }
        else {
          std::stringstream  name;
          name << "calfit" 
            << "_EC" << getEndCap(stationEta)
            << "_sector_" << sector 
            << "_layer_" << wireLayer
            << "_" << (measuresPhi ? "trans" : "prec")
            << "_strip_" 
            << std::setfill ('0') <<  std::setw (measuresPhi ? 2 : 3) 
            << strip;

          sourceGraph->SetName(name.str().c_str());
          ATH_MSG_DEBUG( "CalGraph axis title: " << sourceGraph->GetXaxis()->GetTitle()  );

          ATH_CHECK( regGraph(sourceGraph, calGraphPath, run, ATTRIB_MANAGED) );
        }

        //record fracDev graph
        TProfile * sourceProf = m_fracProfs[idItr];
        if(!sourceProf)
        {
          ATH_MSG_ERROR( "There is no fractional profile available for hash " 
                         << idItr << ". Quitting retrieveHistos()."  );
        }
        else{
          std::stringstream  fracName;
          fracName << "frac_"
            << "_EC" << getEndCap(stationEta)
            << "_sector_" << sector 
            << "_layer_" << wireLayer
            << "_" << (measuresPhi ? "trans" : "prec")
            << "_strip_" 
            << std::setfill ('0') <<  std::setw (measuresPhi ? 2 : 3) 
            << strip;

          sourceProf->SetName(fracName.str().c_str());

          ATH_CHECK( regHist(sourceProf, fracPath, run, ATTRIB_MANAGED) );
        }

        //Bit map histograms
        //copy source histogram into new histogram, and store
        if(bitHists)
        {
          TH1I *  bitHist = const_cast<TH1I*>((*bitHists)[idItr]);
          if(!bitHist)
          {
            ATH_MSG_ERROR( "There is no bit histogram with hashId "
                           << idItr << " Quiting out of detailed histogram loop."  );
          }
          else {

            std::stringstream  name2;
            name2 << "h_bitMap"
              << "_EC" << getEndCap(stationEta)
              << "_sector_" << sector 
              << "_layer_" << wireLayer
              << "_" << (measuresPhi ? "trans" : "prec")
              << "_strip_" 
              << std::setfill ('0') <<  std::setw (measuresPhi ? 2 : 3) 
              << strip;
            // TH1I * newHist2 = (TH1I*)bitHist->Clone(name2.str().c_str());
            bitHist->SetName(name2.str().c_str());
            bitHist->SetFillColor((m_detailedHashIds[idItr] ? m_histColAlert : m_histCol));

            //ATH_CHECK( regHist(newHist2, chanPath, run, ATTRIB_MANAGED) );
            ATH_CHECK( regHist(bitHist, bitHistPath, run, ATTRIB_MANAGED) );
          }
        }//end if bithists*/

      }//end if (m_detailedHashIds)
      else
        ATH_MSG_ERROR("Skipping hash " << idItr << " " <<  m_expectedHashIdsPrec.count(idItr) << " " << m_doAllDetailed);
    }//end idItr loop
  }//if numBad > 0
  return StatusCode::SUCCESS;    
}// end retrieveHistos

//Generate fractional Deviation graphs
StatusCode CscCalibMonToolSlope::makeFracGraphs(const CscCalibReportSlope & slopeReport)
{
  ATH_MSG_DEBUG( "CscCalibMonToolSlope : in makeFracGraphs()"  );
  const DataVector<TGraphErrors> * calGraphs = slopeReport.getCalGraphs();
  if(!calGraphs)
  {
    ATH_MSG_ERROR( "No calGraphs in slopeReport. Not going to make fractional deviation"
                   << " plots."  );
    return StatusCode::RECOVERABLE;
  }

  //Loop through all channels in geometry:
  const std::vector <Identifier> & ids = m_idHelperSvc->cscIdHelper().idVector();
  for(const auto & thisChamberId:ids)
  {
    ATH_MSG_VERBOSE( "in Chamber loop "  );
    unsigned int stationSize = m_idHelperSvc->cscIdHelper().stationName(thisChamberId); //51 = large, 50 = small
    unsigned int stationPhi = m_idHelperSvc->cscIdHelper().stationPhi(thisChamberId);
    int stationEta = m_idHelperSvc->cscIdHelper().stationEta(thisChamberId);
    unsigned int sector = getSector(stationPhi,stationSize);

    std::vector <Identifier> stripVect;
    m_idHelperSvc->cscIdHelper().idChannels(thisChamberId,stripVect);
    for(const auto & thisStrip : stripVect)
    {
      ATH_MSG_VERBOSE( "in strip loop "  );
      IdentifierHash stripHash;
      m_idHelperSvc->cscIdHelper().get_channel_hash(thisStrip,stripHash);
      if(!m_expectedHashIdsPrec.count((int)stripHash)){
        ATH_MSG_VERBOSE( "Skipping hash"  << (int)stripHash  );
        continue;
      }
      ATH_MSG_VERBOSE( "strip hash " << stripHash  );

      const TGraphErrors * graph = (*calGraphs)[stripHash];
      if(!graph)
      {
        ATH_MSG_VERBOSE( "SKipping graph"  );
        continue;
      }
      
      TF1 * func = graph->GetFunction("simpleFunc");
      if(!func){
        ATH_MSG_DEBUG( "Could not retrieve function, skipping this graph"  );
        continue;
      }

      ATH_MSG_VERBOSE( "Address is " << graph );
      ATH_MSG_VERBOSE( "name is "<< graph->GetName()  );
      ATH_MSG_VERBOSE( "Getting n "  );

      int nPoints = graph->GetN();
      ATH_MSG_VERBOSE( "Got n "  );
      //Get identification information
      //Note, we don't ask for measuresPhi because there should be no
      //TGraphs with Y anyways.
      ATH_MSG_VERBOSE( "getting id info " );
      unsigned int layer = m_idHelperSvc->cscIdHelper().wireLayer(thisStrip);
      unsigned int strip = m_idHelperSvc->cscIdHelper().strip(thisStrip);
      ATH_MSG_VERBOSE( "Got strip and layer"  );
      //initialize fractional deviation profile
      ATH_MSG_VERBOSE( "initializing profile "  );

      std::stringstream  nameStream;
      nameStream.setf(std::ios::right,std::ios::adjustfield);
      nameStream << "dev_" 
        << "X"   //orientation
        << "_eta_" << stationEta
        << "_sector_" << std::setw (2) << std::setfill ('0') << sector
        << "_layer_" << layer
        << "_strip_" << strip;
      std::stringstream  titleStream;
      titleStream << "Fractional Deviation of Measured ADC From Fit ADC for Precision Direction"
        << ", Sector " << sector
        << ", Eta "     << stationEta
        << ", Layer " << layer
        << ", Strip " << strip;
      TProfile * fracProf 
        = new TProfile(nameStream.str().c_str(), titleStream.str().c_str(), nPoints, 
            1, nPoints +1);
      fracProf->GetXaxis()->SetTitle("Pulser Attenuator Settings in dB (not auscaled axis) ");
      fracProf->GetYaxis()->SetTitle("ADC response"); 

      ATH_MSG_VERBOSE( "getting parameters "  );
      float intercept = func->GetParameter(0);
      float slope = func->GetParameter(1);
      ATH_MSG_VERBOSE( "got them "  );

      //Loop through points in graph, and compare the adc value and charge
      //to that expected from the fit parameters from the TF1 object
      double adc, db, calcAdc;
      bool isBad = false;
      for(int itr = 0; itr < nPoints; itr++)
      {
        ATH_MSG_VERBOSE( "in point loop on point " << itr );

        //Find Fractional Deviation
        graph->GetPoint(itr,db,adc);
        calcAdc = intercept + slope*std::pow(10,db/-20);
        double frac = (adc - calcAdc)/( calcAdc - intercept);

        //Test to see if this is a bad result
        ATH_MSG_VERBOSE( "bad result test "  );
        if(frac > m_fracDevMax)
        {
          m_h_numBad->Fill(m_fracBadBin);
          isBad = true;
        }

        //test for nan, and put into fracProf 
        if (frac==frac)
        {
          ATH_MSG_VERBOSE( "filling fracProf "  );
          fracProf->Fill(itr +1, frac);     
        }

        //Label bin with db amount
        ATH_MSG_VERBOSE( "labeling bin "  );
        std::stringstream  binLabel;
        binLabel << db;    
        fracProf->GetXaxis()->SetBinLabel(itr+1, binLabel.str().c_str());
      }
      ATH_MSG_VERBOSE( "storing fracProf "  );
      //Store fraction graph for later
      if(stripHash > m_maxHashId)
      {
        ATH_MSG_ERROR( "Tried to assign fracProf with stripHash " << stripHash
                       << " when max stripHash is " << m_maxHashId  );
        return StatusCode::RECOVERABLE;
      }

      m_fracProfs[stripHash] = fracProf; 
      //Record if this was a bad channel
      if(isBad && ((int)m_numBad <= m_maxDetailedChannels) )
      {
        if(!m_detailedHashIds[stripHash])
        {
          m_numBad++;
          m_detailedHashIds[stripHash] = true;
        }
      }

    }//end loop over strips
  }//end loop over chambers
  ATH_MSG_DEBUG( "Finished generating frac graphs"  );

  return StatusCode::SUCCESS;
}//end makeFracGraphs

StatusCode CscCalibMonToolSlope::findDeadChannels(const CscCalibReportSlope & slopeReport)
{//This function has grown a bit unwieldly, and duplicates info in sets and arrays. Can this be merged?

  SG::ReadCondHandle<CscCondDbData> readHandle{m_readKey};
  const CscCondDbData* readCdo{*readHandle};

  //****Find Dead channels
  IdContext channelContext = m_idHelperSvc->cscIdHelper().channel_context(); 

  std::set <int> newDead, newUndead;

  const std::set <int> * pulsedChambers = slopeReport.getPulsedChambers();
  if(!pulsedChambers)
  {
    ATH_MSG_ERROR( "No pulsed chambers stored in slopeReport! Skipping dead channel collecting!"  );
    return StatusCode::RECOVERABLE;
  }

  const std::map <int,TProfile*> * ampProfs = slopeReport.getAmpProfs();
  if(!ampProfs)
  {
    ATH_MSG_ERROR( "There are no amplitude profiles in the slope report! Can't find dead chans." );
    return StatusCode::RECOVERABLE;
  }

  std::map <int,TProfile*>::const_iterator profItr = ampProfs->begin();

  int pulserLevel = profItr->first;
  ATH_MSG_INFO( "Looking for dead channels. Lowest attenuation level is " 
                << pulserLevel  );
  ATH_MSG_INFO( "Dead channel cutoff is: "<< m_deadADCCutoff   );

  if(pulserLevel < m_deadPulserLevelCutoff)
  {//Pulser level is low enough (pulse voltage is high enough) to test for dead channels

    const TProfile * ampProf = profItr->second;
    if(!ampProf)
    {
      ATH_MSG_ERROR( "There is no profile for this attenuation level! Skipping dead channel finding!" );
      return StatusCode::RECOVERABLE;
    }

    ATH_MSG_DEBUG( "There were " << pulsedChambers->size()
                   << " chambers pulsed."  );

    //Prepare dead channel content
    std::vector<float> & currentDeadVals = m_deadNewColl->data;
    std::vector<float> & diffDeadVals = m_deadDiffColl->data;

    setArray(currentDeadVals,0,m_maxHashId+1);
    setArray(diffDeadVals,0,m_maxHashId+1);

    //Loop through each channel in TProfile (bin = hash + 1), and check value to see
    //if its dead
    float adc, ped;
    bool wasDead, isDead;
    int statusWord;
    for(unsigned int hashItr = 0; hashItr <= m_maxHashId ; hashItr++)
    {
      Identifier id;
      m_idHelperSvc->cscIdHelper().get_id(hashItr,id, &channelContext);
      int chamberLayer = m_idHelperSvc->cscIdHelper().chamberLayer(id);
      IdentifierHash chamberHash;
      m_idHelperSvc->cscIdHelper().get_module_hash(id, chamberHash);

      if(chamberLayer == 2 && pulsedChambers->count((int)chamberHash))
      {//This is a good chamber layer and it is a pulsed chamber

        ATH_CHECK(readCdo->readChannelStatus(hashItr, statusWord));
        ATH_CHECK(readCdo->readChannelPed   (hashItr, ped       ));
        wasDead = statusWord & 0x1;

        adc = ampProf->GetBinContent( hashItr + 1 );
        isDead = ( (adc-ped) < m_deadADCCutoff);

        if(isDead && !wasDead)
        {
          ATH_MSG_INFO( "NEW DEAD CHANNEL! Hash: " << hashItr << " ped: " << ped << " adc " << adc << " diff: " << adc-ped  );
          newDead.insert(hashItr);
          m_h_deadOverview->Fill(m_newDeadBin);
          diffDeadVals[hashItr] = 1;
          m_detailedHashIds[hashItr] = true;
        }
        if(!isDead && wasDead)
        {
          ATH_MSG_INFO( "PREVIOUSLY DEAD CHANNEL NOW LIVE! Hash: " << hashItr  );
          newUndead.insert(hashItr);
          m_h_deadOverview->Fill(m_newLiveBin);
          diffDeadVals[hashItr] = -1;
          m_detailedHashIds[hashItr] = true;
        }

        if(isDead)
        {
          currentDeadVals[hashItr] = 1;
          m_h_deadOverview->Fill(m_totalDeadBin);
        }
        else
          m_h_deadOverview->Fill(m_totalLiveBin);

      }//end if chamberLayer

    }//End for(hashItr)

    //***Generate histograms
    ATH_MSG_INFO( "Generating dead histograms"  );
    ATH_CHECK( copyDataToHists(m_deadNewColl) );

    //ATH_CHECK( copyDataToHists(m_deadDiffColl) );
    ATH_MSG_INFO( "Finished generating dead histograms"  );

    //***Fill COOL input file if there are new dead channels******
    if(newDead.size() || newUndead.size())
    {
      ATH_MSG_INFO( "There are " << newDead.size() 
                    << " newly dead channels and " << newUndead.size() 
                    << " newly live channels" );
      std::ofstream  out("deadInfo.cal");
      out <<"00-00 " << newDead.size() + newUndead.size() << " dead_stat END_HEADER\n";

      for(const auto & thisDeadChannel:newDead)
      {
        Identifier id;
        m_idHelperSvc->cscIdHelper().get_id(thisDeadChannel,id, &channelContext);
        IdentifierHash chamHash;
        m_idHelperSvc->cscIdHelper().get_module_hash(id, chamHash);
        out << thisDeadChannel << " " << (int)chamHash << " " 
          << m_idHelperSvc->cscIdHelper().show_to_string(id, &channelContext) << " 1\n";
      }

      for(const auto undeadChannel : newUndead)
      {
        Identifier id;
        m_idHelperSvc->cscIdHelper().get_id(undeadChannel,id, &channelContext);
        IdentifierHash chamHash;
        m_idHelperSvc->cscIdHelper().get_module_hash(id, chamHash);
        out << undeadChannel << " " << (int)chamHash << " " 
          << m_idHelperSvc->cscIdHelper().show_to_string(id, &channelContext)
          << "0\n";
      } 
      out.close();
    }//end if newDead.size()

  }//end if(pulserLevel < cutoff)
  else 
  {
    ATH_MSG_ERROR( "Lowest pulser level isn't low enough to count as a dead channel test. Skipping." );
    return StatusCode::RECOVERABLE;
  }


  return StatusCode::SUCCESS;
}

//I changed the array to a vector. At some point I should replace code like this to make more sense for vector
std::vector<float> & CscCalibMonToolSlope::setArray(std::vector<float>  & array, const float &value, const int &numEntries)
{
  for(int i = 0 ; i < numEntries; i++)
    array[i] = value;
  return array;
}
       
void CscCalibMonToolSlope::genNeighborRatios(const std::vector<float> & source, std::vector<float> & ratios) const{
  size_t nEntries = source.size();
  if(nEntries != ratios.size()){
    ATH_MSG_ERROR( " in genNeighborRatios, source (" << nEntries << ") and ratio (" << ratios.size() 
                   << ") vectors have different numbers of entries!"  );
    return;
  }

  //loop through every channel, taking the ratio between it and the next
  size_t nloops = nEntries -1; //Can't do this to the last channel
  for(size_t cnt=0; cnt < nloops; cnt++){
    ratios[cnt] = source[cnt]/source[cnt+1];
  }
}

