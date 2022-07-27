/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
   @class ElectronChargeEfficiencyCorrectionTool
   @brief Apply the correction for the difference in electron charge mis-ID rates
   in MC/data

   @author Giulia Gonella <giulia.gonella@cern.ch>
   @date   September 2015
*/
// Include this class's header
#include "ElectronEfficiencyCorrection/ElectronChargeEfficiencyCorrectionTool.h"
#include "ElectronEfficiencyCorrection/ElectronEfficiencyHelpers.h"
// xAOD includes
#include "PathResolver/PathResolver.h"
#include "xAODEgamma/Electron.h"
#include "xAODEventInfo/EventInfo.h"
#include "PATInterfaces/SystematicRegistry.h"

// ROOT includes
#include "TFile.h"

// STL includes
#include <cstdlib> /* atoi */

// =============================================================================
// Standard constructor
// =============================================================================
CP::ElectronChargeEfficiencyCorrectionTool::
  ElectronChargeEfficiencyCorrectionTool(const std::string& name)
  : AsgTool(name)
  , m_dataTypeOverwrite(-1)
  , m_eventInfoCollectionName("EventInfo")
  , m_SF_SS()
  , m_SF_OS()
  , m_RunNumbers()
  , m_useRandomRunNumber(true)
  , m_defaultRandomRunNumber(999999)
  , m_filename("")
  , m_workingPoint("")
  , m_eta_lowlimit(0.0)
  , m_eta_uplimit(0.0)
  , m_pt_lowlimit(0.0)
  , m_pt_uplimit(0.0)
  , m_gevmev(0.0)
  , m_filtered_sys_sets()
  , m_mySysConf()
  , m_affectingSys()
  , m_appliedSystematics(nullptr)
  , m_sf_decoration_name("chargeIDEffiSF")
  , m_sfDec(nullptr)
{
  // Declare the needed properties
  declareProperty("CorrectionFileName",
                  m_filename,
                  "Name of the file with charge flipping rates");
  declareProperty(
    "WorkingPoint", m_workingPoint, "Name of working point folder in the file");
  declareProperty("ScaleFactorDecorationName", m_sf_decoration_name);
  declareProperty("ForceDataType",
                  m_dataTypeOverwrite,
                  "Force the DataType of the electron to specified value (to "
                  "circumvent problem of incorrect DataType for forward "
                  "electrons in some old releases)");
  declareProperty("EventInfoCollectionName",
                  m_eventInfoCollectionName,
                  "The EventInfo Collection Name");
  declareProperty("UseRandomRunNumber", m_useRandomRunNumber);
  declareProperty("DefaultRandomRunNumber", m_defaultRandomRunNumber);
}

// =============================================================================
// Standard destructor
// =============================================================================
CP::ElectronChargeEfficiencyCorrectionTool::
  ~ElectronChargeEfficiencyCorrectionTool()
{
  if (m_sfDec)
    delete m_sfDec;
}

// =============================================================================
// Athena initialize method
// =============================================================================
StatusCode
CP::ElectronChargeEfficiencyCorrectionTool::initialize()
{
  ATH_MSG_DEBUG("initializing");

  // initialize the random number generator (used in case of charge flip
  // approach)
  // m_Rndm = new TRandom3(1);

  if (m_sfDec)
    delete m_sfDec;
  m_sfDec = new SG::AuxElement::Decorator<float>(m_sf_decoration_name); // xxxx

  // Resolve the path to the input file for the charge flip rates
  const std::string rootfilename = PathResolverFindCalibFile(m_filename);
  if (m_filename.empty()) {
    ATH_MSG_ERROR(" PathResolver was not able to find the file ... aborting");
    return StatusCode::FAILURE;
  }

  // Getting the root file and histograms
  TFile* rootFile = TFile::Open(rootfilename.c_str());

  // protection against bad file
  if (rootFile == nullptr) {
    ATH_MSG_ERROR(" Was not able to open file: " << rootfilename
                                                 << " ...... aborting");
    return StatusCode::FAILURE;
  }

  //////////////////////////////////////////////////////////////////////////
  //
  // explanation: attempt to loop generally over a file
  // -- if certain SINGALWORD is present -- then this is taken as a signal,
  // that this is another dimension... can be dynamically added.
  // e.g.
  // SFSyst<number>_RunNumber<minRN>-<maxRN>_Nvtx<minNvtx>-<maxNvtx>
  // SFStat_RunNumber<minRN>-<maxRN>_Nvtx<minNvtx>-<maxNvtx>
  // SFCentral_RunNumber<minRN>-<maxRN>_Nvtx<minNvtx>-<maxNvtx>

  //     Then can create a key that will dynamically give us access to a map:
  //    std::map<std::string key, std::vector<TH2 *>> m_SF_SS;     // keys (e.g.
  //    RunNumber223333_319200_Nvtx0_10_Phi1.5_1.6) mapping to vector of SF
  //    histograms --> vector m_SF: 0=nominal, 1=stat, 2,3,4...n=syst
  //     std::map<std::string key, std::vector<TH2 *>> m_SF_OS;     // keys
  //     (e.g. RunNumber223333_319200_Nvtx0_10_Phi1.5_1.6) mapping to vector of
  //     SF histograms --> vector m_SF: 0=nominal, 1=stat, 2,3,4...n=syst
  // TFile*         data/ChMisIDSF_TightLL_FixedCutTight.root
  //  KEY: TH2F     SFCentral_RunNumber296939_311481_SS;1
  //  SFCentral_RunNumber296939_311481_SS KEY: TH2F
  //  SFCentral_RunNumber296939_311481_OS;1 SFCentral_RunNumber296939_311481_OS
  //  KEY: TH2F     STAT_RunNumber296939_311481_SS;1
  //  STAT_RunNumber296939_311481_SS KEY: TH2F STAT_RunNumber296939_311481_OS;1
  //  STAT_RunNumber296939_311481_OS KEY: TH2F
  //  SYST_RunNumber296939_311481_total_SS;1  SYST_RunNumber296939_311481_SS:
  //  total KEY: TH2F     SYST_RunNumber296939_311481_total_OS;1
  //  SYST_RunNumber296939_311481_OS: total

  m_SF_SS.clear();
  m_SF_OS.clear();
  TList* keyListfolder = rootFile->GetListOfKeys();
  std::vector<std::string> names;
  std::set<std::string> set_systematics;

  names.reserve(keyListfolder->GetEntries());
  for (int j = 0; j < keyListfolder->GetEntries(); j++) {
    names.emplace_back((keyListfolder->At(j)->GetName()));
  }
  std::sort(names.begin(), names.end());

  for (unsigned int j = 0; j < names.size(); j++) {

    std::string name = names.at(j);
    ATH_MSG_DEBUG("Got ROOT object with name: " << name);
    if (name.find(Form("SFCentral_")) != std::string::npos) {
      ATH_MSG_VERBOSE("Found name 'SFCentral_' in ROOT object name");
      // Check for opposite-sign (=opposite-charge)
      bool isOS = false;
      if (name.find(Form("_OS")) != std::string::npos) {
        isOS = true;
        ATH_MSG_VERBOSE("Found name '_OS' in ROOT object name");
      }
      if (isOS) {
        std::string histid = (names.at(j));
        histid.erase(0, 10);
        histid.erase(histid.size() - 3, 3); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using histid: " << histid);

        if (histid.find("RunNumber") != std::string::npos) {
          ATH_MSG_VERBOSE("Found name 'RunNumber' in histid");
          std::string runlow = histid;
          runlow.erase(histid.find(Form("RunNumber")), 9);
          runlow.erase(runlow.find('_'), runlow.size());
          m_RunNumbers.push_back(
            static_cast<unsigned int>(atoi(runlow.c_str())));
          std::string runhigh = histid;
          runhigh.erase(histid.find(Form("RunNumber")), 9);
          runhigh.erase(0, runhigh.find('_') + 1);
          m_RunNumbers.push_back(
            static_cast<unsigned int>(atoi(runhigh.c_str())));
        }
        ATH_MSG_VERBOSE("Using histid (OS hid): " << histid);
        m_SF_OS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      } else {
        std::string histid = (names.at(j));
        histid.erase(0, 10);
        histid.erase(histid.size() - 3, 3); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using histid (do we this in ? SS): " << histid);
        m_SF_SS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      }
    } ///// if ( name.find(Form("SFCentral_") ) != std::string::npos)

    /// STAT ERROR
    if (name.find(Form("STAT_")) != std::string::npos) {
      ATH_MSG_VERBOSE("Found name 'STAT_' in ROOT object name");
      bool isOS = false;
      if (name.find(Form("_OS")) != std::string::npos) {
        isOS = true;
        ATH_MSG_VERBOSE("Found name '_OS' in ROOT object name");
      }
      if (isOS) {
        std::string histid = (names.at(j));
        histid.erase(0, 5);
        histid.erase(histid.size() - 3, 3); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using histid: " << histid);

        if (histid.find("RunNumber") != std::string::npos) {
          ATH_MSG_VERBOSE("Found name 'RunNumber' in histid");
          std::string runlow = histid;
          runlow.erase(histid.find(Form("RunNumber")), 9);
          runlow.erase(runlow.find('_'), runlow.size());
          //          m_RunNumbers.push_back( static_cast<unsigned
          //          int>(atoi(runlow.c_str())) );
          std::string runhigh = histid;
          runhigh.erase(histid.find(Form("RunNumber")), 9);
          runhigh.erase(0, runhigh.find('_') + 1);
          //          m_RunNumbers.push_back( static_cast<unsigned
          //          int>(atoi(runhigh.c_str())) );
        }
        ATH_MSG_VERBOSE("Using histid (OS hid): " << histid);
        m_SF_OS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      } else {
        std::string histid = (names.at(j));
        ATH_MSG_VERBOSE("Found  histid: " << histid);
        histid.erase(0, 5);
        histid.erase(histid.size() - 3, 3); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using histid (do we this in ? SS): " << histid);
        m_SF_SS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      }

    } ///// if ( name.find(Form("SYST") ) != std::string::npos)

    /// STAT ERROR
    if (name.find(Form("SYST")) != std::string::npos) {
      ATH_MSG_VERBOSE("Found name 'SYST' in ROOT object name");
      bool isOS = false;
      if (name.find(Form("_OS")) != std::string::npos) {
        isOS = true;
        ATH_MSG_VERBOSE("Found name '_OS' in ROOT object name");
      }
      if (isOS) {
        std::string histid = (names.at(j));
        histid.erase(0, 4);
        histid.erase(histid.size() - 3, 3); // remove _SS, _OS

        std::string sysname = histid;
        sysname.erase(sysname.find('_'), sysname.size());
        set_systematics.insert(sysname);

        histid.erase(0, histid.find('_') + 1); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using syst histid: " << histid);

        if (histid.find("RunNumber") != std::string::npos) {
          std::string runlow = histid;
          runlow.erase(histid.find(Form("RunNumber")), 9);
          runlow.erase(runlow.find('_'), runlow.size());
          //        m_RunNumbers.push_back( static_cast<unsigned
          //        int>(atoi(runlow.c_str())) );
          std::string runhigh = histid;
          runhigh.erase(histid.find(Form("RunNumber")), 9);
          runhigh.erase(0, runhigh.find('_') + 1);
          //      m_RunNumbers.push_back( static_cast<unsigned
          //      int>(atoi(runhigh.c_str())) );
        }
        ATH_MSG_VERBOSE("Using histid (OS hid): " << histid);
        m_SF_OS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      } else {
        std::string histid = (names.at(j));
        histid.erase(0, 4);
        histid.erase(histid.size() - 3, 3);    // remove _SS, _OS
        histid.erase(0, histid.find('_') + 1); // remove _SS, _OS
        ATH_MSG_VERBOSE("Using histid (sys ? SS): " << histid);
        m_SF_SS[histid].push_back((TH2*)rootFile->Get(names.at(j).c_str()));
      }

    } /// end // if ( name.find(Form("SYST") ) != std::string::npos)
  }

  /////////// checks ... --> same vector length... all files there?

  if (m_SF_OS.empty() || m_SF_SS.empty() || m_SF_SS.size() != m_SF_OS.size()) {
    ATH_MSG_ERROR(
      "OS/SS SF vectors not filled or of different size. -- Problem with "
      "files. -- Report to <hn-atlas-EGammaWG@cern.ch>");
    return StatusCode(CP::CorrectionCode::Error);
  }

  m_systematics.insert(m_systematics.end(), set_systematics.begin(), set_systematics.end());

  std::sort(m_RunNumbers.begin(), m_RunNumbers.end());
  ///////////////////////////////////////////////////////////////////////////////////
  // Determine the limits of validity

  /// here: need to use iterator over map!!!
  ATH_MSG_DEBUG("Having m_SF_OS.size() = " << m_SF_OS.size());
  std::map<std::string, std::vector<TH2*>>::iterator it = m_SF_OS.begin();

  // Get the kinematic limits
  m_eta_lowlimit = (*it).second.at(0)->GetYaxis()->GetXmin();
  m_eta_uplimit = (*it).second.at(0)->GetYaxis()->GetXmax();
  ATH_MSG_VERBOSE("|eta| limits " << m_eta_lowlimit << ", " << m_eta_uplimit);

  m_pt_lowlimit = (*it).second.at(0)->GetXaxis()->GetXmin();
  m_pt_uplimit = (*it).second.at(0)->GetXaxis()->GetXmax();
  ATH_MSG_VERBOSE("pt limits " << m_pt_lowlimit << ", " << m_pt_uplimit);

  // Check if the input file is in GeV or MeV
  if (m_pt_uplimit > 1500) {
    ATH_MSG_VERBOSE("Rates in input file are in MeV");
    m_gevmev = 1.;
  } else {
    ATH_MSG_VERBOSE("Rates in input file are in GeV");
    m_gevmev = 0.001;
  }

  // Systematics // dynamic too?
  m_affectingSys = affectingSystematics();

  // Add the recommended systematics to the registry
  if (registerSystematics() != StatusCode::SUCCESS) {
    ATH_MSG_ERROR("(registerSystematics() != CP::SystematicCode::Ok)");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

// =============================================================================
// Athena finalize method
// =============================================================================
StatusCode
CP::ElectronChargeEfficiencyCorrectionTool::finalize()
{
  return StatusCode::SUCCESS;
}

//---------------------------------------------------------------------------------------
// Get the scale factor for the electron
//---------------------------------------------------------------------------------------

//

CP::CorrectionCode
CP::ElectronChargeEfficiencyCorrectionTool::getEfficiencyScaleFactor(
  const xAOD::Electron& ele,
  double& sf) const
{

  // initialize the SF at 1
  sf = 1.0;

  // checking on the truth electron: up to now if this is not a good ele it's
  // returning
  bool goodEle = false;
  CP::CorrectionCode goodEle_result =
    ElectronEfficiencyHelpers::isGoodEle( ele, goodEle); 
  if (goodEle_result != CP::CorrectionCode::Ok) {
    sf = -999.0;
    ATH_MSG_DEBUG("This is the check of goodeleCC in getscalefactor. Scale "
                  "factor set to -999");
    return goodEle_result;
  }

  if (!goodEle) {
    // electron is background electron and should not be corrected
    return CP::CorrectionCode::Ok;
    ATH_MSG_DEBUG("Here goodele is false but CC ok");
  }

  // taking reconstructed variables
  int reco_ele_charge = ele.charge();
  const double ele_pt = ele.pt() * m_gevmev;
  const double ele_eta = std::abs(ele.caloCluster()->etaBE(2));

  // getting the truth charge
  int truth_ele_charge = 9999;
  CP::CorrectionCode charge_result = 
    ElectronEfficiencyHelpers::getEleTruthCharge( ele, truth_ele_charge);
  if (charge_result != CP::CorrectionCode::Ok) {
    sf = -9999.0;
    ATH_MSG_VERBOSE("This is check of geteletruthchargeCC in getscalefactor. "
                    "Scale factor set to -9999");
    return charge_result;
  }

  if (truth_ele_charge == 0) {
    ATH_MSG_DEBUG("Here truth charge is =0!!");
    return CP::CorrectionCode::Ok;
  }

  ATH_MSG_DEBUG("Reco charge = " << reco_ele_charge
                                 << "; Truth charge = " << truth_ele_charge);

  // getting the rates from file....
  float retVal(0.0);

  //////////////////////////////////////
  // here determine, WHICH of the [histid] to choose (after cuuts on runnumber
  // etc....)
  std::string cutRunNumber = "all";

  if (!m_RunNumbers.empty()) {
    unsigned int runnumber = m_defaultRandomRunNumber;
    ATH_MSG_DEBUG("RandomRunNumber: " << runnumber << " "
                                      << m_useRandomRunNumber);
    if (m_useRandomRunNumber) {
      const xAOD::EventInfo* eventInfo =
        evtStore()->retrieve<const xAOD::EventInfo>(m_eventInfoCollectionName);
      if (!eventInfo) {
        ATH_MSG_ERROR("Could not retrieve EventInfo object!");
        sf = 1.0;
        return CP::CorrectionCode::Error;
      }
      static const SG::AuxElement::Accessor<unsigned int> randomrunnumber(
        "RandomRunNumber");
      if (!randomrunnumber.isAvailable(*eventInfo)) {
        sf = 1.0;
        ATH_MSG_WARNING(
          "Pileup tool not run before using ElectronEfficiencyTool! SFs do not "
          "reflect PU distribution in data");
        return CP::CorrectionCode::Error;
      }
      runnumber = randomrunnumber(*(eventInfo));
    }
    ATH_MSG_DEBUG("Number of RunNumbers in file: " << m_RunNumbers.size());
    for (std::size_t r = 0; r < m_RunNumbers.size(); r++) {
      ATH_MSG_DEBUG( " - " << m_RunNumbers.at(r));
    }
    ATH_MSG_VERBOSE("DONE");

    bool isInRunNumberRange = false;
    for ( std::size_t r=0; r<m_RunNumbers.size()-1; r+=2 ){
      // increment by two, run numbers always come in pairs (upper and lower bound specified in the histogram name)

      if ( runnumber >= (unsigned int)m_RunNumbers.at(r) && 
	   runnumber <= (unsigned int)m_RunNumbers.at(r+1) ) {
        cutRunNumber.clear();
        cutRunNumber =
          Form("RunNumber%d_%d", m_RunNumbers.at(r), m_RunNumbers.at(r + 1));
        ATH_MSG_DEBUG("Random run number lies in range " << m_RunNumbers.at(r) << " " << m_RunNumbers.at(r+1));
	isInRunNumberRange = true;
      }
    }

    if (runnumber < m_RunNumbers.at(0) ||
        (runnumber > m_RunNumbers.at(m_RunNumbers.size() - 1))) {
      ATH_MSG_DEBUG("RunNumber " << runnumber << " is not in valid RunNumber Range ");
      sf = 1.0;
      return CP::CorrectionCode::OutOfValidityRange;
    }
    
    if ( !isInRunNumberRange ) {
      return CP::CorrectionCode::OutOfValidityRange;
    }
  }

  // check if electron is within recommendations in eta/Et
  if ( ele_eta < m_eta_lowlimit || ele_eta > m_eta_uplimit ) {

    ATH_MSG_DEBUG("Got an electron outside of the range of eta validity " << ele_eta);
    return CP::CorrectionCode::OutOfValidityRange;
  }

  if ( ele_pt < m_pt_lowlimit ) {

    ATH_MSG_DEBUG("Got an electron outside of the range of pt validity: pt lower than lower limit");
    return CP::CorrectionCode::OutOfValidityRange;
  }

  // Determine WHICH histograms to use here
  const std::vector<TH2*>& SShistograms = m_SF_SS.at(cutRunNumber.c_str());
  const std::vector<TH2*>& OShistograms = m_SF_OS.at(cutRunNumber.c_str());

  // here check OS or SS
  bool isOS = false;

  if (truth_ele_charge * reco_ele_charge > 0)
    isOS = true;

  if (isOS) {
    retVal = this->getChargeFlipRate(ele_eta, ele_pt, OShistograms.at(0), sf);
    if (retVal != 0) {
      sf = -9999.0;
      return CP::CorrectionCode::OutOfValidityRange;
    }
  } else {
    ATH_MSG_DEBUG("Get SS his");
    retVal = this->getChargeFlipRate(ele_eta, ele_pt, SShistograms.at(0), sf);
    if (retVal != 0) {
      sf = -9999.0;
      return CP::CorrectionCode::OutOfValidityRange;
    }
  }

  ATH_MSG_DEBUG("eta: " << ele_eta << "  pt: " << ele_pt);
  ATH_MSG_DEBUG("SF Rates---- . SF: " << sf);

  // Systematics
  // ------------------------------------------------------------------------------------------------------
  double val_stat;

  /// STAT
  if (isOS) {
    retVal =
      this->getChargeFlipRate(ele_eta, ele_pt, OShistograms.at(1), val_stat);
    if (retVal != 0) {
      sf = -9999.0;
      return CP::CorrectionCode::OutOfValidityRange;
    }
  } else {
    ATH_MSG_DEBUG("Get SS his");
    retVal =
      this->getChargeFlipRate(ele_eta, ele_pt, SShistograms.at(1), val_stat);
    if (retVal != 0) {
      sf = -9999.0;
      return CP::CorrectionCode::OutOfValidityRange;
    }
  }

  std::vector<float> systs;
  double val_sys{ 0.0 };
  /// STAT
  for (unsigned int s = 2; s < OShistograms.size(); s++) {
    if (isOS) {
      retVal =
        this->getChargeFlipRate(ele_eta, ele_pt, OShistograms.at(s), val_sys);
      if (retVal != 0) {
        val_sys = -9999.0;
        return CP::CorrectionCode::OutOfValidityRange;
      }
    } else {
      ATH_MSG_DEBUG("Get SS his");
      retVal =
        this->getChargeFlipRate(ele_eta, ele_pt, SShistograms.at(s), val_sys);
      if (retVal != 0) {
        val_sys = -9999.0;
        return CP::CorrectionCode::OutOfValidityRange;
      }
    }
    systs.push_back(static_cast<float>(val_sys));
  }

  ATH_MSG_DEBUG(" ... nominal SF: " << sf);

  if (m_mySysConf.empty()) {
    ATH_MSG_DEBUG(" ... nominal SF: " << sf);
  } else if (*(m_mySysConf.begin()) ==
             SystematicVariation("EL_CHARGEID_STAT", 1)) {
    sf = (sf + (val_stat));
    ATH_MSG_DEBUG("SF after STATup = " << sf);
  } else if (*(m_mySysConf.begin()) ==
             SystematicVariation("EL_CHARGEID_STAT", -1)) {
    sf = (sf - (val_stat));
    ATH_MSG_DEBUG("SF after STATdown = " << sf);
  } else {

    for (unsigned int i = 0; i < m_systematics.size(); i++) {
      if (*(m_mySysConf.begin()) ==
          SystematicVariation(
            Form("EL_CHARGEID_SYS%s", m_systematics.at(i).c_str()), 1)) {
        sf = (sf + (val_sys));
        ATH_MSG_DEBUG("SF after SYSup = " << sf);
      }

      if (*(m_mySysConf.begin()) ==
          SystematicVariation(
            Form("EL_CHARGEID_SYS%s", m_systematics.at(i).c_str()), -1)) {
        sf = (sf - (val_sys));
        ATH_MSG_DEBUG("SF after SYSdown = " << sf);
      }
    }

  }

  return CP::CorrectionCode::Ok;
}

//---------------------------------------------------------------------------------------
// Decorate the electron with the scale factor
//---------------------------------------------------------------------------------------

CP::CorrectionCode
CP::ElectronChargeEfficiencyCorrectionTool::applyEfficiencyScaleFactor(
  const xAOD::Electron& part) const
{
  ATH_MSG_DEBUG(
    "In "
    "CP::ElectronChargeEfficiencyCorrectionTool::applyEfficiencyScaleFactor("
    "const xAOD::IParticle& part) const");
  double sf = 0.0;
  CP::CorrectionCode result = this->getEfficiencyScaleFactor(part, sf);
  // Decorate the electron
  (*m_sfDec)(part) = static_cast<float>(sf);
  return result;
}

// Get the correction rate given pt (E), eta, histogram
float
CP::ElectronChargeEfficiencyCorrectionTool::getChargeFlipRate(
  double eta,
  double pt,
  TH2* hrates,
  double& flipRate) const
{
  ATH_MSG_VERBOSE(" -> in: getChargeFlipRate(" << pt << ", " << eta
                                               << " TH2, double&)");

  if (pt > m_pt_uplimit)
    pt = m_pt_uplimit * 0.999;

  int bin2D = hrates->FindBin(pt, eta);
  flipRate = hrates->GetBinContent(bin2D);

  ATH_MSG_VERBOSE(" -> flipRate is " << flipRate << ", for histogram "
                                     << hrates->GetName());

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
// returns whether this tool is affected by the given systematics
bool
CP::ElectronChargeEfficiencyCorrectionTool::isAffectedBySystematic(
  const SystematicVariation& systematic) const
{

  CP::SystematicSet sys = affectingSystematics();
  return sys.find(systematic) != sys.end();
}

///////////////////////////////////////////////////////////////////////////////////////////
// returns the list of all systematics this tool can be affected by

CP::SystematicSet
CP::ElectronChargeEfficiencyCorrectionTool::affectingSystematics() const
{
  CP::SystematicSet result;
  result.insert(SystematicVariation("EL_CHARGEID_STAT", 1));
  result.insert(SystematicVariation("EL_CHARGEID_STAT", -1));

  for (unsigned int i = 0; i < m_systematics.size(); i++) {
    result.insert(SystematicVariation(
      Form("EL_CHARGEID_SYS%s", m_systematics.at(i).c_str()), 1));
    result.insert(SystematicVariation(
      Form("EL_CHARGEID_SYS%s", m_systematics.at(i).c_str()), -1));
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
// returns the list of all systematics this tool recommends to use
CP::SystematicSet
CP::ElectronChargeEfficiencyCorrectionTool::recommendedSystematics() const
{

  return affectingSystematics();
}

///////////////////////////////////////////////////////////////////////////////////////////
// Gets a SystematicSet and filters it
StatusCode
CP::ElectronChargeEfficiencyCorrectionTool::applySystematicVariation(
  const SystematicSet& systConfig)
{

  if (!SystematicSet::filterForAffectingSystematics(
        systConfig, m_affectingSys, m_mySysConf)) {
    ATH_MSG_ERROR(
      "Unsupported combination of systematics passed to the tool! ");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// Register the systematics with the registry and add them to the recommended list
StatusCode CP::ElectronChargeEfficiencyCorrectionTool::registerSystematics() {
  CP::SystematicRegistry &registry = CP::SystematicRegistry::getInstance();

  if (registry.registerSystematics(*this) != StatusCode::SUCCESS) {
    ATH_MSG_ERROR("Failed to add systematic to list of recommended systematics.");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

/////////////////////////////////////////////////////////////////////////////////////////////
