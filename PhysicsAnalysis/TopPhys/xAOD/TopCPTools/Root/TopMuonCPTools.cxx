/*
   Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TopCPTools/TopMuonCPTools.h"

#include <map>
#include <string>

// Top includes
#include "TopConfiguration/TopConfig.h"
#include "TopEvent/EventTools.h"

// PathResolver include(s):
#include "PathResolver/PathResolver.h"

// Muon include(s):
#include "MuonMomentumCorrections/MuonCalibTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "MuonEfficiencyCorrections/MuonTriggerScaleFactors.h"
#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"

namespace top {
  MuonCPTools::MuonCPTools(const std::string& name) :
    asg::AsgTool(name) {
    declareProperty("config", m_config);

    declareProperty("MuonMomentumCalibrationTool", m_muonMomentumCalibrationTool);

    declareProperty("MuonSelectionTool", m_muonSelectionTool);
    declareProperty("MuonSelectionToolLoose", m_muonSelectionToolLoose);
    declareProperty("MuonSelectionToolVeryLooseVeto", m_muonSelectionToolVeryLooseVeto);

    declareProperty("MuonEfficiencyCorrectionsTool", m_muonEfficiencyCorrectionsTool);
    declareProperty("MuonEfficiencyCorrectionsToolLoose", m_muonEfficiencyCorrectionsToolLoose);
    declareProperty("MuonEfficiencyCorrectionsToolIso", m_muonEfficiencyCorrectionsToolIso);
    declareProperty("MuonEfficiencyCorrectionsToolLooseIso", m_muonEfficiencyCorrectionsToolLooseIso);
    declareProperty("MuonEfficiencyCorrectionsToolBadMuonVeto", m_muonEfficiencyCorrectionsToolBadMuonVeto);

    declareProperty("SoftMuonSelectionTool", m_softmuonSelectionTool);
    declareProperty("SoftMuonEfficiencyCorrectionsTool", m_softmuonEfficiencyCorrectionsTool);

  }

  StatusCode MuonCPTools::initialize() {
    ATH_MSG_INFO("top::MuonCPTools initialize...");

    if (m_config->isTruthDxAOD()) {
      ATH_MSG_INFO("top::MuonCPTools: no need to initialise anything on truth DxAOD");
      return StatusCode::SUCCESS;
    }

    if (!m_config->useMuons() && !m_config->useSoftMuons()) {
      ATH_MSG_INFO("top::MuonCPTools: no need to initialise anything since not using muons");
      return StatusCode::SUCCESS;
    }

    if (m_config->makeAllCPTools()) {// skiping calibrations on mini-xAODs
      top::check(setupCalibration(), "Failed to setup Muon calibration tools");
    }
    if (m_config->isMC()) {// scale-factors are only for MC
      top::check(setupScaleFactors(), "Failed to setup Muon scale-factor tools");
    }
    return StatusCode::SUCCESS;
  }

  StatusCode MuonCPTools::setupCalibration() {
    ///-- Selection --///
    m_muonSelectionTool = setupMuonSelectionTool("MuonSelectionTool",
                                                 m_config->muonQuality(),
                                                 m_config->muonEtacut(),
                                                 m_config->muonUseMVALowPt(),
                                                 m_config->muonUse2stationMuonsHighPt());
    m_muonSelectionToolLoose = setupMuonSelectionTool("MuonSelectionToolLoose",
                                                      m_config->muonQualityLoose(),
                                                      m_config->muonEtacut(),
                                                      m_config->muonUseMVALowPtLoose(),
                                                      m_config->muonUse2stationMuonsHighPtLoose());
    // the following is needed to make sure all muons for which d0sig is calculated are at least Loose
    m_muonSelectionToolVeryLooseVeto = setupMuonSelectionTool("MuonSelectionToolVeryLooseVeto",
                                                              "Loose",
                                                              2.5,
                                                              m_config->muonUseMVALowPt(),
                                                              m_config->muonUse2stationMuonsHighPt());
    ///-- Calibration and smearing --///  ---> now passing the flags (true/false) to CalibAndSmearingTool
    m_muonMomentumCalibrationTool = setupMuonCalibrationAndSmearingTool("MuonMomentumCalibrationTool", 
                                                                        m_config->muonMuonDoExtraSmearingHighPt(),
                                                                        m_config->muonMuonDoSmearing2stationHighPt());
    //now the soft muon part
    if (m_config->useSoftMuons()) {
      m_softmuonSelectionTool = setupMuonSelectionTool("SoftMuonSelectionTool",
                                                       m_config->softmuonQuality(),
                                                       m_config->softmuonEtacut(),
                                                       m_config->softmuonUseMVALowPt(),
                                                       false);
    }

    return StatusCode::SUCCESS;
  }

  StatusCode MuonCPTools::setupScaleFactors() {
    // Setup muon SF tools
    // However if we are running on data- we don't need these,
    // so carry on.
    if (!m_config->isMC()) return StatusCode::SUCCESS;

    /************************************************************
    *
    * Muon Scale Factors:
    *    muonSF = trigSF*effSF*isoSF*TTVASF
    *
    ************************************************************/

    /************************************************************
    * Trigger Scale Factors:
    *    setup trigger SFs for nominal and 'loose' muon WPs
    *    recommendation for EOYE not to pass any isolation to tool
    *    as SFs very similar for all WPs.
    ************************************************************/

    // In R21 now, we only need one instance of the tool
    // and do not need to set the year as it is handled
    // internally with PRW tool
    m_muonTriggerScaleFactors_R21
      = setupMuonTrigSFTool("MuonTriggerScaleFactors_R21",
                            m_config->muonQuality());
    m_muonTriggerScaleFactorsLoose_R21
      = setupMuonTrigSFTool("MuonTriggerScaleFactorsLoose_R21",
                            m_config->muonQualityLoose());

    /************************************************************
    * Efficiency Scale Factors:
    *    setup muon efficiency SFs for the nominal and
    *    'loose' muon WPs.
    ************************************************************/
    
    //if !Use2stationMuonsHighPt, HighPt -> HighPt3Layers
    //if UseMVALowPt, LowPt -> LowPtMVA
    std::string muonQuality_name = m_config->muonQuality();
    if (m_config->muonQuality() == "HighPt" && !(m_config->muonUse2stationMuonsHighPt()) ) muonQuality_name = "HighPt3Layers";
    if (m_config->muonQuality() == "LowPt" && m_config->muonUseMVALowPt()) muonQuality_name = "LowPtMVA";
    m_muonEfficiencyCorrectionsTool
      = setupMuonSFTool("MuonEfficiencyScaleFactorsTool",
                        muonQuality_name, false);

    std::string muonQualityLoose_name = m_config->muonQualityLoose();
    if (m_config->muonQualityLoose() == "HighPt" && !(m_config->muonUse2stationMuonsHighPtLoose()) ) muonQualityLoose_name = "HighPt3Layers";
    if (m_config->muonQualityLoose() == "LowPt" && m_config->muonUseMVALowPtLoose()) muonQualityLoose_name = "LowPtMVA";
    m_muonEfficiencyCorrectionsToolLoose
      = setupMuonSFTool("MuonEfficiencyScaleFactorsToolLoose",
                        muonQualityLoose_name, false);

    if (m_config->muonQuality() == "HighPt" || m_config->muonQualityLoose() == "HighPt") {
      m_muonEfficiencyCorrectionsToolLoose
        = setupMuonSFTool("MuonEfficiencyScaleFactorsToolBadMuonVeto",
                          "BadMuonVeto_HighPt", false);
    }

    //now the soft muon part
    std::string softmuonQuality_name = m_config->softmuonQuality();
    if (m_config->softmuonQuality() == "LowPt" && m_config->softmuonUseMVALowPt()) softmuonQuality_name = "LowPtMVA";
    if (m_config->useSoftMuons()) {
      m_softmuonEfficiencyCorrectionsTool
        = setupMuonSFTool("SoftMuonEfficiencyScaleFactorsTool",
	                  softmuonQuality_name, false);
    }
    
    /************************************************************
     * Isolation Scale Factors:
    *    setup muon isolation SFs for the nominal and 'loose'
    *    muons
    *
    *    Note: if isolation WP is None, then don't setup the tool
    ************************************************************/
    // If we don't want isolation then we don't need the tool
    if (m_config->muonIsolationSF() != "None") {
      // Add iso as a suffix (see above for consistency between tools :) )
      std::string muon_isolation = m_config->muonIsolationSF() + "Iso";
      m_muonEfficiencyCorrectionsToolIso =
        setupMuonSFTool("MuonEfficiencyScaleFactorsToolIso",
                        muon_isolation, true);
    }

    // Do we have isolation on our loose muons? If not no need for the tool...
    if (m_config->muonIsolationSFLoose() != "None") {
      // Add iso as a suffix (see above for consistency between tools :) )
      std::string muon_isolation = m_config->muonIsolationSFLoose() + "Iso";
      m_muonEfficiencyCorrectionsToolLooseIso =
        setupMuonSFTool("MuonEfficiencyScaleFactorsToolLooseIso",
                        muon_isolation, true);
    }

    /************************************************************
    * Muon TTVA SF:
    *    Track-to-vertex association. This depends on whether or
    *    not we apply the tracking groups recommended impact
    *    parameter cuts to associate muon to vertex.
    ************************************************************/
    m_muonEfficiencyCorrectionsToolTTVA
      = setupMuonSFTool("MuonEfficiencyScaleFactorsToolTTVA",
                        "TTVA", false);

    // WARNING - The PromptLeptonIsolation scale factors are only derived with respect to the loose PID
    //         - Hence we need to fail if this has occured
    if ((m_config->muonQuality() != "Loose" && m_config->muonIsolationSF() == "PromptLepton")
        || (m_config->muonQualityLoose() != "Loose" && m_config->muonIsolationSFLoose() == "PromptLepton")) {
      ATH_MSG_ERROR(
        "Cannot use PromptLeptonIsolation on muons without using Loose quality - Scale factors are not available");
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  CP::IMuonSelectionTool*
  MuonCPTools::setupMuonSelectionTool(const std::string& name, const std::string& quality, double max_eta, const bool& useMVALowPt, const bool& use2stationMuonsHighPt) {
    std::map<std::string, int> muon_quality_map = {
      {"Tight", 0}, {"Medium", 1}, {"Loose", 2}, {"VeryLoose", 3}, {"HighPt", 4}, {"LowPt", 5}
    };
    int qual_int;
    try {
      qual_int = muon_quality_map.at(quality);
    } catch (const std::out_of_range& oor_exc) {
      ATH_MSG_ERROR("\n Invalid muon quality ("
                    + quality + ") for " + name
                    + ". Valid options are: "
                      " \n\t- Tight"
                      " \n\t- Medium"
                      " \n\t- Loose"
                      " \n\t- VeryLoose"
                      " \n\t- HighPt"
                      " \n\t- LowPt");
      throw;  // Re-throw
    }

    CP::IMuonSelectionTool* tool = nullptr;
    if (asg::ToolStore::contains<CP::IMuonSelectionTool>(name)) {
      tool = asg::ToolStore::get<CP::IMuonSelectionTool>(name);
    } else {
      tool = new CP::MuonSelectionTool(name);
      top::check(asg::setProperty(tool, "IsRun3Geo", m_config->isRun3()),
                 "Failed to set IsRun3Geo for " + name);
      top::check(asg::setProperty(tool, "MuQuality", qual_int),
                 "Failed to set MuQuality for " + name);
      top::check(asg::setProperty(tool, "MaxEta", max_eta),
                 "Failed to set MaxEta for " + name);
      top::check(asg::setProperty(tool, "UseMVALowPt", useMVALowPt),
                 "Failed to set UseMVALowPt for " + name + " tool");
      top::check(asg::setProperty(tool, "Use2stationMuonsHighPt", use2stationMuonsHighPt),
                 "Failed to set Use2stationMuonsHighPt for " + name + " tool");
      top::check(tool->initialize(), "Failed to initialize " + name);
    }
    return tool;
  }

  CP::IMuonTriggerScaleFactors*
  MuonCPTools::setupMuonTrigSFTool(const std::string& name, const std::string& quality) {
    CP::IMuonTriggerScaleFactors* tool = nullptr;
    if (asg::ToolStore::contains<CP::IMuonTriggerScaleFactors>(name)) {
      tool = asg::ToolStore::get<CP::IMuonTriggerScaleFactors>(name);
    } else {
      tool = new CP::MuonTriggerScaleFactors(name);
      top::check(asg::setProperty(tool, "MuonQuality", quality),
                 "Failed to set MuonQuality for " + name);
      top::check(asg::setProperty(tool, "AllowZeroSF", true),
                 "Failed to set AllowZeroSF for " + name);
      if (m_config->muonSFCustomInputFolderTrigger() != " ") {
        top::check(asg::setProperty(tool, "CustomInputFolder", m_config->muonSFCustomInputFolderTrigger()),
                   "Failed to set CustomInputFolder property for MuonTriggerScaleFactors tool");
      }
      if (m_config->muonForcePeriod() != " ") {
        top::check(asg::setProperty(tool, "forcePeriod", m_config->muonForcePeriod()),
                   "Failed to set forcePeriod property for MuonTriggerScaleFactors tool");
      }
      if (m_config->muonForceYear() != -1) {
        top::check(asg::setProperty(tool, "forceYear", m_config->muonForceYear()),
                   "Failed to set forceYear property for MuonTriggerScaleFactors tool");
      }
      top::check(tool->initialize(), "Failed to init. " + name);
    }
    return tool;
  }

  CP::IMuonEfficiencyScaleFactors*
  MuonCPTools::setupMuonSFTool(const std::string& name, const std::string& WP, const bool isIso) {
    CP::IMuonEfficiencyScaleFactors* tool = nullptr;
    if (asg::ToolStore::contains<CP::IMuonEfficiencyScaleFactors>(name)) {
      tool = asg::ToolStore::get<CP::MuonEfficiencyScaleFactors>(name);
    } else {
      tool = new CP::MuonEfficiencyScaleFactors(name);
      top::check(asg::setProperty(tool, "WorkingPoint", WP),
                 "Failed to set WP for " + name + " tool");
      top::check(asg::setProperty(tool, "CloseJetDRDecorator", "dRMuJet_AT_usingWeirdNameToAvoidUsingOnTheFlyCalculation"), 
                 "Failed to set WP for " + name + " tool"); //in this way we'll only read the dR(mu,jet) from the derivation, IF the variable is there, but we'll not use on-the-fly calculation, which is tricky in AT
      if (m_config->muonSFCustomInputFolder() != " ") {
        top::check(asg::setProperty(tool, "CustomInputFolder", m_config->muonSFCustomInputFolder()),
                   "Failed to set CustomInputFolder property for MuonEfficiencyScaleFactors tool");
      }
      if (m_config->isRun3()) {
        top::check(asg::setProperty(tool, "CalibrationRelease", "220817_Preliminary_r22run3"),
                   "Failed to set CalibrationRelease property for MuonEfficiencyScaleFactors tool");
      }

      if (!isIso) {
        top::check(asg::setProperty(tool, "BreakDownSystematics", m_config->muonBreakDownSystematics()), 
                  "Failed to set BreakDownSystematics for " + name + " tool");
      }
      top::check(tool->initialize(),
                 "Failed to set initialize " + name);
    }
    return tool;
  }


  CP::IMuonCalibrationAndSmearingTool*
  MuonCPTools::setupMuonCalibrationAndSmearingTool(const std::string& name, const bool& doExtraSmearingHighPt, const bool& do2StationsHighPt) {
    CP::IMuonCalibrationAndSmearingTool* tool = nullptr;
    if (asg::ToolStore::contains<CP::IMuonCalibrationAndSmearingTool>(name)) {
      tool = asg::ToolStore::get<CP::IMuonCalibrationAndSmearingTool>(name);
    } else {
      tool = new CP::MuonCalibTool(name);

      top::check(asg::setProperty(tool, "doExtraSmearing", doExtraSmearingHighPt),
                 "Failed to set doExtraSmearing for " + name + " tool");
      top::check(asg::setProperty(tool, "do2StationsHighPt", do2StationsHighPt),
                 "Failed to set do2StationsHighPt for " + name + " tool");
      top::check(asg::setProperty(tool, "IsRun3Geo", m_config->isRun3()),
                 "Failed to set IsRun3Geo for " + name + " tool");
      if (m_config->isMC() && m_config->forceRandomRunNumber() > 0) {
        top::check(asg::setProperty(tool, "useRandomRunNumber", false),
                   "Failed to set useRandomRunNumber for " + name + " tool");
      }
      if (m_config->muonCalibMode() == "correctData_CB")
        top::check(asg::setProperty(tool, "calibMode", CP::MuonCalibTool::CalibMode::correctData_CB), "Failed to set calibrationMode for " + name + " tool");
      else if (m_config->muonCalibMode() == "correctData_IDMS")
        top::check(asg::setProperty(tool, "calibMode", CP::MuonCalibTool::CalibMode::correctData_IDMS), "Failed to set calibrationMode for " + name + " tool");
      else if (m_config->muonCalibMode() == "notCorrectData_IDMS")
        top::check(asg::setProperty(tool, "calibMode", CP::MuonCalibTool::CalibMode::notCorrectData_IDMS), "Failed to set calibrationMode for " + name + " tool");
      top::check(tool->initialize(),
                 "Failed to set initialize " + name);
    }
    return tool;
  }
}  // namespace top
