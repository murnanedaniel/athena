/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TAUANALYSISTOOLS_COMMONSMEARINGTOOL_H
#define TAUANALYSISTOOLS_COMMONSMEARINGTOOL_H

/*
  author: Dirk Duschinger
  mail: dirk.duschinger@cern.ch
  documentation in: ../README.rst
*/

// Framework include(s):
#include "AsgTools/AsgMetadataTool.h"
#include "AsgTools/AnaToolHandle.h"

// EDM include(s):
#include "xAODTau/TauJet.h"
#include "xAODTruth/TruthParticle.h"
#include "PATInterfaces/CorrectionCode.h"

#include "tauRecTools/TauCombinedTES.h"

// Local include(s):
#include "TauAnalysisTools/Enums.h"
#include "TauAnalysisTools/ITauSmearingTool.h"
#include "TauAnalysisTools/HelperFunctions.h"

// ROOT include(s):
#include "TFile.h"
#include "TH1.h"

// tauRecTools include(s)
#include "tauRecTools/ITauToolBase.h"

namespace TauAnalysisTools
{

class CommonSmearingTool
  : public virtual ITauSmearingTool
  , public asg::AsgMetadataTool
{
  /// Create a proper constructor for Athena
  ASG_TOOL_CLASS( CommonSmearingTool, TauAnalysisTools::ITauSmearingTool )

public:

  CommonSmearingTool(const std::string& sName);

  ~CommonSmearingTool();

  virtual StatusCode initialize();

  // CommonSmearingTool pure virtual public functionality
  //__________________________________________________________________________

  /// Apply the correction on a modifyable object
  virtual CP::CorrectionCode applyCorrection( xAOD::TauJet& xTau ) const;

  /// Create a corrected copy from a constant tau
  virtual CP::CorrectionCode correctedCopy( const xAOD::TauJet& xTau,
      xAOD::TauJet*& xTauCopy) const;

  /// returns: whether this tool is affected by the given systematics
  virtual bool isAffectedBySystematic( const CP::SystematicVariation& systematic ) const;

  /// returns: the list of all systematics this tool can be affected by
  virtual CP::SystematicSet affectingSystematics() const;

  /// returns: the list of all systematics this tool recommends to use
  virtual CP::SystematicSet recommendedSystematics() const;

  /// configure this tool for the given list of systematic variations.  any
  /// requested systematics that are not affecting this tool will be silently ignored
  virtual StatusCode applySystematicVariation ( const CP::SystematicSet& sSystematicSet);

protected:

  std::string ConvertProngToString(const int iProngness) const;

  std::map<std::string, TH1*> m_mSF;
  std::unordered_map < CP::SystematicSet, std::string > m_mSystematicSets;
  const CP::SystematicSet* m_sSystematicSet;
  std::map<std::string, int> m_mSystematics;
  std::map<std::string, std::string> m_mSystematicsHistNames;

  // std::function<double> m_fX;
  double (*m_fX)(const xAOD::TauJet& xTau);
  double (*m_fY)(const xAOD::TauJet& xTau);

  template<class T>
  void ReadInputs(TFile* fFile, std::map<std::string, T>& mMap);

  virtual CP::CorrectionCode getValue(const std::string& sHistName,
                                      const xAOD::TauJet& xTau,
                                      double& dEfficiencyScaleFactor) const;

  void generateSystematicSets();

  std::string m_sInputFilePath;
  bool m_bIsData;
  bool m_bIsConfigured;
  bool m_bSkipTruthMatchCheck;
  bool m_bApplyFading;
  bool m_bApplyMVATESQualityCheck;
  bool m_bApplyInsituCorrection;

  asg::AnaToolHandle<ITauToolBase> m_tTauCombinedTES;

  TruthMatchedParticleType m_eCheckTruth;
  CP::SystematicSet m_sAffectingSystematics;
  CP::SystematicSet m_sRecommendedSystematics;
  
private:

  // Execute at each event
  virtual StatusCode beginEvent();

};
} // namespace TauAnalysisTools

#endif // TAUANALYSISTOOLS_COMMONSMEARINGTOOL_H
