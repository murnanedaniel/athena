/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

// Dear emacs, this is -*-c++-*-

/**
   @class AsgElectronIsEMSelector
   @brief Electron isEM selector

   @author Jovan Mitrevski (UCSC) Karsten Koeneke (CERN)
   @date   Dec 2011 - Fab 2012

   11-MAR-2014 convert to ASG tool (Jovan Mitrevski)

*/

#include "ElectronPhotonSelectorTools/AsgElectronIsEMSelector.h"
#include "AsgTools/CurrentContext.h"
#include "EGSelectorConfigurationMapping.h"
#include "EgammaAnalysisHelpers/AsgEGammaConfigHelper.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "PathResolver/PathResolver.h"
#include "TElectronIsEMSelector.h"
#include "TEnv.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/Photon.h"
#include "xAODTracking/TrackParticle.h"
#include <cstdint>

//=============================================================================
// Standard constructor
//=============================================================================
AsgElectronIsEMSelector::AsgElectronIsEMSelector(const std::string& myname)
  : AsgTool(myname)
  , m_configFile("")
  , m_rootTool(nullptr)
  , m_useF3core(false)
{

  m_rootTool = new Root::TElectronIsEMSelector(myname.c_str());

  declareProperty("WorkingPoint", m_WorkingPoint = "", "The Working Point");
  declareProperty("ConfigFile",
                  m_configFile = "",
                  "The config file to use (if not setting cuts one by one)");

  // Name of the quality to use
  declareProperty(
    "isEMMask",
    m_rootTool->m_isEMMask =
      egammaPID::EgPidUndefined, // All pass by default, if not specified
    "The mask to use");

  declareProperty("useF3core", m_useF3core = false, "Cut on f3 or f3core?");

  // for the trigger needs:
  declareProperty("caloOnly",
                  m_caloOnly = false,
                  "Flag to tell the tool if its a calo only cutbase");
  declareProperty("trigEtTh", m_trigEtTh = -999., "Trigger threshold");
}

//=============================================================================
// Standard destructor
//=============================================================================
AsgElectronIsEMSelector::~AsgElectronIsEMSelector()
{
  delete m_rootTool;
}

StatusCode
AsgElectronIsEMSelector::initialize()
{
  // The standard status code
  StatusCode sc = StatusCode::SUCCESS;

  if (!m_WorkingPoint.empty()) {
    m_configFile = AsgConfigHelper::findConfigFile(
      m_WorkingPoint, EgammaSelectors::ElectronCutPointToConfFile);
  }

  // find the file and read it in
  std::string filename = PathResolverFindCalibFile(m_configFile);
  if (filename.empty()) {
    ATH_MSG_ERROR("Could not locate " << m_configFile);
    sc = StatusCode::FAILURE;
    return sc;
  }
  ATH_MSG_INFO("Configfile to use  " << m_configFile);
  TEnv env;
  env.ReadFile(filename.c_str(), kEnvLocal);

  ///------- Read in the TEnv config ------///

  // Override the mask via the config only if it is not set
  if (m_rootTool->m_isEMMask == egammaPID::EgPidUndefined) {
    int default_mask = static_cast<int>(egammaPID::EgPidUndefined);
    int mask(env.GetValue("isEMMask", default_mask));
    m_rootTool->m_isEMMask = static_cast<unsigned int>(mask);
  }
  //
  // From here on the conf ovverides all other properties
  bool useTRTOutliers(env.GetValue("useTRTOutliers", true));
  m_rootTool->m_useTRTOutliers = useTRTOutliers;
  bool useTRTXenonHits(env.GetValue(" useTRTXenonHits", false));
  m_rootTool->m_useTRTXenonHits = useTRTXenonHits;

  ///------- Use helpers to read in the cut arrays ------///
  m_rootTool->m_cutBinEta = AsgConfigHelper::HelperFloat("CutBinEta", env);
  m_rootTool->m_cutBinET = AsgConfigHelper::HelperFloat("CutBinET", env);
  m_rootTool->m_cutF1 = AsgConfigHelper::HelperFloat("CutF1", env);
  m_rootTool->m_cutHadLeakage =
    AsgConfigHelper::HelperFloat("CutHadLeakage", env);
  m_rootTool->m_cutReta37 = AsgConfigHelper::HelperFloat("CutReta37", env);
  m_rootTool->m_cutRphi33 = AsgConfigHelper::HelperFloat("CutRphi33", env);
  m_rootTool->m_cutWeta2c = AsgConfigHelper::HelperFloat("CutWeta2c", env);
  m_rootTool->m_cutDeltaEmax2 =
    AsgConfigHelper::HelperFloat("CutDeltaEmax2", env);
  m_rootTool->m_cutDeltaE = AsgConfigHelper::HelperFloat("CutDeltaE", env);
  m_rootTool->m_cutDEmaxs1 = AsgConfigHelper::HelperFloat("CutDEmaxs1", env);
  m_rootTool->m_cutDeltaE = AsgConfigHelper::HelperFloat("CutDeltaE", env);
  m_rootTool->m_cutWtot = AsgConfigHelper::HelperFloat("CutWtot", env);
  m_rootTool->m_cutWeta1c = AsgConfigHelper::HelperFloat("CutWeta1c", env);
  m_rootTool->m_cutFracm = AsgConfigHelper::HelperFloat("CutFracm", env);
  m_rootTool->m_cutF3 = AsgConfigHelper::HelperFloat("CutF3", env);
  m_rootTool->m_cutBL = AsgConfigHelper::HelperInt("CutBL", env);
  m_rootTool->m_cutPi = AsgConfigHelper::HelperInt("CutPi", env);
  m_rootTool->m_cutSi = AsgConfigHelper::HelperInt("CutSi", env);
  m_rootTool->m_cutA0 = AsgConfigHelper::HelperFloat("CutA0", env);
  m_rootTool->m_cutA0Tight = AsgConfigHelper::HelperFloat("CutA0Tight", env);
  m_rootTool->m_cutDeltaEta = AsgConfigHelper::HelperFloat("CutDeltaEta", env);
  m_rootTool->m_cutDeltaEtaTight =
    AsgConfigHelper::HelperFloat("CutDeltaEtaTight", env);
  m_rootTool->m_cutminDeltaPhi =
    AsgConfigHelper::HelperFloat("CutminDeltaPhi", env);
  m_rootTool->m_cutmaxDeltaPhi =
    AsgConfigHelper::HelperFloat("CutmaxDeltaPhi", env);
  m_rootTool->m_cutminEp = AsgConfigHelper::HelperFloat("CutminEp", env);
  m_rootTool->m_cutmaxEp = AsgConfigHelper::HelperFloat("CutmaxEp", env);
  m_rootTool->m_cutBinEta_TRT =
    AsgConfigHelper::HelperFloat("CutBinEta_TRT", env);
  m_rootTool->m_cutBinET_TRT =
    AsgConfigHelper::HelperFloat("CutBinET_TRT", env);
  m_rootTool->m_cutNumTRT = AsgConfigHelper::HelperFloat("CutNumTRT", env);
  m_rootTool->m_cutTRTRatio = AsgConfigHelper::HelperFloat("CutTRTRatio", env);
  m_rootTool->m_cutTRTRatio90 =
    AsgConfigHelper::HelperFloat("CutTRTRatio90", env);
  m_rootTool->m_cutEProbabilityHT =
    AsgConfigHelper::HelperFloat("CutEProbabilityHT", env);

  ATH_MSG_INFO("operating point : " << this->getOperatingPointName()
                                    << " with mask: "
                                    << m_rootTool->m_isEMMask);

  // Get the message level and set the underlying ROOT tool message level
  // accordingly
  m_rootTool->msg().setLevel(this->msg().level());

  // We need to initialize the underlying ROOT TSelectorTool
  if (m_rootTool->initialize().isFailure()) {
    ATH_MSG_ERROR("Could not initialize the TElectronIsEMSelector!");
    sc = StatusCode::FAILURE;
    return sc;
  }

  return sc;
}

//=============================================================================
// return the accept info object
//=============================================================================

const asg::AcceptInfo&
AsgElectronIsEMSelector::getAcceptInfo() const
{
  return m_rootTool->getAcceptInfo();
}

//=============================================================================
// The main accept method: the actual cuts are applied here
//=============================================================================
asg::AcceptData
AsgElectronIsEMSelector::accept(const xAOD::IParticle* part) const
{
  return accept(Gaudi::Hive::currentContext(), part);
}

asg::AcceptData
AsgElectronIsEMSelector::accept(const EventContext& ctx,
                                const xAOD::IParticle* part) const
{

  if (part->type() == xAOD::Type::Electron ||
      part->type() == xAOD::Type::Photon) {
    return accept(ctx, static_cast<const xAOD::Egamma*>(part));
  }

  ATH_MSG_ERROR(
    "AsgElectronIsEMSelector::could not convert argument to Electron/Photon");
  return m_rootTool->accept();
}

asg::AcceptData
AsgElectronIsEMSelector::accept(const EventContext& ctx,
                                const xAOD::Egamma* eg) const
{

  if (eg) {
    unsigned int isEM = ~0;
    StatusCode sc = execute(ctx, eg, isEM);
    if (sc.isFailure()) {
      ATH_MSG_ERROR("could not calculate isEM");
      return m_rootTool->accept();
    }
    return m_rootTool->fillAccept(isEM);
  }

  ATH_MSG_ERROR("AsgElectronIsEMSelector::accept was given a bad argument");
  return m_rootTool->accept();
}

asg::AcceptData
AsgElectronIsEMSelector::accept(const EventContext& ctx,
                                const xAOD::Electron* el) const
{
  return accept(ctx, static_cast<const xAOD::Egamma*>(el));
}

asg::AcceptData
AsgElectronIsEMSelector::accept(const EventContext& ctx,
                                const xAOD::Photon* ph) const
{
  return accept(ctx, static_cast<const xAOD::Egamma*>(ph));
}

//=============================================================================
/// Get the name of the current operating point
//=============================================================================
std::string
AsgElectronIsEMSelector::getOperatingPointName() const
{

  if (!m_WorkingPoint.empty()) {
    return m_WorkingPoint;
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronLoosePP) {
    return "Loose";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronMediumPP) {
    return "Medium";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronTightPP) {
    return "Tight";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronLoose1) {
    return "Loose1";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronMedium1) {
    return "Medium1";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronTight1) {
    return "Tight1";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronLooseHLT) {
    return "LooseHLT";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronMediumHLT) {
    return "MediumHLT";
  }
  if (m_rootTool->m_isEMMask == egammaPID::ElectronTightHLT) {
    return "TightHLT";
  }
  if (m_rootTool->m_isEMMask == 0) {
    return "0 No cuts applied";
  }

  ATH_MSG_INFO("Didn't recognize the given operating point with mask: "
               << m_rootTool->m_isEMMask);
  return "";
}

// ==============================================================
StatusCode
AsgElectronIsEMSelector::execute(const EventContext& ctx,
                                 const xAOD::Egamma* eg,
                                 unsigned int& isEM) const
{
  //
  // Particle identification for electrons based on cuts
  //
  (void)ctx;
  // initialisation
  isEM = 0;
  // protection against null pointer
  if (eg == nullptr) {
    // if object is bad then use the bit for "bad eta"
    ATH_MSG_ERROR("exiting because el is NULL");
    isEM = (0x1 << egammaPID::ClusterEtaRange_Electron);
    return StatusCode::SUCCESS;
  }
  // retrieve associated cluster
  const xAOD::CaloCluster* cluster = eg->caloCluster();
  if (cluster == nullptr) {
    // if object is bad then use the bit for "bad eta"
    ATH_MSG_ERROR("exiting because cluster is NULL");
    isEM = (0x1 << egammaPID::ClusterEtaRange_Electron);
    return StatusCode::SUCCESS;
  }
  // eta position in second sampling
  const float eta2 = fabsf(cluster->etaBE(2));
  // energy in calorimeter
  const double energy = cluster->e();
  // transverse energy of the electron (using the track eta)
  // const double et = el->pt();
  double et = (cosh(eta2) != 0.) ? energy / cosh(eta2) : 0.;
  ;

  // see if we have an electron, with track, for eta
  const xAOD::Electron* el = nullptr;
  if (eg->type() == xAOD::Type::Electron) {
    el = static_cast<const xAOD::Electron*>(eg);
  }
  if (el && el->trackParticle() && !m_caloOnly) {
    et = (cosh(el->trackParticle()->eta()) != 0.)
           ? energy / cosh(el->trackParticle()->eta())
           : 0.;
  }

  // Call the calocuts using the egamma object
  isEM = calocuts_electrons(eg, eta2, et, m_trigEtTh, 0);

  // Call the calo cuts using the el , if available and we want to apply them
  if (el && el->trackParticle() && !m_caloOnly) {
    isEM = TrackCut(el, eta2, et, energy, isEM);
  }

  return StatusCode::SUCCESS;
}

// ======================================================================
unsigned int
AsgElectronIsEMSelector::calocuts_electrons(const xAOD::Egamma* eg,
                                            float eta2,
                                            double et,
                                            double trigEtTh,
                                            unsigned int iflag) const
{

  //
  // apply cut-based selection based on calo information
  // eg : xAOD::Electron object
  // trigETthr : threshold in ET to apply the cuts at trigger level
  // iflag: the starting isEM
  //

  float Reta(0);
  float Rphi(0);
  float Rhad1(0);
  float Rhad(0);
  float e277(0);
  float weta1c(0);
  float weta2c(0);
  float f1(0);
  float emax2(0);
  float Eratio(0);
  float DeltaE(0);
  float wtot(0);
  float fracm(0);
  float f3(0);

  bool allFound = true;
  // Reta
  allFound =
    allFound && eg->showerShapeValue(Reta, xAOD::EgammaParameters::Reta);
  // Rphi
  allFound =
    allFound && eg->showerShapeValue(Rphi, xAOD::EgammaParameters::Rphi);
  // transverse energy in 1st scintillator of hadronic calorimeter
  allFound =
    allFound && eg->showerShapeValue(Rhad1, xAOD::EgammaParameters::Rhad1);
  // transverse energy in hadronic calorimeter
  allFound =
    allFound && eg->showerShapeValue(Rhad, xAOD::EgammaParameters::Rhad);
  // E(7*7) in 2nd sampling
  allFound =
    allFound && eg->showerShapeValue(e277, xAOD::EgammaParameters::e277);
  // shower width in 3 strips in 1st sampling
  allFound =
    allFound && eg->showerShapeValue(weta1c, xAOD::EgammaParameters::weta1);
  // shower width in 2nd sampling
  allFound =
    allFound && eg->showerShapeValue(weta2c, xAOD::EgammaParameters::weta2);
  // fraction of energy reconstructed in the 1st sampling
  allFound = allFound && eg->showerShapeValue(f1, xAOD::EgammaParameters::f1);
  // E of 2nd max between max and min in strips
  allFound =
    allFound && eg->showerShapeValue(Eratio, xAOD::EgammaParameters::Eratio);
  // E of 1st max in strips
  allFound =
    allFound && eg->showerShapeValue(DeltaE, xAOD::EgammaParameters::DeltaE);
  // total shower width in 1st sampling
  allFound =
    allFound && eg->showerShapeValue(wtot, xAOD::EgammaParameters::wtots1);
  // E(+/-3)-E(+/-1)/E(+/-1)
  allFound =
    allFound && eg->showerShapeValue(fracm, xAOD::EgammaParameters::fracs1);

  if (m_useF3core) {
    allFound =
      allFound && eg->showerShapeValue(f3, xAOD::EgammaParameters::f3core);
  } else {
    allFound = allFound && eg->showerShapeValue(f3, xAOD::EgammaParameters::f3);
  }

  if (!allFound) {
    // if object is bad then use the bit for "bad eta"
    ATH_MSG_WARNING("Have some variables missing.");
    iflag = (0x1 << egammaPID::ClusterEtaRange_Electron);
    return iflag;
  }

  // For cut-based triggers above 20 GeV threshold, the online cut values on the
  // discriminant variables are always taken from the 20 GeV optimisation.
  //  if(et > 20000 )  { if(trigEtTh > 0) et = trigEtTh*1.01; }

  return m_rootTool->calocuts_electrons(eta2,
                                        et,
                                        Reta,  // replacing e233
                                        Rphi,  // replacing e237,
                                        Rhad1, // replacing ethad1,
                                        Rhad,  // replacing ethad,
                                        e277,
                                        weta1c,
                                        weta2c,
                                        f1,
                                        emax2,  // emax2
                                        Eratio, // emax
                                        DeltaE, // emin,
                                        wtot,
                                        fracm,
                                        f3,
                                        iflag,
                                        trigEtTh);
}

//================================================================
unsigned int
AsgElectronIsEMSelector::TrackCut(const xAOD::Electron* eg,
                                  float eta2,
                                  double et,
                                  double energy,
                                  unsigned int iflag) const
{
  // apply track cuts for electron identification
  //  - Track quality cuts
  //  - (eta,phi) and E/p matching between ID and ECAL
  //  - use of TRT
  // eg : egamma object
  // iflag: the starting isEM to use
  //
  // retrieve associated track
  const xAOD::TrackParticle* t = eg->trackParticle();

  // protection against bad pointers
  if (t == nullptr) {
    ATH_MSG_ERROR("Something is bad with the variables as passed");
    // if object is bad then use the bit for "bad eta"
    iflag = (0x1 << egammaPID::ClusterEtaRange_Electron);
    return iflag;
  }

  // Track quality cuts
  uint8_t nSiHitsPlusDeadSensors =
    ElectronSelectorHelpers::numberOfSiliconHitsAndDeadSensors(*t);
  uint8_t nPixHitsPlusDeadSensors =
    ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(*t);
  bool passBLayerRequirement =
    ElectronSelectorHelpers::passBLayerRequirement(*t);

  // TRT information
  uint8_t nTRThigh = 0;
  uint8_t nTRThighOutliers = 0;
  uint8_t nTRT = 0;
  uint8_t nTRTOutliers = 0;
  uint8_t nTRTXenonHits = 0;
  float TRT_PID = 0.0;

  bool allFound = true;

  allFound =
    allFound && t->summaryValue(nTRThigh, xAOD::numberOfTRTHighThresholdHits);
  allFound =
    allFound &&
    t->summaryValue(nTRThighOutliers, xAOD::numberOfTRTHighThresholdOutliers);
  allFound = allFound && t->summaryValue(nTRT, xAOD::numberOfTRTHits);
  allFound =
    allFound && t->summaryValue(nTRTOutliers, xAOD::numberOfTRTOutliers);
  allFound =
    allFound && t->summaryValue(nTRTXenonHits, xAOD::numberOfTRTXenonHits);
  allFound = allFound && t->summaryValue(TRT_PID, xAOD::eProbabilityHT);

  const float trackd0 = fabsf(t->d0());

  // Delta eta,phi matching
  float deltaeta;
  float deltaphi;

  allFound = allFound && eg->trackCaloMatchValue(
                           deltaeta, xAOD::EgammaParameters::deltaEta1);
  allFound = allFound && eg->trackCaloMatchValue(
                           deltaphi, xAOD::EgammaParameters::deltaPhi2);

  // E/p
  const double ep = energy * fabs(t->qOverP());

  if (!allFound) {
    // if object is bad then use the bit for "bad eta"
    ATH_MSG_WARNING("Have some variables missing.");
    iflag = (0x1 << egammaPID::ClusterEtaRange_Electron);
    return iflag;
  }

  return m_rootTool->TrackCut(eta2,
                              et,
                              passBLayerRequirement,
                              nPixHitsPlusDeadSensors,
                              nSiHitsPlusDeadSensors,
                              nTRThigh,
                              nTRThighOutliers,
                              nTRT,
                              nTRTOutliers,
                              nTRTXenonHits,
                              TRT_PID,
                              trackd0,
                              deltaeta,
                              deltaphi,
                              ep,
                              iflag);
}

