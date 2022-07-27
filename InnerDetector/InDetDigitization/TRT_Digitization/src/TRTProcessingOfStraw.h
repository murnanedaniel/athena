/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRT_DIGITIZATION_TRTPROCESSINGOFSTRAW_H
#define TRT_DIGITIZATION_TRTPROCESSINGOFSTRAW_H

#include "AthenaBaseComps/AthMessaging.h"

//Hit classes
#include "HitManagement/TimedHitCollection.h"
//Particle Table
#include "HepPDT/ParticleDataTable.hh"

#include "InDetIdentifier/TRT_ID.h"

#include "InDetSimEvent/TRTUncompressedHit.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MagField cache
#include "MagFieldElements/AtlasFieldCache.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TRT_ConditionsServices/ITRT_StrawStatusSummaryTool.h"
#include "TRT_ConditionsServices/ITRT_CalDbTool.h"

#include "TRTElectronicsProcessing.h"

#include "GeoPrimitives/GeoPrimitives.h"

#include "CLHEP/Random/RandomEngine.h"
#include "AtlasCLHEP_RandomGenerators/RandBinomialFixedP.h"

#include <memory>
#include <vector>

class TRTDigit;
class TRTTimeCorrection;
class TRTNoise;
class TRTDigCondBase;

class TRTUncompressedHit;
class ITRT_PAITool;
class ITRT_SimDriftTimeTool;

namespace InDetDD { class TRT_DetectorManager; }

class TRTDigSettings;

/**
 * TRT Digitization: Processing of a TRT Straws. @n
 * The main controlling function
 * in @c ProcessStraw(). See detailed description of this.
 */
class TRTProcessingOfStraw : public AthMessaging {
public:
  /** Constructor: Calls Initialize method */
  TRTProcessingOfStraw( const TRTDigSettings*,
                        const InDetDD::TRT_DetectorManager*,
                        ITRT_PAITool*,
                        ITRT_SimDriftTimeTool*,
                        TRTElectronicsProcessing * ep,
                        TRTNoise * noise,
                        TRTDigCondBase* digcond,
                        const HepPDT::ParticleDataTable*,
                        const TRT_ID*,
                        ITRT_PAITool* = NULL,
                        ITRT_PAITool* = NULL,
                        const ITRT_CalDbTool* = NULL);
  /** Destructor */
  ~TRTProcessingOfStraw();

  typedef TimedHitCollection<TRTUncompressedHit>::const_iterator
  hitCollConstIter;

  /**
   * Process this straw all the way from Geant4 @e hit to output @e digit.@n
   * Steps:
   * -# Loop over the simhits in this straw
   *    and produce a list of primary ionisation clusters.
   * -# Use the cluster list along with
   *    gas and wire properties, to create a list of energy
   *    deposits (i.e. potential fluctuations) reaching the
   *    frontend electronics.
   * -# Simulate how the FE
   *    turns the results into an output digit. This
   *    includes the shaping/amplification and subsequent
   *    discrimination as well as addition of noise.
   *
   * @param i:        Geant4 hit collection iterator
   * @param e:        last hit in collection
   * @param outdigit: The 27 bit digit
   * (bits: 8 low + 1 high + 8 low + 1 high + 8 low + 1 high)
   */
  void ProcessStraw (MagField::AtlasFieldCache& fieldCache,
                     hitCollConstIter i,
		     hitCollConstIter e,
		     TRTDigit& outdigit,
		     bool & m_alreadyPrintedPDGcodeWarning,
		     double m_cosmicEventPhase, //const ComTime* m_ComTime,
                     int strawGasType,
                     bool emulationArflag,
                     bool emulationKrflag,
                     CLHEP::HepRandomEngine* rndmEngine,
                     CLHEP::HepRandomEngine* elecProcRndmEngine,
                     CLHEP::HepRandomEngine* elecNoiseRndmEngine,
                     CLHEP::HepRandomEngine* paiRndmEngine );

private:

  //NB copy-constructor and assignment operator declared, but not defined.
  TRTProcessingOfStraw(const TRTProcessingOfStraw&);
  TRTProcessingOfStraw& operator= (const TRTProcessingOfStraw&);

  /** Initialize */
  void Initialize(const ITRT_CalDbTool *);

  const TRTDigSettings* m_settings;
  const InDetDD::TRT_DetectorManager* m_detmgr;
  ITRT_PAITool* m_pPAItoolXe;
  ITRT_PAITool* m_pPAItoolAr;
  ITRT_PAITool* m_pPAItoolKr;
  ITRT_SimDriftTimeTool* m_pSimDriftTimeTool;

  /** Time to be corrected for flight and wire propagation delays false when beamType='cosmics' */
  bool m_timeCorrection = false;

  double m_signalPropagationSpeed = 0.0;
  double m_attenuationLength = 0.0;

  bool m_useAttenuation = false;
  bool m_useMagneticFieldMap = false;

  double m_maxCrossingTime = 0.0;
  double m_minCrossingTime = 0.0;
  double m_shiftOfZeroPoint = 0.0;

  double m_innerRadiusOfStraw = 0.0;
  double m_outerRadiusOfWire = 0.0;

  double m_solenoidFieldStrength = 0.0;


  TRTTimeCorrection*        m_pTimeCorrection;
  TRTElectronicsProcessing* m_pElectronicsProcessing;
  TRTNoise*                 m_pNoise;
  TRTDigCondBase*           m_pDigConditions;

  const HepPDT::ParticleDataTable* m_pParticleTable;

  /** Primary ionisation cluster */
  class cluster {
  public:
    cluster(double e, double t, double x, double y, double z) :
      energy(e), time(t), xpos(x), ypos(y), zpos(z) {};
    double energy;
    double time;
    double xpos;
    double ypos;
    double zpos;
  };

  std::vector<cluster> m_clusterlist;
  std::vector<TRTElectronicsProcessing::Deposit> m_depositList;

  /**
   * This is the main function for re-simulation of the ionisation
   * in the active gas via the
   * PAI model. From the given G4 step, ionisation clusters are distributed
   * randomly (mean free path given by @c TRT_PAI_Process::GetMeanFreePath())
   * along path.
   * Cluster energies are given by @c TRT_PAI_Process::GetEnergyTransfer().
   * @param scaledKineticEnergy: The kinetic energy a proton would have had
   *                             if it had the same Lorentz gamma factor as
   *                             the particle in question.
   * @param particleCharge:      Particle charge
   * @param timeOfHit:           Time of hit
   * @param prex:                PreStepPoint @a x coordinate
   * @param prey:                PreStepPoint @a y coordinate
   * @param prez:                PreStepPoint @a z coordinate
   * @param postx:               PostStepPoint @a x coordinate
   * @param posty:               PostStepPoint @a y coordinate
   * @param postz:               PostStepPoint @a z coordinate
   * @param clusterlist:         List of ionisation clusters along step
   */
  void addClustersFromStep ( const double& scaledKineticEnergy,
                             const double& particleCharge,
                             const double& timeOfHit,
                             const double& prex,
                             const double& prey,
                             const double& prez,
                             const double& postx,
                             const double& posty,
                             const double& postz,
                             std::vector<cluster>& clusterlist,
                             int strawGasType,
                             CLHEP::HepRandomEngine* rndmEngine,
                             CLHEP::HepRandomEngine* paiRndmEngine);
  /**
   * Transform the ioniation clusters along the particle trajectory inside a
   * straw to energy deposits (i.e. potential fluctuations) reaching the
   * front-end electronics.
   * Effects taken into effect are:
   * - number of primary electrons from deposited energy
   * - stocastic recapture of drift electrons
   * - electron drift in the magnetic field (different direction
   *   end-cap/barrel). Optionally the drift time includes spread.
   * - Optionally signal wire propagation dealys are included.
   *   Each drift electron then gives rise to two signal: direct and reflected.
   *
   * @param hitID: Id of straw
   * @param clusters: ionisation clusters along particle trajectory
   * @param deposits: energy deposits on wire
   */
  void ClustersToDeposits (MagField::AtlasFieldCache& fieldCache, const int& hitID,
			   const std::vector<cluster>& clusters,
			   std::vector<TRTElectronicsProcessing::Deposit>& deposits,
			   Amg::Vector3D TRThitGlobalPos,
                           double m_cosmicEventPhase, // const ComTime* m_ComTime
                           int strawGasType,
                           CLHEP::HepRandomEngine* rndmEngine);

  double setClusterZ(double cluster_z_in, bool isLong, bool isShort, bool isEC) const;

  std::vector<double> m_drifttimes;     // electron drift times
  std::vector<double> m_expattenuation; // tabulation of exp()
  unsigned int  m_maxelectrons = 0U;         // maximum number of them (minmum is 100 for the Gaussian approx to be ok);

  bool m_alreadywarnedagainstpdg0;

  Amg::Vector3D getGlobalPosition( int hitID, const TimedHitPtr<TRTUncompressedHit> *theHit );

  std::unique_ptr<CLHEP::RandBinomialFixedP> m_randBinomialXe{};
  std::unique_ptr<CLHEP::RandBinomialFixedP> m_randBinomialKr{};
  std::unique_ptr<CLHEP::RandBinomialFixedP> m_randBinomialAr{};

protected:
  const TRT_ID* m_id_helper;

};

#endif
