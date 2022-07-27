/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRT_DIGITIZATION_TRTELECTRONICSPROCESSING_H
#define TRT_DIGITIZATION_TRTELECTRONICSPROCESSING_H

#include "AthenaBaseComps/AthMessaging.h"
#include "GaudiKernel/ServiceHandle.h"

#include "CLHEP/Random/RandomEngine.h"

#include <atomic>
#include <vector>

class TRTDigSettings;
class TRTDigit;
class TRTElectronicsNoise;

/**
 * Electronics Processing
 */
class TRTElectronicsProcessing : public AthMessaging {
public:
  TRTElectronicsProcessing( const TRTDigSettings* digset,
                            TRTElectronicsNoise * electronicsnoise );
  ~TRTElectronicsProcessing();

  /**
   * TRT energy deposit
   */
  class Deposit {
  public:
    /**
     * Create deposit
     */
    Deposit(const double& e, const double& t): m_energy(e), m_time(t) {}

    double energy() const { return m_energy; } /**< Get energy of deposit */
    double time() const { return m_time; }     /**< Get time of deposit */


  private:
    double m_energy;  /**< Energy */
    double m_time;    /**< Time   */
  };

  /**
   * Process deposits in a straw
   * -# fill cluster energies into appropriate time bins
   * -# signal shaping
   * -# add analog electronics noise
   * -# discriminator response
   * -# encode into a single digit
   *
   * @param deposits:       all deposits, this straw
   * @param hitID:          ID of straw
   * @param outdigit:       output TRTDigit
   * @param lowthreshold:   discriminator value, low threshold
   * @param noiseamplitude: noise amplitude
   * @param highthreshold:  discriminator value, high threshold
   */
  void ProcessDeposits( const std::vector<Deposit>&,
                        const int& hitID,
                        TRTDigit& outdigit,
                        double lowthreshold,
                        const double& noiseamplitude,
                        int strawGasType,
                        CLHEP::HepRandomEngine* rndmEngine,
                        CLHEP::HepRandomEngine* elecNoiseRndmEngine,
                        double highthreshold = -1.0
                        );

  double getHighThreshold ( int hitID, int strawGasType );
  void LTt0Shift          ( int hitID, int strawGasType );
  void HTt0Shift          ( int hitID );

private:
  //NB copy-constructor and assignment operator declared, but not defined.
  TRTElectronicsProcessing(const TRTElectronicsProcessing&);
  TRTElectronicsProcessing& operator= (const TRTElectronicsProcessing&);

  //--- Initialization methods

  /**
   * Initialize constants. Call TabulateSignalShape()
   */
  void Initialize();

  /**
   * Tabulate low and high level signal shapes from functions
   * LowThresholdSignalShape() and HighThresholdSignalShape(), respectively
   */
  void TabulateSignalShape();

  //---

  /**
   * Shape electron drift signal according to appropriate signal shapes.
   * - input: @c m_energyDistribution[]
   * - @param strawGasType: straw type (1:Xe, 2:Kr, 3:Ar)
   * - output: @c m_lowThresholdSignal[], @c m_highThresholdSignal[]
   */
  void SignalShaping(int strawGasType);

  /**
   * Simulate discriminator response
   * - input: @c m_lowThresholdSignal[], @c m_highThresholdSignal[]
   * - output: @c m_lowThresholdDiscriminator[],
   *           @c m_highThresholdDiscriminator[],
   * @param lowthreshold:   low threshold discriminator setting
   * @param highthreshold: high threshold discriminator setting
   */
  void DiscriminatorResponse(const double& lowthreshold,
                             const double& highthreshold) const;

  /**
   * Encode 27 bit digit from discriminator response
   * [8 low + 1 high + 8 low + 1 high + 8 low + 1 high]:
   * - input: @c m_lowThresholdDiscriminator[],
   *          @c m_highThresholdDiscriminator[]
   *
   * @return digit
   */
  unsigned EncodeDigit() const;

  const TRTDigSettings* m_settings;

  TRTElectronicsNoise * m_pElectronicsNoise;

  double m_timeInterval = 0.0;       /**< Time interval covered by digit [75 ns] */

  int m_totalNumberOfBins = 0;     /**< Total no. internal bins             */
  int m_numberOfPreZeroBins = 0;   /**< No. internal bins before main event */
  int m_numberOfPostZeroBins = 0;  /**< No. internal bins after main event  */

  double m_binWidth = 0.0;   /**< Time interval corresponding to each internal bin*/

  int m_minDiscriminatorWidthInBinWidths = 0; /**< Min. discriminator width [internal bins]*/
  int m_discriminatorSettlingTimeInBinWidths = 0; /**< Discriminator settling time [bins]  */
  int m_discriminatorDeadTimeInBinWidths = 0;     /**< Discriminator dead time [int. bins] */
  int m_minWidthMinusSettlingTimeInBinWidths = 0; /**< Min. discriminator time minus settling * time [internal bins] */

  // 0:Xe,1:Kr,2:Ar
  double m_lowThresholdBar[3]{};
  double m_lowThresholdEC[3]{};
  std::vector<double> m_lowThresholdSignalShape[3];
  std::vector<double> m_highThresholdSignalShape[3];

  // unsigned int m_maskA;   /**< mask - ever used? */
  // unsigned int m_maskB;   /**< mask - ever used? */
  // unsigned int m_maskC;   /**< mask - ever used? */
  // unsigned int m_maskHT;  /**< mask - ever used? */

  // Deposit energy in timed bins before shaping.
  double* m_energyDistribution = nullptr;

  // LT and HT signals after shaping
  std::vector<double> m_lowThresholdSignal;
  std::vector<double> m_highThresholdSignal;

  // LT and HT discriminated signals (bits).
  int* m_lowThresholdDiscriminator = nullptr;  /**< Signal after discrimination */
  int* m_highThresholdDiscriminator = nullptr; /**< Signal after discrimination */

  mutable std::atomic<bool> m_first{true};
};

#endif
