// This file's extension implies that it's C, but it's really -*- C++ -*-.
/*
 * Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration.
 */
/**
 * @file LumiBlockComps/src/LuminosityCondAlg.h
 * @author scott snyder <snyder@bnl.gov>, from existing LuminosityTool.
 * @date May, 2019
 * @brief Conditions algorithm for luminosity data.
 */


#ifndef LUMIBLOCKCOMPS_LUMINOSITYCONDALG_H
#define LUMIBLOCKCOMPS_LUMINOSITYCONDALG_H


#include "LumiBlockData/LuminosityCondData.h"
#include "CoolLumiUtilities/OnlineLumiCalibrationCondData.h"
#include "CoolLumiUtilities/FillParamsCondData.h"
#include "CoolLumiUtilities/BunchLumisCondData.h"
#include "CoolLumiUtilities/BunchGroupCondData.h"
#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "CoralBase/Blob.h"
#include "xAODEventInfo/EventInfo.h"


/**
 * @brief Conditions algorithm for luminosity data.
 */
class LuminosityCondAlg
  : public AthReentrantAlgorithm
{
public:
  /// Forward base class ctor.
  using AthReentrantAlgorithm::AthReentrantAlgorithm;


  /// Gaudi initialize method.
  virtual StatusCode initialize() override;


  /// Algorithm execute method.
  virtual StatusCode execute (const EventContext& ctx) const override;
  virtual bool isReEntrant() const override final { return false; }


private:
  /**
   * @brief Unpack luminosity data from the attribute list.
   * @param lumiData Input luminosity data.
   * @param lumi Output luminosity data being filled.
   * @param preferredChannel[out] Preferred luminosity channel to use.
   * @param calibChannel[out] Calibration luminosity channel to use.
   * @param bunchInstLumiBlob[out] Packed per-bunch luminosity data.
   *                               Set to null for Run 1.
   *
   * Unpacks luminosity data from the attribute list.
   * Fills in the average luminosity fields in @c lumi,
   * and determines the luminosity channels to use.
   * For Run 2 and later, returns the packed luminosity data.
   */
  StatusCode
  updateAvgLumi (const CondAttrListCollection& lumiData,
                 LuminosityCondData& lumi,
                 unsigned int& preferredChannel,
                 unsigned int& calibChannel,
                 const coral::Blob*& bunchInstLumiBlob) const;


  /**
   * @brief Fill in per-bunch luminosity data.
   * @param ctx Event context.
   * @param bunchInstLumiBlob Packed per-bunch luminosity data.
   *                          Null for Run 1.
   * @param preferredChannel Preferred luminosity channel to use.
   * @param calibChannel Calibration luminosity channel to use.
   * @param wHdl WriteCondHandle conditions data being filled.
   *              Range is updated if needed.
   * @param lumi Output luminosity data being filled.
   */
  StatusCode
  updatePerBunchLumi (const EventContext& ctx,
                      const coral::Blob* bunchInstLumiBlob,
                      unsigned int preferredChannel,
                      unsigned int calibChannel,
                      SG::WriteCondHandle<LuminosityCondData>& wHdl,
                      LuminosityCondData& lumi) const;


  /**
   * @brief Fill in mu-to-lumi calibration.
   * @param ctx Event context.
   * @param calibChannel Calibration luminosity channel to use.
   * @param wHdl WriteCondHandle of the conditions data being filled.
   *              Range is updated if needed.
   * @param lumi Output luminosity data being filled.
   * @param isValid Set to false if data are not valid.
   */
  StatusCode
  updateMuToLumi (const EventContext& ctx,
                  unsigned int calibChannel,
                  SG::WriteCondHandle<LuminosityCondData>& wHdl,
                  LuminosityCondData& lumi,
                  bool& isValid) const;


  /**
   * @brief Fill in per-bunch luminosity data, run 2 and later.
   * @param bunchInstLumiBlob Packed per-bunch luminosity data.
   * @param preferredChannel Preferred luminosity channel to use.
   * @param lumi Output luminosity data being filled.
   */
  StatusCode
  updatePerBunchLumiRun2 (const coral::Blob& bunchInstLumiBlob,
                          unsigned int preferredChannel,
                          LuminosityCondData& lumi) const;


  /**
   * @brief Fill in per-bunch luminosity data, run 1.
   * @param preferredChannel Preferred luminosity channel to use.
   * @param wHdl WriteHandle of the conditions data being filled.
   *              Range is updated if needed.
   * @param lumi Output luminosity data being filled.
   */
  StatusCode
  updatePerBunchLumiRun1 (const EventContext& ctx,
                          unsigned int preferredChannel,
                          SG::WriteCondHandle<LuminosityCondData>& wHdl,
                          LuminosityCondData& lumi) const;


  Gaudi::Property<unsigned long> m_lumiChannel
    { this, "LumiChannelNumber", 0, "Luminosity channel to read.  0 means to determine from the data." };

  Gaudi::Property<unsigned long> m_calibBackupChannel
    { this, "CalibBackupChannel", 112, "Backup channel in case calibChannel doesn't exist in online calibration folder" };

  Gaudi::Property<bool> m_skipInvalid
    { this, "SkipInvalid", true, "Flag to control whether invalid data is skipped:  True (default), returning a zero luminosity; false, returning available luminosity values anyway." };

  Gaudi::Property<bool> m_expectInvalid
    { this, "ExpectInvalid", false, "Flag to control printouts when invalid data are encountered: True suppresses messages, False (default) leaves them in" };

  // following properties are for MC
  Gaudi::Property<bool> m_isMC
    { this, "IsMC", false, "Set to true when running on MC instead of data" };

  Gaudi::Property<float> m_muToLumi
    { this, "MCMuToLumi", 0.140569, "mu to lumi conversion factor in MC (80 mb/LHC rev freq)" };


  SG::ReadCondHandleKey<CondAttrListCollection> m_luminosityFolderInputKey
  { this, "LuminosityFolderInputKey", "/TRIGGER/OFLLUMI/LBLESTOFL",
    "Input luminosity COOL folder." };

  SG::ReadCondHandleKey<OnlineLumiCalibrationCondData> m_onlineLumiCalibrationInputKey
  { this, "OnlineLumiCalibrationInputKey", "OnlineLumiCalibrationCondData",
    "Input luminosity calibration." };

  SG::ReadCondHandleKey<BunchLumisCondData> m_bunchLumisInputKey
  { this, "BunchLumisInputKey", "",
    "Input raw luminosities.  Only used for Run 1." };

  SG::ReadCondHandleKey<BunchGroupCondData> m_bunchGroupInputKey
  { this, "BunchGroupInputKey", "",
    "Input filled bunch data.  Only used for Run 1." };

  SG::ReadCondHandleKey<FillParamsCondData> m_fillParamsInputKey
  { this, "FillParamsInputKey", "",
    "Input luminous bunch data.  Only used for Run 1." };

  // the following are for MC
  SG::ReadCondHandleKey<AthenaAttributeList> m_mcDigitizationInputKey
  { this, "DigitizationFolderInputKey", "/Digitization/Parameters",
    "Digitization parameters metadata folder." };

  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey 
  { this, "EventInfoKey", "EventInfo", "EventInfo key, used to read in simulated mu in MC" };

  SG::ReadDecorHandleKey<xAOD::EventInfo> m_actualMuKey 
  { this, "actualMuKey", "EventInfo.actualInteractionsPerCrossing",
    "Decoration for Actual Interaction Per Crossing, for MC" };

  SG::ReadDecorHandleKey<xAOD::EventInfo> m_averageMuKey 
  { this, "averageMuKey", "EventInfo.averageInteractionsPerCrossing",
    "Decoration for Average Interaction Per Crossing" };

  /// Output conditions object.
  SG::WriteCondHandleKey<LuminosityCondData> m_luminosityOutputKey
  { this, "LuminosityOutputKey", "LuminosityCondData",
    "Output luminosity data." };
};


#endif // not LUMIBLOCKCOMPS_LUMINOSITYCONDALG_H
