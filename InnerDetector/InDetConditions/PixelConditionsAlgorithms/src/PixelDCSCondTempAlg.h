/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file PixelConditionsAlgorithms/PixelDCSCondTempAlg.h
 * @author Soshi Tsuno <Soshi.Tsuno@cern.ch>
 * @date November, 2019
 * @brief Created pixel temperature data in PixelDCSTempData.
 */

#ifndef PIXELDCSCONDTEMPALG
#define PIXELDCSCONDTEMPALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"

#include "StoreGate/WriteCondHandleKey.h"
#include "PixelConditionsData/PixelModuleData.h"
#include "PixelConditionsData/PixelDCSTempData.h"

#include "InDetIdentifier/PixelID.h"

#include "Gaudi/Property.h"

class PixelDCSCondTempAlg : public AthReentrantAlgorithm {
  public:
    PixelDCSCondTempAlg(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~PixelDCSCondTempAlg() = default;

    virtual StatusCode initialize() override final;
    virtual StatusCode execute(const EventContext& ctx) const override final;
    virtual bool isReEntrant() const override final { return false; }

  private:
    const PixelID* m_pixelID{nullptr};

    SG::ReadCondHandleKey<PixelModuleData> m_moduleDataKey
    {this, "PixelModuleData", "PixelModuleData", "Pixel module data"};

    SG::ReadCondHandleKey<CondAttrListCollection> m_readKey
    {this, "ReadKey", "/PIXEL/DCS/TEMPERATURE", "Key of input (raw) temperature conditions folder"};

    SG::WriteCondHandleKey<PixelDCSTempData> m_writeKey
    {this, "WriteKey", "PixelDCSTempCondData", "Key of output (derived) temperature conditions folder"};

};

#endif
