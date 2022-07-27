/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file PixelConditionsAlgorithms/PixelDeadMapCondAlg.h
 * @author Soshi Tsuno <Soshi.Tsuno@cern.ch>
 * @date December, 2019
 * @brief Store pixel module/FE status in PixelDeadMapCondData.
 */

#ifndef PIXELDEADMAPCONDALG_H
#define PIXELDEADMAPCONDALG_H

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"

#include "StoreGate/WriteCondHandleKey.h"
#include "PixelConditionsData/PixelDeadMapCondData.h"

#include "Gaudi/Property.h"

class PixelDeadMapCondAlg : public AthReentrantAlgorithm {
  public:
    PixelDeadMapCondAlg(const std::string& name, ISvcLocator* pSvcLocator);

    virtual StatusCode initialize() override final;
    virtual StatusCode execute(const EventContext& ctx) const override final;
    virtual bool isReEntrant() const override final { return false; }

  private:
    SG::ReadCondHandleKey<CondAttrListCollection> m_readKey
    {this, "ReadKey", "/PIXEL/PixelModuleFeMask", "Input deadmap folder"};

    SG::WriteCondHandleKey<PixelDeadMapCondData> m_writeKey
    {this, "WriteKey", "PixelDeadMapCondData", "Output deadmap data"};

};

#endif
