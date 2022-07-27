/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file PixelConditionsAlgorithms/PixelDCSCondStatusAlg.h
 * @author Soshi Tsuno <Soshi.Tsuno@cern.ch>
 * @date November, 2019
 * @brief Created pixel DCS module status in PixelDCSStatusData.
 */

#ifndef PIXELDCSCONDSTATUSALG
#define PIXELDCSCONDSTATUSALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"

#include "StoreGate/WriteCondHandleKey.h"
#include "PixelConditionsData/PixelDCSStatusData.h"

#include "InDetIdentifier/PixelID.h"

#include "Gaudi/Property.h"

class PixelDCSCondStatusAlg : public AthReentrantAlgorithm {
  public:
    PixelDCSCondStatusAlg(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~PixelDCSCondStatusAlg() = default;

    virtual StatusCode initialize() override final;
    virtual StatusCode execute(const EventContext& ctx) const override final;
    virtual bool isReEntrant() const override final { return false; }

  private:
    const PixelID* m_pixelID{nullptr};

    SG::ReadCondHandleKey<CondAttrListCollection> m_readKeyStatus
    {this, "ReadKeyStatus", "", "Key of input DCS status conditions folder"};

    SG::WriteCondHandleKey<PixelDCSStatusData> m_writeKeyStatus
    {this, "WriteKeyStatus", "PixelDCSStatusCondData", "Key of output DCS status data"};

};

#endif
