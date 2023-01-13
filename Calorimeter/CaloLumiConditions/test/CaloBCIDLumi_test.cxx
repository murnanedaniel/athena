/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file CaloLumiConditions/CaloBCIDLumi_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Mar, 2020
 * @brief Tests for CaloBCIDLumi.
 */

#undef NDEBUG
#include "AthenaKernel/getMessageSvc.h"
#include "CaloLumiConditions/CaloBCIDLumi.h"
#include "CaloLumiConditions/CaloBCIDCoeffs.h"
#include "LArElecCalib/ILArOFC.h"
#include "LArElecCalib/ILArShape.h"
#include "LArElecCalib/ILArMinBiasAverage.h"
#include "LArIdentifier/LArOnlineID.h"
#include "LumiBlockData/BunchCrossingCondData.h"
#include "LumiBlockData/LuminosityCondData.h"
#include "IdDictParser/IdDictParser.h"
#include "Identifier/HWIdentifier.h"
#include "TestTools/random.h"
#include "TestTools/FLOATassert.h"
#include <vector>
#include <unordered_map>
#include <cassert>
#include <iostream>


static constexpr size_t NCELL = 1833;


//************************************************************************

class LArOnlineIDTest
{
public:
  LArOnlineIDTest();
  ~LArOnlineIDTest();

  const LArOnlineID&  onlineID() const { return m_helper; }


private:
  std::unique_ptr<IdDictParser> m_parser;
  LArOnlineID                   m_helper;
};


LArOnlineIDTest::LArOnlineIDTest()
  : m_parser (std::make_unique<IdDictParser>())
{
  m_parser->register_external_entity("LArCalorimeter", "IdDictLArCalorimeter_DC3-05.xml");
  IdDictMgr& idd = m_parser->parse("IdDictParser/ATLAS_IDS.xml");
  m_helper.set_quiet (true);
  if (m_helper.initialize_from_dictionary(idd) != 0) {
    std::abort();
  }
}


LArOnlineIDTest::~LArOnlineIDTest()
= default;


//************************************************************************


class TestData
{
public:
  TestData (const uint32_t seed_in,
            const size_t ndata,
            const std::vector<HWIdentifier>& hwids,
            float maxval = 1,
            float offset = 0);
  LArVectorProxy data (const HWIdentifier& id) const;

  std::unordered_map<HWIdentifier, std::vector<float> > m_data;
};


TestData::TestData (uint32_t seed_in,
                    const size_t ndata,
                    const std::vector<HWIdentifier>& hwids,
                    float maxval /*= 1*/,
                    float offset /*= 0*/)
{
  for (HWIdentifier id : hwids) {
    std::vector<float>& data = m_data[id];
    data.resize (ndata);
    uint32_t seed = id.get_identifier32().get_compact() + seed_in;
    for (size_t i = 0; i < ndata; ++i) {
      data[i] = Athena_test::randf_seed (seed, maxval) - offset;
    }
  }
}


LArVectorProxy TestData::data (const HWIdentifier& id) const
{
  auto it = m_data.find (id);
  if (it != m_data.end()) {
    return {it->second.data(), it->second.data()+it->second.size()};
  }
  return {};
}


//************************************************************************


class TestOFC : public ILArOFC
{
public:
  using OFCRef_t = LArVectorProxy;

  TestOFC (size_t nofc, const std::vector<HWIdentifier>& hwids)
    : m_data (1234, nofc, hwids)
  {
  }

  virtual OFCRef_t OFC_a (const HWIdentifier& id,
                          int,
                          int=0) const override
  {
    return m_data.data (id);
  }
  
  virtual OFCRef_t OFC_b(const HWIdentifier&,
                         int,
                         int=0) const override
  { std::abort(); }
  
  virtual float timeOffset(const HWIdentifier&, int) const override
  { std::abort(); }

  virtual unsigned nTimeBins(const HWIdentifier&, int) const override
  { std::abort(); }
 
  virtual float timeBinWidth(const HWIdentifier&, int) const override
  { std::abort(); }

  TestData m_data;
};


//************************************************************************


// I'm kind of pulling a fast one here.
// The actual shapes have both positive and negative entries,
// averaging to zero.  That means that we're summing a series
// where the terms have varying sign, and so the result will depend
// on the order in which the sums are done.  Thus, in general, we won't
// get exactly the same result if we calculate the offset using the
// original method or the vectorized method.  However, here we choose
// the entries to be entirely non-negative.  In that case,
// the order of summation won't matter, so we should get the same result.
class TestShape : public ILArShape
{
public:
  TestShape (const std::vector<HWIdentifier>& hwids)
    : m_data (4321, 32, hwids)
  {
  }

  virtual ShapeRef_t Shape   (const HWIdentifier& id,
                              int,
                              int = 0,
                              int = 0 )  const override
  {
    return m_data.data (id);
  }

  
  virtual ShapeRef_t ShapeDer(const HWIdentifier&,
                              int,
                              int = 0,
                              int = 0 )  const override
  { std::abort(); }

  TestData m_data;
};


//************************************************************************


class TestMinBiasAverage : public ILArMinBiasAverage
{
public:
  TestMinBiasAverage (const std::vector<HWIdentifier>& hwids)
    : m_data (563424, 1, hwids, 100, 10)
  {
  }

  virtual const float& minBiasAverage(const HWIdentifier& id)  const override
  {
    return *m_data.data(id).data();
  }

  TestData m_data;
};


//************************************************************************


std::vector<HWIdentifier> get_hwids (const LArOnlineID& online_id)
{
  std::vector<HWIdentifier> hwids;

  size_t nhec_left = 100;
  for (HWIdentifier hwid : online_id.channel_range()) {
    if (hwids.size() == NCELL) {
      break;
    }
    else if (online_id.isHECchannel (hwid) && nhec_left > 0) {
      hwids.push_back (hwid);
      --nhec_left;
    }
    else if ((NCELL - hwids.size()) > nhec_left) {
      hwids.push_back (hwid);
    }
  }
  
  return hwids;
}


//************************************************************************


class BunchCrossingCondAlg
{
public:
  static void fill (BunchCrossingCondData& bccd);
};


void BunchCrossingCondAlg::fill (BunchCrossingCondData& bccd)
{
  static constexpr int MAX_BCID = BunchCrossingCondData::m_MAX_BCID;
  uint32_t seed = 12734;
  for (size_t i = 0; i < MAX_BCID; i++) {
    bccd.m_luminous[i] = Athena_test::rng_seed (seed) & 0x10000;
  }
}


BunchCrossingCondData get_bccd()
{
  BunchCrossingCondData bccd;
  BunchCrossingCondAlg::fill (bccd);
  return bccd;
}


std::vector<float> get_lumivec()
{
  static constexpr int MAX_BCID = BunchCrossingCondData::m_MAX_BCID;
  uint32_t seed = 42345;
  std::vector<float> lumivec;
  lumivec.reserve (200);
  for (size_t i = 0; i < MAX_BCID; i++) {
    lumivec.push_back (Athena_test::randf_seed (seed, 10));
  }
  return lumivec;
}


//************************************************************************


std::vector<float> calc_old (const std::vector<HWIdentifier>& hwids,
                             const LArOnlineID& online_id,
                             const ILArOFC& ofcs,
                             const ILArShape& shapes,
                             const ILArMinBiasAverage& minbias,
                             const std::vector<float>& lumiVec,
                             const size_t bcid)
{
  static constexpr int MAX_BCID = BunchCrossingCondData::m_MAX_BCID;
  std::vector<float> out;
  out.reserve (hwids.size());
  for (HWIdentifier hwid : hwids) {
    float eOFC = 0;

    const auto& ofc = ofcs.OFC_a(hwid,0);
    const unsigned int nsamples = ofc.size();

    const auto& samp = shapes.Shape(hwid,0);
    const unsigned int nshapes = samp.size();

    // choice of first sample : i.e sample on the pulse shape to which first OFC sample is applied
    unsigned int ifirst= 0;
    if (nsamples==4) {
      if (online_id.isHECchannel(hwid)) {
        ifirst=1/*m_firstSampleHEC*/;
      }
    }
    unsigned int ishift = ifirst + bcid; // correct for first sample
    for (unsigned int i=0;i<nsamples;i++) {
      float sumShape=0.;
      /*unsigned*/ int k = ishift;
      if (k >= MAX_BCID) k-= MAX_BCID;
      for (unsigned int j=0;j<nshapes;j++) {
        const float& lumi = lumiVec[k];
        sumShape += samp[j]*lumi;
        k--;
        if (k < 0) k += MAX_BCID;
      } 
      eOFC += sumShape*(ofc[i]);
      ishift++;
    }

    float MinBiasAverage = minbias.minBiasAverage(hwid);
    if (MinBiasAverage<0.) MinBiasAverage=0.;
    eOFC = eOFC * MinBiasAverage;
    out.push_back (eOFC);
  }
  return out;
}


std::vector<float> calc_old (const std::vector<HWIdentifier>& hwids,
                             const LArOnlineID& online_id,
                             const ILArOFC& ofcs,
                             const ILArShape& shapes,
                             const ILArMinBiasAverage& minbias,
                             const BunchCrossingCondData& bcData,
                             const float averageInteractionsPerCrossing,
                             const size_t bcid)
{
  static constexpr int MAX_BCID = BunchCrossingCondData::m_MAX_BCID;
  std::vector<float> lumiVec;
  lumiVec.assign(MAX_BCID,0.0);
  const float xlumiMC = averageInteractionsPerCrossing*0.158478605;
  int ii = bcid-38;
  if (ii < 0) ii += MAX_BCID;
  for (int i=bcid-38; i<(int)bcid+38; ++i) {
    lumiVec[ii]=bcData.isFilled(ii)*xlumiMC;
    ++ii;
    if (ii >= MAX_BCID) ii -= MAX_BCID;
  }

  return calc_old (hwids,
                   online_id,
                   ofcs,
                   shapes,
                   minbias,
                   lumiVec,
                   bcid);
}


void test1 (const size_t nofc, const LArOnlineID& online_id)
{
  std::cout << "test1 " << nofc << "\n";

  std::vector<HWIdentifier> hwids = get_hwids (online_id);
  TestOFC ofcs (nofc, hwids);
  TestShape shapes (hwids);
  TestMinBiasAverage minbias (hwids);

  CaloBCIDCoeffs coeffs (hwids, online_id, ofcs, shapes, minbias);

  {
    BunchCrossingCondData bccd = get_bccd();
    CaloBCIDLumi bcidlumi (coeffs, bccd);
    CxxUtils::vec_aligned_vector<float> out;
    for (size_t bcid = 0; bcid < BunchCrossingCondData::m_MAX_BCID; ++bcid) {
      bcidlumi.calc (bcid, 2, out);

      std::vector<float> out_old = calc_old (hwids,
                                             online_id,
                                             ofcs,
                                             shapes,
                                             minbias,
                                             bccd,
                                             2,
                                             bcid);
      assert (out.size() == out_old.size());
      for (size_t i = 0; i < out.size(); i++) {
        assert (Athena_test::isEqual (out[i], out_old[i]));
      }
    }
  }

  {
    LuminosityCondData lcd;
    lcd.setLbLuminosityPerBCIDVector (get_lumivec());
    CaloBCIDLumi bcidlumi (coeffs, lcd);
    CxxUtils::vec_aligned_vector<float> out;
    for (size_t bcid = 0; bcid < BunchCrossingCondData::m_MAX_BCID; ++bcid) {
      bcidlumi.calc (bcid, 1, out);

      std::vector<float> out_old = calc_old (hwids,
                                             online_id,
                                             ofcs,
                                             shapes,
                                             minbias,
                                             lcd.lbLuminosityPerBCIDVector(),
                                             bcid);

      assert (out.size() == out_old.size());
      for (size_t i = 0; i < out.size(); i++) {
        assert (Athena_test::isEqual (out[i], out_old[i]));
      }
    }
  }
}


int main()
{
  std::cout << "CaloLumiConditions/CaloBCIDLumi_test\n";
  Athena::getMessageSvcQuiet = true;
  LArOnlineIDTest larhelpers;
  test1 (5, larhelpers.onlineID());
  test1 (4, larhelpers.onlineID());

  return 0;
}
