/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/MM_DigitContainerCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2016
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/MM_DigitContainerCnv_p1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "TestTools/initGaudi.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


#include "make_dd.icc"


void compare (const MmDigit& p1,
              const MmDigit& p2)
{
  assert (p1.identify() == p2.identify());
  assert (p1.stripResponseTime() == p2.stripResponseTime());
  assert (p1.stripResponseCharge() == p2.stripResponseCharge());
  assert (p1.stripResponsePosition() == p2.stripResponsePosition());
  assert (p1.chipResponseTime() == p2.chipResponseTime());
  assert (p1.chipResponsePosition() == p2.chipResponsePosition());
  assert (p2.stripTimeForTrigger().size() == 1);
  assert (p1.stripTimeForTrigger()[0] == p2.stripTimeForTrigger()[0]);
  assert (p2.stripPositionForTrigger().size() == 1);
  assert (p1.stripPositionForTrigger()[0] == p2.stripPositionForTrigger()[0]);
  assert (p2.stripChargeForTrigger().size() == 1);
  assert (p1.stripChargeForTrigger()[0] == p2.stripChargeForTrigger()[0]);
  assert (p2.MMFE_VMM_idForTrigger().size() == 1);
  assert (p1.MMFE_VMM_idForTrigger()[0] == p2.MMFE_VMM_idForTrigger()[0]);
  assert (p2.VMM_idForTrigger().size() == 1);
  assert (p1.VMM_idForTrigger()[0] == p2.VMM_idForTrigger()[0]);
}

void compare (const MmDigitContainer& p1,
              const MmDigitContainer& p2)
{
  MmDigitContainer::const_iterator it1 = p1.begin();
  MmDigitContainer::const_iterator it1e = p1.end();
  MmDigitContainer::const_iterator it2 = p2.begin();
  MmDigitContainer::const_iterator it2e = p2.end();
  while (it1 != it1e && it2 != it2e) {
    assert (it1.hashId() == it2.hashId());
    const MmDigitCollection& coll1 = **it1;
    const MmDigitCollection& coll2 = **it2;
    assert (coll1.size() == coll2.size());
    for (size_t j = 0; j < coll1.size(); j++) {
      compare (*coll1[j], *coll2[j]);
    }
    ++it1;
    ++it2;
  }
  assert (it1 == it1e && it2 == it2e);
}


void testit (const MmDigitContainer& trans1)
{
  MsgStream log (nullptr, "test");
  Muon::MM_DigitContainerCnv_p1 cnv;
  Muon::MM_DigitContainer_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  auto trans2 = std::unique_ptr<MmDigitContainer> (cnv.createTransient (&pers, log));

  compare (trans1, *trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE (const MuonGM::MuonDetectorManager& muo_dd)
{
  std::cout << "test1\n";

  MmDigitContainer trans1 (5);
  for (int hash=2; hash <= 3; hash++) {
    auto coll = std::make_unique<MmDigitCollection>
      (muo_dd.mmIdHelper()->elementID (55, 1, hash), // MML
       IdentifierHash(hash));

    for (int i=0; i < 10; i++) {
      int offs = i*10 + hash*100;
      Identifier id = muo_dd.mmIdHelper()->channelID (55, 1, hash,
                                                      1, 1, 2+i);
      auto dig = std::make_unique<MmDigit>
        (id,
         std::vector<float> {1.5f+offs, 2.5f+offs},
         std::vector<int> {3+offs, 4+offs, 5+offs},
         std::vector<float> {5.5f+offs, 6.5f+offs,  7.5f+offs},
         std::vector<float> {8.5f+offs, 9.5f+offs, 10.5f+offs, 11.5f+offs},
         std::vector<int> {12+offs, 13+offs},
         std::vector<float> {14.5f+offs},
         std::vector<float> {15.5f+offs, 16.5f+offs, 17.5f+offs, 18.5f+offs, 19.5f+offs},
         std::vector<int> {20+offs, 21+offs, 22+offs, 23+offs},
         std::vector<float> {24.5f+offs, 26.5f+offs},
         std::vector<int> {27+offs},
         std::vector<int> {28+offs, 29+offs}
         );
      coll->push_back (std::move (dig));
    }
    assert(trans1.addCollection (coll.release(), hash));
  }

  testit (trans1);
  Athena_test::Leakcheck check;
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  ISvcLocator* pSvcLoc;
  if (!Athena_test::initGaudi("MuonEventTPCnv/MuonEventTPCnv_test.txt", pSvcLoc))
  {
    std::cerr << "This test can not be run" << std::endl;
    return 0;
  }

  const MuonGM::MuonDetectorManager& muo_dd = make_dd();

  test1 (muo_dd);
  return 0;
}
