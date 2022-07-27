/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/STGC_DigitContainerCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2016
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/STGC_DigitContainerCnv_p1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "TestTools/initGaudi.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


#include "make_dd.icc"


void compare (const sTgcDigit& p1,
              const sTgcDigit& p2)
{
  assert (p1.identify() == p2.identify());
  assert (p1.bcTag() == p2.bcTag());
  assert (p1.time() == p2.time());
  assert (p1.charge() == p2.charge());
}

void compare (const sTgcDigitContainer& p1,
              const sTgcDigitContainer& p2)
{
  sTgcDigitContainer::const_iterator it1 = p1.begin();
  sTgcDigitContainer::const_iterator it1e = p1.end();
  sTgcDigitContainer::const_iterator it2 = p2.begin();
  sTgcDigitContainer::const_iterator it2e = p2.end();
  while (it1 != it1e && it2 != it2e) {
    assert (it1.hashId() == it2.hashId());
    const sTgcDigitCollection& coll1 = **it1;
    const sTgcDigitCollection& coll2 = **it2;
    assert (coll1.size() == coll2.size());
    for (size_t j = 0; j < coll1.size(); j++) {
      compare (*coll1[j], *coll2[j]);
    }
    ++it1;
    ++it2;
  }
  assert (it1 == it1e && it2 == it2e);
}


void testit (const sTgcDigitContainer& trans1)
{
  MsgStream log (nullptr, "test");
  Muon::STGC_DigitContainerCnv_p1 cnv;
  Muon::STGC_DigitContainer_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  auto trans2 = std::unique_ptr<sTgcDigitContainer> (cnv.createTransient (&pers, log));

  compare (trans1, *trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE (const MuonGM::MuonDetectorManager& muo_dd)
{
  std::cout << "test1\n";

  sTgcDigitContainer trans1 (5);
  for (int hash=2; hash <= 3; hash++) {
    auto coll = std::make_unique<sTgcDigitCollection>
      (muo_dd.stgcIdHelper()->elementID (57, 1, hash), // STS
       IdentifierHash(hash));

    for (int i=0; i < 10; i++) {
      int offs = i*10 + hash*100;
      Identifier id = muo_dd.stgcIdHelper()->channelID (57, 1, hash,
                                                        1, 2, 1, 2+i);
      auto dig = std::make_unique<sTgcDigit>
        (id,
         123+offs, 4.5+offs, 6.5+offs, false, false);
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
