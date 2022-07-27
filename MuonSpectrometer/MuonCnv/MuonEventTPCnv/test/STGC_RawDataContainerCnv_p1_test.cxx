/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/STGC_RawDataContainerCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2016
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonRDO/STGC_RawDataContainerCnv_p1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "TestTools/initGaudi.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


#include "make_dd.icc"

using namespace Muon;

void compare (const Muon::STGC_RawData& p1,
              const Muon::STGC_RawData& p2)
{
  assert (p1.identify() == p2.identify());
}

void compare (const Muon::STGC_RawDataContainer& p1,
              const Muon::STGC_RawDataContainer& p2)
{
  STGC_RawDataContainer::const_iterator it1 = p1.begin();
  STGC_RawDataContainer::const_iterator it1e = p1.end();
  STGC_RawDataContainer::const_iterator it2 = p2.begin();
  STGC_RawDataContainer::const_iterator it2e = p2.end();
  unsigned int count=0;
  while (it1 != it1e && it2 != it2e) {
    std::cout <<"Collection #"<<++count<<std::endl;
    assert (it1.hashId() == it2.hashId());
    const STGC_RawDataCollection& coll1 = **it1;
    const STGC_RawDataCollection& coll2 = **it2;
    assert (coll1.size() == coll2.size());
    for (size_t j = 0; j < coll1.size(); j++) {
      compare (*coll1[j], *coll2[j]);
    }
    ++it1;
    ++it2;
  }
  assert (it1 == it1e && it2 == it2e);
}


void testit (const STGC_RawDataContainer& trans1, const MuonGM::MuonDetectorManager& muon_dd)
{
  MsgStream log (nullptr, "test");
  Muon::STGC_RawDataContainerCnv_p1 cnv;
  cnv.initialize(muon_dd.stgcIdHelper());
  Muon::STGC_RawDataContainer_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  auto trans2 = std::unique_ptr<STGC_RawDataContainer> (cnv.createTransient (&pers, log));
  std::cout<<" Before has size: "<<trans1.numberOfCollections()<<", after has size: "<<trans2->numberOfCollections()<<std::endl;
  compare (trans1, *trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE (const MuonGM::MuonDetectorManager& muo_dd)
{
  std::cout << "test1\n";

  STGC_RawDataContainer trans1 (5);
  for (int hash=2; hash <= 3; hash++) {
    auto coll = std::make_unique<STGC_RawDataCollection>(hash);
    // (muo_dd.stgcIdHelper()->elementID (1, 1, hash), // STS
    //  IdentifierHash(hash));

    // TODO This needs to be made reasonable. Ed
    for (int i=0; i < 10; i++) {
      // int offs = i*10 + hash*100;
      // Identifier id = muo_dd.stgcIdHelper()->channelID (1, 1, hash,
      //                                                   1, 2, 1, 2+i);
      auto dig = std::make_unique<STGC_RawData>
        ( Identifier(0x712a0054) );
      coll->push_back (std::move (dig));
    }
    assert(trans1.addCollection (coll.release(), hash));
  }

  testit (trans1, muo_dd);
  Athena_test::Leakcheck check;
  testit (trans1, muo_dd);
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

