/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/TgcPrepDataCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2016
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonPrepRawData/TgcPrepDataCnv_p1.h"
#include "MuonEventTPCnv/TgcPrepDataContainerCnv_tlp1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


void compare (const Trk::PrepRawData& p1,
              const Trk::PrepRawData& p2)
{
  assert (p1.identify() == p2.identify());
  assert (p1.localPosition()[0] == p2.localPosition()[0]);
  assert (p1.localCovariance() == p2.localCovariance());
  assert (p1.rdoList() == p2.rdoList());
}


void compare (const Muon::MuonCluster& p1,
              const Muon::MuonCluster& p2)
{
  compare (static_cast<const Trk::PrepRawData&>(p1),
           static_cast<const Trk::PrepRawData&>(p2));
  //assert (p1.globalPosition() == p2.globalPosition());
}


void compare (const Muon::TgcPrepData& p1,
              const Muon::TgcPrepData& p2)
{
  compare (static_cast<const Muon::MuonCluster&>(p1),
           static_cast<const Muon::MuonCluster&>(p2));
  assert (p1.detectorElement() == p2.detectorElement());
  //assert (p1.getBcBitMap() == p2.getBcBitMap());
  //assert (p1.globalPosition() == p2.globalPosition());
}

void testit (const Muon::TgcPrepData& trans1)
{
  MsgStream log (nullptr, "test");
  TgcPrepDataCnv_p1 cnv;
  TgcPrepDataContainerCnv_tlp1 tlcnv;
  cnv.setTopConverter (&tlcnv, TPObjRef::typeID_t());
  Muon::TgcPrepData_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  Muon::TgcPrepData trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE ()
{
  std::cout << "test1\n";

  Amg::Vector2D locpos (2.5, 3.5);

  Amg::MatrixX cov(1,1);
  cov(0,0) = 101;

  std::vector<Identifier> rdoList { Identifier(5432),
                                    Identifier(5361),
                                    Identifier(6456) };

  Muon::TgcPrepData trans1 (Identifier (1234),
                            IdentifierHash (1234),
                            locpos,
                            rdoList,
                            cov,
                            nullptr,
                            123);
                            
  testit (trans1);
  Athena_test::Leakcheck check;
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  test1();
  return 0;
}
