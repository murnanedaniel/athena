/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/MuonMCDataCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Dec, 2015
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/MuonMCDataCnv_p1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


void compare (const MuonMCData& p1,
              const MuonMCData& p2)
{
  assert (p1.firstEntry() == p2.firstEntry());
  assert (p1.secondEntry() == p2.secondEntry());
}

void testit (const MuonMCData& trans1)
{
  MsgStream log (nullptr, "test");
  MuonMCDataCnv_p1 cnv;
  Muon::MuonMCData_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  MuonMCData trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE ()
{
  std::cout << "test1\n";
  Athena_test::Leakcheck check;

  MuonMCData trans1 (2.5, 3.5);
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  test1();
  return 0;
}
