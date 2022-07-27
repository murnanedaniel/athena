/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/MM_DigitCnv_p2_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2016
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/MM_DigitCnv_p2.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


void compare (const MmDigit& p1,
              const MmDigit& p2)
{
  assert (p1.identify() == p2.identify());
  assert (p1.stripResponseTime() == p2.stripResponseTime());
  assert (p1.stripResponseCharge() == p2.stripResponseCharge());
  assert (p1.stripResponsePosition() == p2.stripResponsePosition());
  assert (p1.chipResponseTime() == p2.chipResponseTime());
  assert (p1.chipResponsePosition() == p2.chipResponsePosition());
  assert (p1.stripTimeForTrigger() == p2.stripTimeForTrigger());
  assert (p1.stripPositionForTrigger() == p2.stripPositionForTrigger());
  assert (p1.stripChargeForTrigger() == p2.stripChargeForTrigger());
  assert (p1.MMFE_VMM_idForTrigger() == p2.MMFE_VMM_idForTrigger());
  assert (p1.VMM_idForTrigger() == p2.VMM_idForTrigger());
}

void testit (const MmDigit& trans1)
{
  MsgStream log (nullptr, "test");
  MM_DigitCnv_p2 cnv;
  Muon::MM_Digit_p2 pers;
  cnv.transToPers (&trans1, &pers, log);
  MmDigit trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE ()
{
  std::cout << "test1\n";
  Athena_test::Leakcheck check;

  MmDigit trans1 (Identifier (1234),
                  std::vector<float> {1.5, 2.5},
                  std::vector<int> {3, 4, 5},
                  std::vector<float> {5.5, 6.5, 7.5},
                  std::vector<float> {8.5, 9.5, 10.5, 11.5},
                  std::vector<int> {12, 13},
                  std::vector<float> {14.5},
                  std::vector<float> {15.5, 16.5, 17.5, 18.5, 19.5},
                  std::vector<int> {20, 21, 22, 23},
                  std::vector<float> {24.5, 26.5},
                  std::vector<int> {27},
                  std::vector<int> {28, 29}
                  );
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  test1();
  return 0;
}
