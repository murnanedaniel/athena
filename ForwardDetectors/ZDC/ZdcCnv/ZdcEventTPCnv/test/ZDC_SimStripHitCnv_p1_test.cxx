/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file ZdcEventTPCnv/test/ZDC_SimStripHitCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Mar, 2016
 * @brief Tests for ZDC_SimStripHitCnv_p1.
 */


#undef NDEBUG
#include "ZdcEventTPCnv/ZDC_SimStripHitCnv_p1.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include <cassert>
#include <iostream>


void compare (const ZDC_SimStripHit& p1,
              const ZDC_SimStripHit& p2)
{
  assert (p1.GetSide() == p2.GetSide());
  assert (p1.GetMod() == p2.GetMod());
  assert (p1.GetEdep() == p2.GetEdep());
  assert (p1.GetNPhotons() == p2.GetNPhotons());
}


void testit (const ZDC_SimStripHit& trans1)
{
  MsgStream log (nullptr, "test");
  ZDC_SimStripHitCnv_p1 cnv;
  ZDC_SimStripHit_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  ZDC_SimStripHit trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE ()
{
  std::cout << "test1\n";
  Athena_test::Leakcheck check;

  ZDC_SimStripHit trans1 (123, 234, 21, 12345.5);
    
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  test1();
  return 0;
}
