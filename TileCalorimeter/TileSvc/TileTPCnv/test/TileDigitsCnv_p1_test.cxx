/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file TileTPCnv/test/TileDigitsCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2016
 * @brief Tests for TileDigitsCnv_p1.
 */


#undef NDEBUG
#include "TileTPCnv/TileDigitsCnv_p1.h"
#include "TileConditions/TileCablingService.h"
#include "TileIdentifier/TileHWID.h"
#include "IdDictParser/IdDictParser.h"
#include "TestTools/initGaudi.h"
#include "TestTools/leakcheck.h"
#include "CxxUtils/checker_macros.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


TileHWID hwid;
TileTBID tbid;
TileID   tileid;

class TileCablingSvc
{
public:
  static
  void init_idhelpers ATLAS_NOT_THREAD_SAFE (IdDictParser& parser)
  {
    tileid.set_do_neighbours (false);
    IdDictMgr& idd = parser.parse ("IdDictParser/ATLAS_IDS.xml");
    assert (hwid.initialize_from_dictionary (idd) == 0);
    assert (tbid.initialize_from_dictionary (idd) == 0);
    assert (tileid.initialize_from_dictionary (idd) == 0);
    TileCablingService* svc = TileCablingService::getInstance_nc();
    svc->setTileHWID (&hwid);
    svc->setTileTBID (&tbid);
    svc->setTileID (&tileid);
  }
};


void compare (const TileDigits& p1,
              const TileDigits& p2)
{
  assert (p1.identify() == p2.identify());
  assert (p1.samples() == p2.samples());
}


void testit (const TileDigits& trans1)
{
  MsgStream log (0, "test");
  TileDigitsCnv_p1 cnv;
  TileDigits_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  TileDigits trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE ()
{
  std::cout << "test1\n";
  Athena_test::Leakcheck check;

  std::vector<float> digits {4.5, 7.5, 10.5, 12.5, 3.5 };
  TileDigits trans1 (hwid.adc_id (1, 2, 3, 0), digits);

  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  IdDictParser parser;
  TileCablingSvc::init_idhelpers(parser);
  test1();
  return 0;
}
