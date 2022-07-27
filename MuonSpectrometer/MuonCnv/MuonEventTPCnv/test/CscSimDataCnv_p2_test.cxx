/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/CscSimDataCnv_p2_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Dec, 2015
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/CscSimDataCnv_p2.h"
#include "CxxUtils/checker_macros.h"
#include "TestTools/leakcheck.h"
#include "GaudiKernel/MsgStream.h"
#include <cassert>
#include <iostream>


#include "GeneratorObjectsTPCnv/initMcEventCollection.h"
#include "AtlasHepMC/GenEvent.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/Operators.h"


void compare (const HepMcParticleLink& p1,
              const HepMcParticleLink& p2)
{
  assert ( p1.isValid() == p2.isValid() );
  assert ( HepMC::barcode(p1) == HepMC::barcode(p2) );
  assert ( p1.eventIndex() == p2.eventIndex() );
  assert ( p1.getEventCollectionAsChar() == p2.getEventCollectionAsChar() );
  assert ( p1.cptr() == p2.cptr() );
  assert ( p1 == p2 );
}


void compare (const CscMcData& p1,
              const CscMcData& p2)
{
  assert (p1.energy() == p2.energy());
  assert (p1.ypos() == p2.ypos());
  assert (p1.zpos() == p2.zpos());
  assert (p1.charge() == p2.charge());
}


void compare (const CscSimData& p1,
              const CscSimData& p2)
{
  assert (p1.word() == p2.word());
  const std::vector< CscSimData::Deposit >& dep1 = p1.getdeposits();
  const std::vector< CscSimData::Deposit >& dep2 = p2.getdeposits();
  assert (dep1.size() == dep2.size());
  for (size_t i = 0; i < dep1.size(); i++) {
    compare (dep1[i].first, dep2[i].first);
    assert (dep1[i].first == dep2[i].first);
    compare (dep1[i].second, dep2[i].second);
  }
}


void testit (const CscSimData& trans1)
{
  MsgStream log (nullptr, "test");
  CscSimDataCnv_p2 cnv;
  Muon::CscSimData_p2 pers;
  cnv.transToPers (&trans1, &pers, log);
  CscSimData trans2;
  cnv.persToTrans (&pers, &trans2, log);

  compare (trans1, trans2);
}


void test1 ATLAS_NOT_THREAD_SAFE (std::vector<HepMC::GenParticlePtr>& genPartVector)
{
  std::cout << "test1\n";
  auto particle = genPartVector.at(0);
  // Create HepMcParticleLink outside of leak check.
  HepMcParticleLink dummyHMPL(HepMC::barcode(particle),particle->parent_event()->event_number());
  assert(dummyHMPL.cptr()==particle);
  Athena_test::Leakcheck check;

  std::vector<CscSimData::Deposit> deps;
  HepMcParticleLink trkLink1(HepMC::barcode(genPartVector.at(0)),genPartVector.at(0)->parent_event()->event_number());
  deps.emplace_back (trkLink1, CscMcData ( 2.5,  3.5,  4.5));
  HepMcParticleLink trkLink2(HepMC::barcode(genPartVector.at(1)),genPartVector.at(1)->parent_event()->event_number());
  deps.emplace_back (trkLink2, CscMcData (12.5, 13.5, 14.5));
  HepMcParticleLink trkLink3(HepMC::barcode(genPartVector.at(2)),genPartVector.at(2)->parent_event()->event_number());
  deps.emplace_back (trkLink3, CscMcData (22.5, 23.5, 24.5));
  deps[0].second.setCharge ( 5.5);
  deps[1].second.setCharge (15.5);
  deps[2].second.setCharge (25.5);
  CscSimData trans1 (deps, 4321);
  testit (trans1);
}


int main ATLAS_NOT_THREAD_SAFE ()
{
  ISvcLocator* pSvcLoc = nullptr;
  std::vector<HepMC::GenParticlePtr> genPartVector;
  if (!Athena_test::initMcEventCollection(pSvcLoc,genPartVector)) {
    std::cerr << "This test can not be run" << std::endl;
    return 0;
  }

  test1(genPartVector);
  return 0;
}
