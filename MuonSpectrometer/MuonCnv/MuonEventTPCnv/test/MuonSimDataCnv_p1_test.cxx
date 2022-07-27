/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file MuonEventTPCnv/test/MuonSimDataCnv_p1_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Dec, 2015
 * @brief Regression tests.
 */

#undef NDEBUG
#include "MuonEventTPCnv/MuonDigitContainer/MuonSimDataCnv_p1.h"
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


void compare (const MuonMCData& p1,
              const MuonMCData& p2)
{
  assert (p1.firstEntry() == p2.firstEntry());
  assert (p1.secondEntry() == p2.secondEntry());
}


void compare (const MuonSimData& p1,
              const MuonSimData& p2)
{
  assert (p1.word() == p2.word());
  assert (p1.globalPosition() == p2.globalPosition());
  const std::vector< MuonSimData::Deposit >& dep1 = p1.getdeposits();
  const std::vector< MuonSimData::Deposit >& dep2 = p2.getdeposits();
  assert (dep1.size() == dep2.size());
  for (size_t i = 0; i < dep1.size(); i++) {
    compare (dep1[i].first, dep2[i].first);
    assert (dep1[i].first == dep2[i].first);
    compare (dep1[i].second, dep2[i].second);
  }
}


void testit (const MuonSimData& trans1)
{
  MsgStream log (nullptr, "test");
  MuonSimDataCnv_p1 cnv;
  Muon::MuonSimData_p1 pers;
  cnv.transToPers (&trans1, &pers, log);
  MuonSimData trans2;
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

  std::vector<MuonSimData::Deposit> deps;
  HepMcParticleLink trkLink1(HepMC::barcode(genPartVector.at(0)),genPartVector.at(0)->parent_event()->event_number());
  deps.emplace_back (trkLink1, MuonMCData ( 2.5,  3.5));
  HepMcParticleLink trkLink2(HepMC::barcode(genPartVector.at(1)),genPartVector.at(1)->parent_event()->event_number());
  deps.emplace_back (trkLink2, MuonMCData (12.5, 13.5));
  HepMcParticleLink trkLink3(HepMC::barcode(genPartVector.at(2)),genPartVector.at(2)->parent_event()->event_number());
  deps.emplace_back (trkLink3, MuonMCData (22.5, 23.5));
  MuonSimData trans1 (deps, 4321);
  trans1.setPosition (Amg::Vector3D ( 4.5,  5.5,  6.5));
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
