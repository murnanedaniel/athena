// This file's extension implies that it's C, but it's really -*- C++ -*-.

/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file DataModelTestDataRead/versions/G_v1.h
 * @file G_v1.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Nov, 2014
 * @brief Test for xAOD auto schema evolution.
 */


#ifndef DATAMODELTESTDATAREAD_G_V1_H
#define DATAMODELTESTDATAREAD_G_V1_H


#include "AthContainers/AuxElement.h"
#include "AthenaKernel/BaseInfo.h"


namespace DMTest {


class G_v1
  : public SG::AuxElement
{
public:
  int anInt() const;
  void setAnInt (int i);
  float gFloat() const;
  void setgFloat (float d);
  const std::vector<float>& gvFloat() const;
  void setgvFloat (const std::vector<float>& v);
};


} // namespace DMTest


SG_BASE (DMTest::G_v1, SG::AuxElement);


#endif // not DATAMODELTESTDATAREAD_G_V1_H
