/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonSimEvent/sTgcHitIdHelper.h"

#include <iostream>
#include <iomanip> // for std::array

sTgcHitIdHelper* sTgcHitIdHelper::m_help = nullptr;

namespace {
    const static std::array<char, 1> v1 = {'T'};
    const static std::array<char, 2> v2 = {'S','L'};
    const static std::array<char, 5> v3 = {'0','1','2','3','4'};
    const static std::array<char, 2> v4 = {'1','2'};
}

//private constructor
sTgcHitIdHelper::sTgcHitIdHelper() : HitIdHelper(){
  InitializeStationName();
  Initialize();
}

sTgcHitIdHelper* sTgcHitIdHelper::GetHelper(){
  if (!m_help) m_help = new sTgcHitIdHelper();
  return m_help;
}

void sTgcHitIdHelper::Initialize(){

  InitializeField("PhiSector",1,16);
  InitializeField("ZSector",0,4);
  InitializeField("MultiLayer",1,2);
  InitializeField("Layer",1,4);
  InitializeField("Side",-1,1);

}
void sTgcHitIdHelper::InitializeStationName()
{
  InitializeField("Station[1]",0,sizeof(v1));
  InitializeField("Station[2]",0,sizeof(v2));
  InitializeField("Station[3]",0,sizeof(v3));
  InitializeField("Station[4]",0,sizeof(v4));
}
void sTgcHitIdHelper::SetStationName(std::string name, int& hid) const
{
  for (unsigned int i=0;i<sizeof(v1);i++)
    if (v1[i]==name[0]) SetFieldValue("Station[1]",i,hid);
  for (unsigned int i=0;i<sizeof(v2);i++)
    if (v2[i]==name[1]) SetFieldValue("Station[2]",i,hid);
  for (unsigned int i=0;i<sizeof(v3);i++)
    if (v3[i]==name[2]) SetFieldValue("Station[3]",i,hid);
  for (unsigned int i=0;i<sizeof(v4);i++)
    if (v4[i]==name[3]) SetFieldValue("Station[4]",i,hid);
}
std::string sTgcHitIdHelper::GetStationName(const int& hid) const
{
  char v[5];
  v[0]=v1[this->GetFieldValue("Station[1]", hid)];
  v[1]=v2[this->GetFieldValue("Station[2]", hid)];
  v[2]=v3[this->GetFieldValue("Station[3]", hid)];
  v[3]=v4[this->GetFieldValue("Station[4]", hid)];
  v[4]='\0';
  std::string temp=v;
  return temp;
}

int sTgcHitIdHelper::GetPhiSector(const int& hid) const {
  return this->GetFieldValue("PhiSector", hid);
}

int sTgcHitIdHelper::GetZSector(const int& hid) const {
  return this->GetFieldValue("ZSector", hid);
}

//----sTgc
int sTgcHitIdHelper::GetMultiLayer(const int& hid) const {
  return this->GetFieldValue("MultiLayer", hid);
}
int sTgcHitIdHelper::GetLayer(const int& hid) const {
  return this->GetFieldValue("Layer", hid);
}
int sTgcHitIdHelper::GetSide(const int& hid) const {
  return this->GetFieldValue("Side", hid);
}


//packing method
int sTgcHitIdHelper::BuildsTgcHitId(const std::string statName, const int phiSect,
                                    const int zSect, const int multiLayer, const int layer, const int side) const {

  int theID(0);
  this->SetStationName(statName, theID);
  this->SetFieldValue("PhiSector", phiSect, theID);
  this->SetFieldValue("ZSector", zSect, theID);
  this->SetFieldValue("MultiLayer", multiLayer, theID);
  this->SetFieldValue("Layer", layer, theID);
  this->SetFieldValue("Side", side, theID);
  return theID;
}
