/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// ====================================================================
/*
        TGCTMDB.cc
                                      
*/
// ====================================================================

#include <iostream>
#include <iomanip>

#include "TrigT1TGC/TGCTMDB.h"
#include "TrigT1TGC/TGCTMDBOut.h"

namespace LVL1TGCTrigger {

using namespace std;

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////
TGCTMDB::TGCTMDB()
//////////////////////
{
  for (int side=0; side < 2; side++) {
    for (int mod=0; mod < NumberOfTileModules(); mod++) {
      buffer[side*NumberOfTileModules() + mod] = new TGCTMDBOut(side, mod); 
    }
  }
}

//////////////////////
TGCTMDB::~TGCTMDB()
//////////////////////
{
  for (int idx=0; idx<2*NumberOfTileModules(); idx++){
    delete buffer[idx]; 
  }
}

///////////////////////////////////////////////////////////////
TGCTMDB::TGCTMDB(const TGCTMDB& right)
/////////////////////////////////////////////////////////////
{
  *this= right;
}


/////////////////////////////////////////////////////////////
TGCTMDB& TGCTMDB::operator=(const TGCTMDB& right)
/////////////////////////////////////////////////////////////
{
  if (this != &right) {
    for (int idx=0; idx<2*NumberOfTileModules(); idx++){
      (*buffer[idx]) = *(right.buffer[idx]);
    }
  }
  return *this;
}
  
/////////////////////////////////////////////////////////////
const TGCTMDBOut* TGCTMDB::getOutput(int side, int mod) const
/////////////////////////////////////////////////////////////
{
  if ( (side<0)||(side>1) ) return 0;
  if ( (mod<0)||(mod>=NumberOfTileModules()) ) return 0;
  return buffer[side*NumberOfTileModules() + mod] ;
}

/////////////////////////////////////////////////////////////
const TGCTMDBOut* TGCTMDB::getOutput(int side, int sector, int mod) const
/////////////////////////////////////////////////////////////
{
  if ((side<0)||(side>1)) return 0;
  if ((sector<0)||(sector>47)) return 0;
  if ((mod<0)||(mod>3)) return 0;
  int octant = sector / 6;
  int sec = sector % 6;
  int offset = 0;
  if      (sec==0) offset = -4;
  else if (sec==1) offset = -3;
  else if (sec==2) offset = -1;
  else if (sec==3) offset =  0;
  else if (sec==4) offset =  2;
  else if (sec==5) offset =  3;
  int moduleID = (octant*(NumberOfTileModules()/8) + offset + NumberOfTileModules()) % NumberOfTileModules();;
  moduleID = (moduleID + mod) % NumberOfTileModules();
  return buffer[side*NumberOfTileModules() + moduleID];
}

/////////////////////////////////////////////////////////////
void  TGCTMDB::setOutput(int side, int module, int hit56, int hit6)
/////////////////////////////////////////////////////////////
{
  if ( (side<0)||(side>1) ) return;
  if ( (module<0)||(module>=NumberOfTileModules()) ) return;
  buffer[side*NumberOfTileModules() +module]->SetHit56(hit56);
  buffer[side*NumberOfTileModules() +module]->SetHit6(hit6);
}

/////////////////////////////////////////////////////////////
void  TGCTMDB::eraseOutput()
/////////////////////////////////////////////////////////////
{
  for (int idx=0; idx<2*NumberOfTileModules(); idx++){
    buffer[idx]->Clear(); 
  }
}

/////////////////////////////
void TGCTMDB::Print() const
/////////////////////////////
{
  for (int idx=0; idx<2*NumberOfTileModules(); idx++){
    buffer[idx]->Print(); 
  }
}
  

} //end of namespace bracket
