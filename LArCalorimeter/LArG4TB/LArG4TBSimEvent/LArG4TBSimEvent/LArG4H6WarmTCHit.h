/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//=================================================
// LArG4H6WarmTCHit class
//
// information stored
//=================================================

#ifndef LArG4H6WarmTCHit_h
#define LArG4H6WarmTCHit_h 1

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
class LArG4H6WarmTCHit 
{
 public:
     LArG4H6WarmTCHit(){;};
     LArG4H6WarmTCHit(int ad, double en){addr = ad; energy = en; };
     LArG4H6WarmTCHit(const LArG4H6WarmTCHit &right);
//     ~LArG4H6WarmTCHit(){;};

     LArG4H6WarmTCHit& operator=(const LArG4H6WarmTCHit &right);
     int operator==(const LArG4H6WarmTCHit &right);
//     inline void *operator new(size_t, void* aHit = 0);
//     inline void operator delete(void *aHit);

     void Draw(){;};
     void Print(){;};

//Set- Get- methods
 // energy deposit
     inline void SetEnergy(double ed) { energy = ed; }
     inline void AddEnergy(double ed) { energy += ed; }
     inline double GetEnergy() { return energy; }

 // address
     inline void SetAddr(int d) { addr = d; }
     inline int GetAddr() { return addr; }

 private:
     int     addr;   // Address: 
     double  energy;
};

typedef std::vector<LArG4H6WarmTCHit> LArG4H6WarmTCHitsCollection;

// extern G4Allocator<LArG4H6WarmTCHit> LArG4H6WarmTCHitAllocator;

/*
inline void* LArG4H6WarmTCHit::operator new(size_t, void* aHit)
{
// void *aHit;
 aHit = (void*) LArG4H6WarmTCHitAllocator.MallocSingle();
 return aHit;
}

inline void LArG4H6WarmTCHit::operator delete(void *aHit)
{
 LArG4H6WarmTCHitAllocator.FreeSingle( (LArG4H6WarmTCHit*) aHit);
}
*/
#endif
