/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TileHitCollection_H       
#define TileHitCollection_H       
                                         
#include "TileEvent/TileHit.h" 
#include "TileEvent/TileRawDataCollection.h"
                                         
class TileHitCollection : public TileRawDataCollection<TileHit>
{
  friend class TileHitContainerCnv;
public:

  typedef TileRawDataCollection<TileHit> MyBase;
  typedef MyBase::ID ID;
  typedef TileHit _TElement;

  TileHitCollection ( ID id, SG::OwnershipPolicy ownPolicy=SG::OWN_ELEMENTS )
    : MyBase(id,ownPolicy) { }
  TileHitCollection ( SG::OwnershipPolicy ownPolicy=SG::OWN_ELEMENTS )
    : MyBase(ownPolicy) { }


  /**
   * @brief Copy constructor.
   * @param other Collection to copy.
   * Does a deep copy of the contents.
   */
  TileHitCollection (const TileHitCollection& other);

  
  ~TileHitCollection() { }
};

CLASS_DEF(TileHitCollection,2924,0)
           
#endif

