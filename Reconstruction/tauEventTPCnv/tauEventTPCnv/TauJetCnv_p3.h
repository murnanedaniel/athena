/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/// @file TauJetCnv_p3.h
/// 
/// Implementation of transient<->persistent converter 
/// for Analysis::TauCommonDetails v3
/// 
/// Package: Reconstruction/tauEventTPCnv
/// 
/// @author Thomas Burgess  
/// @author Lukasz Janyst
/// 
/// Created by Lukasz Janyst in May 2008
///
/// $Id: TauJetCnv_p3.h,v 1.2 2009-01-27 17:38:43 tburgess Exp $


#ifndef tauEventTPCnv_TAUJETCNV_P3_H
#define tauEventTPCnv_TAUJETCNV_P3_H

#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "tauEventTPCnv/TauJet_p3.h"

class MsgStream;

/// Transient<->persistent converter for Analysis::TauCommonDetails v3
class TauJetCnv_p3: 
    public T_AthenaPoolTPCnvBase<Analysis::TauJet, TauJet_p3>
{
public:
    ///Constructor
    TauJetCnv_p3() {}

    ///Destructor
    virtual ~TauJetCnv_p3() {}

    ///Fills a transient object from a persistent object
    ///
    ///@param pers Persistent object
    ///@param trans Transient object
    ///@param msg Message stream
    virtual void persToTrans( 
	const TauJet_p3 *pers,
	Analysis::TauJet *trans,
	MsgStream &msg );

    ///Fills a persistent object from a transient object   
    ///
    ///@param trans Transient object
    ///@param pers Persistent object
    ///@param msg Message stream
    virtual void transToPers( 
	const Analysis::TauJet *trans,
	TauJet_p3 *pers,
	MsgStream &msg );
private:
};

#endif // tauEventTPCnv_TAUJETCNV_P3_H
