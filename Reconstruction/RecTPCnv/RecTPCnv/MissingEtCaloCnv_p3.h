/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef RECTPCNV_MISSINGETCALOCNV_P3_H
#define RECTPCNV_MISSINGETCALOCNV_P3_H

/********************************************************************
NAME:     MissingEtCaloCnv_p3.h
PURPOSE:  Transient/Persisten converter for MissingEtCalo class
********************************************************************/

#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "RecTPCnv/MissingEtCalo_p3.h"
#include "MissingETEvent/MissingEtCalo.h"

class MsgStream;

class MissingEtCaloCnv_p3: public T_AthenaPoolTPCnvConstBase<MissingEtCalo, MissingEtCalo_p3>
{
    public:
        using base_class::transToPers;
        using base_class::persToTrans;

            
	MissingEtCaloCnv_p3() {};
        
 	virtual void persToTrans( const MissingEtCalo_p3    *persObj,
                                  MissingEtCalo             *transObj,
                                  MsgStream                 &msg ) const override;

        virtual void transToPers( const MissingEtCalo       *transObj,
                                  MissingEtCalo_p3          *persObj,
                                  MsgStream                 &msg ) const override;
};


#endif 
