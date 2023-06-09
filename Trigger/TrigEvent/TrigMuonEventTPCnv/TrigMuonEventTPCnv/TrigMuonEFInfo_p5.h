/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**********************************************************************************
 * @Project: Trigger
 * @Package: TrigMuonEventTPCnv
 * @Class  : TrigMuonEFInfo_p5
 *
 * @brief persistent partner for TrigMuonEFInfo
 *
 * @author Andrew Hamilton  <Andrew.Hamilton@cern.ch>  - U. Geneva
 * @author Francesca Bucci  <F.Bucci@cern.ch>          - U. Geneva
 * @author Sergio Grancagnolo  <Sergio.Grancagnolo@le.infn.it> - U.Salento/INFN Le
 * @author Alexander Oh  <alexander.oh@cern.ch>        - U. Manchester
 *
 **********************************************************************************/
#ifndef TRIGMUONEVENTTPCNV_TRIGMUONEFINFO_P5_H
#define TRIGMUONEVENTTPCNV_TRIGMUONEFINFO_P5_H

#include <stdint.h>
#include <string>
#include "AthenaPoolUtilities/TPObjRef.h"

class TrigMuonEFInfo_p5
{
	friend class TrigMuonEFInfoCnv;

public:

	TrigMuonEFInfo_p5() : m_etaPreviousLevel(0.0), m_phiPreviousLevel(0.0) {}
	virtual ~TrigMuonEFInfo_p5(){}

        //private:
    // unsigned short int m_roi;
    // unsigned short int m_nSegments;
    // unsigned short int m_nMdtHits;
    // unsigned short int m_nRpcHits;
    // unsigned short int m_nTgcHits;
    // unsigned short int m_nCscHits;
	// this array holds all the unsigned ints from above. in that order.
    unsigned short int m_allTheInts[6];
    
	float m_etaPreviousLevel;
	float m_phiPreviousLevel;
    // 
    // TPObjRef m_spectrometerTrack; // probably not needed
    // TPObjRef m_extrapolatedTrack; // probably not needed
    // TPObjRef m_combinedTrack;     // probably not needed
    // 
	TPObjRef m_trackContainer;

};

#endif
