/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/*
 * Author: Noel Dawe <Noel-dot-Dawe-AT-cern-dot-ch>
 */

#ifndef TAUDETAILS_H
#define TAUDETAILS_H

#include <string>

// Add new details to the proper block below (4 blocks: int/float tau details and int/float event details)
// Use all caps, but note that in the end, detail names are case insensitive
// since they are always compared in the uppercase state when building the discriminants at init
// Using this "ENUM_OR_STRING" macro allows as to use the same list to init enums (below)
// as well as arrays of strings (in the same order) (see TauDetailsManager.cxx)

#define INT_TAU_DETAILS \
    ENUM_OR_STRING( AUTHOR ), \
    ENUM_OR_STRING( TAU_PI0_N ), \
    ENUM_OR_STRING( CHARGE ), \
    ENUM_OR_STRING( NUMTRACK ), \
    ENUM_OR_STRING( NUMWIDETRACK ), \
    ENUM_OR_STRING( NSTRIP ), \
    ENUM_OR_STRING( NISOLLOOSETRK ), \
    ENUM_OR_STRING( NPI0 ), \
    ENUM_OR_STRING( NPI0CL), \
    ENUM_OR_STRING( NISOLTRK ), \
    ENUM_OR_STRING( NUMTOPOCLUSTERS ), \
    ENUM_OR_STRING( __IntTauDetail__END__ )

#define FLOAT_TAU_DETAILS \
    ENUM_OR_STRING( ETA ), \
    ENUM_OR_STRING( ABS_ETA ), \
    ENUM_OR_STRING( ABS_ETA_LEAD_TRACK ), \
    ENUM_OR_STRING( TAU_ABSDELTAETA ), \
    ENUM_OR_STRING( TAU_ABSDELTAPHI ), \
    ENUM_OR_STRING( PHI ), \
    ENUM_OR_STRING( E ), \
    ENUM_OR_STRING( ET ), \
    ENUM_OR_STRING( PT ), \
    ENUM_OR_STRING( EMFRACTIONATEMSCALE ), \
    ENUM_OR_STRING( EMRADIUS ), \
    ENUM_OR_STRING( HADRADIUS ), \
    ENUM_OR_STRING( CALRADIUS ), \
    ENUM_OR_STRING( ISOLFRAC ), \
    ENUM_OR_STRING( CENTFRAC ), \
    ENUM_OR_STRING( STRIPWIDTH2 ), \
    ENUM_OR_STRING( TRFLIGHTPATHSIG ), \
    ENUM_OR_STRING( IPSIGLEADTRK ), \
    ENUM_OR_STRING( IPSIGLEADLOOSETRK ), \
    ENUM_OR_STRING( ETOVERPTLEADTRK ), \
    ENUM_OR_STRING( PTLEADTRKOVERET ), \
    ENUM_OR_STRING( IPZ0SINTHETASIGLEADTRK ), \
    ENUM_OR_STRING( MASSTRKSYS ), \
    ENUM_OR_STRING( TRKWIDTH2 ), \
    ENUM_OR_STRING( TRKAVGDIST ), \
    ENUM_OR_STRING( ETEFLOWOVERET ), \
    ENUM_OR_STRING( MEFLOW ), \
    ENUM_OR_STRING( SUMPT3TRK ), \
    ENUM_OR_STRING( SUMPT ), \
    ENUM_OR_STRING( DRMIN ), \
    ENUM_OR_STRING( DRMAX ), \
    ENUM_OR_STRING( SUMPT_OVER_ET ), \
    ENUM_OR_STRING( SUMPT3TRK_OVER_ET ), \
    ENUM_OR_STRING( ETHAD_EM_OVER_SUMPT3TRK ), \
    ENUM_OR_STRING( ETEM_EM_OVER_SUMPT3TRK ), \
    ENUM_OR_STRING( ETHAD_EM_OVER_SUMPT ), \
    ENUM_OR_STRING( ETEM_EM_OVER_SUMPT ), \
    ENUM_OR_STRING( TRT_NHT_OVER_NLT ), \
    ENUM_OR_STRING( M ), \
    ENUM_OR_STRING( TOPOINVMASS ), \
    ENUM_OR_STRING( EFFTOPOINVMASS ), \
    ENUM_OR_STRING( JVF ), \
    ENUM_OR_STRING( PT_PILEUP ), \
    ENUM_OR_STRING( TRACK_ISO ), \
    ENUM_OR_STRING( CALO_ISO ), \
    ENUM_OR_STRING( CALO_ISO_CORRECTED ), \
    ENUM_OR_STRING( LEAD2CLUSTEREOVERALLCLUSTERE ), \
    ENUM_OR_STRING( LEAD3CLUSTEREOVERALLCLUSTERE ), \
    ENUM_OR_STRING( HADLEAKET ),\
    ENUM_OR_STRING( SUMEMCELLETOVERLEADTRKPT ),\
    ENUM_OR_STRING( SECMAXSTRIPET ),\
    ENUM_OR_STRING( BDTJETSCORE ),\
    ENUM_OR_STRING( CHPIEMEOVERCALOEME ),\
    ENUM_OR_STRING( PSSFRACTION ),\
    ENUM_OR_STRING( EMPOVERTRKSYSP ),\
    ENUM_OR_STRING( EMEOVERTRKSYSE ),		\
    ENUM_OR_STRING( CORRCENTFRAC ),		\
    ENUM_OR_STRING( CORRFTRK ),		\
    ENUM_OR_STRING( TAU_PI0_VISTAU_M ), \
    ENUM_OR_STRING( TAU_PTRATIO ), \
    ENUM_OR_STRING( INTERAXIS_ETA ), \
    ENUM_OR_STRING( INTERAXIS_PHI ), \
    ENUM_OR_STRING( PI0CL1_PT ), \
    ENUM_OR_STRING( PI0CL1_ETA ), \
    ENUM_OR_STRING( PI0CL1_PHI ), \
    ENUM_OR_STRING( PI0CL2_PT ), \
    ENUM_OR_STRING( PI0CL2_ETA ), \
    ENUM_OR_STRING( PI0CL2_PHI ), \
    ENUM_OR_STRING( VISTAU_PI0CL_PT ),	\
    ENUM_OR_STRING( VISTAU_PI0CL_ETA ), \
    ENUM_OR_STRING( VISTAU_PI0CL_PHI ), \
    ENUM_OR_STRING( VISTAU_PI0CL_M ), \
    ENUM_OR_STRING( EMFRACTIONATEMSCALE_MOVEE3 ), \
    ENUM_OR_STRING( TAU_SEEDTRK_SECMAXSTRIPETOVERPT ), \
    ENUM_OR_STRING( __FloatTauDetail__END__ )
    
#define INT_EVENT_DETAILS \
    ENUM_OR_STRING( NUMVERTICES ), \
    ENUM_OR_STRING( NUMGOODVERTICES ),\
    ENUM_OR_STRING( NUM_PILEUP_AND_PRIMARY_VERTICES ),\
    ENUM_OR_STRING( __IntEventDetail__END__ )

#define FLOAT_EVENT_DETAILS \
    ENUM_OR_STRING( __FloatEventDetail__END__ )

#undef ENUM_OR_STRING
#define ENUM_OR_STRING( x ) x

namespace Details
{
    enum IntTauDetail
    {
        INT_TAU_DETAILS
    };
    
    enum FloatTauDetail
    {
        FLOAT_TAU_DETAILS
    };
    
    enum IntEventDetail
    {
        INT_EVENT_DETAILS
    };
    
    enum FloatEventDetail
    {
        FLOAT_EVENT_DETAILS
    };
}

#endif
