/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//****************************************************************************//
//*                                                                          *//
//*            Original code from Andrea Salamon                             *//
//*                                                                          *//
//****************************************************************************//
//
// 20/01/2004     A. Nisati      First code cleaning + use of BaseObject
// 17/07/2009     A. Salamon     Added output from all BCs + fixed one bug
//                               (mask for bcid output in SectorLogic::output)
// 21/02/2011     M. Corradi &   Fixed overlap resolution algorithm for new
//                A. Salamon     RoI geometry
//
//****************************************************************************//
//
#include "TrigT1RPChardware/SectorLogic.h"

using namespace std;

//
// Some usefull functions
//****************************************************************************//
// returns bits from lsb to msb right aligned
CMAword getbits(CMAword x, int msb, int lsb) { return ~(~0u << (msb - lsb + 1)) & (x >> lsb); }
//****************************************************************************//
// returns bits from lsb to msb right aligned
// sets bits between lsb and msb of x to first (msb - lsb + 1) right
// bits of y, remaining bits are left unchanged
unsigned setbits(unsigned x, int msb, int lsb, unsigned y) {
    unsigned mask = 0u;
    mask = (~0u << (msb + 1)) | ~(~0u << lsb);
    return (x & mask) | (y << lsb);
}
//****************************************************************************//
// returns bits from lsb to msb right aligned
// ### DataFromPad ###

// constructor of the class
// sets all the fields to 0

//****************************************************************************//
// prints the mask for DataFromPad
ostream &dfpa(ostream &stream, int indent, int whitesp, int ntimes) {
    int iind = 0, iwhit = 0, itime = 0;
    for (iind = 0; iind <= indent - 1; iind++) stream << " ";
    for (itime = 0; itime <= ntimes - 1; itime++) {
        for (iwhit = 0; iwhit <= whitesp - 1; iwhit++) stream << " ";
        stream << "b o r o o p m n s";
    }
    stream << std::endl;
    for (iind = 0; iind <= indent - 1; iind++) stream << " ";
    for (itime = 0; itime <= ntimes - 1; itime++) {
        for (iwhit = 0; iwhit <= whitesp - 1; iwhit++) stream << " ";
        stream << "c p s e f t t t g";
    }
    stream << std::endl;
    return stream;
}
//****************************************************************************//
// overload of the << operator for DataFromPad
ostream &operator<<(ostream &stream, DataFromPad o) {
    stream.width(3);
    stream.fill('x');
    stream << o.bcid << " ";
    stream.width(1);
    stream.fill(' ');
    stream << o.opl << " ";
    stream << o.r << " ";
    stream << o.oveta << " ";
    stream << o.ovphi << " ";
    stream << o.pt << " ";
    stream << o.roi << " ";
    stream << o.ntrig << " ";
    stream << o.sign;
    return stream;
}

//****************************************************************************//
// prints the mask for OutputFromSectorLogic
ostream &ofsla(ostream &stream, int indent, int whitesp, int ntimes) {
    int iind = 0, iwhit = 0, itime = 0;
    for (iind = 0; iind <= indent - 1; iind++) stream << " ";
    for (itime = 0; itime <= ntimes - 1; itime++) {
        for (iwhit = 0; iwhit <= whitesp - 1; iwhit++) stream << " ";
        stream << "ss  bnnppoorproorprn";
    }
    stream << std::endl;
    for (iind = 0; iind <= indent - 1; iind++) stream << " ";
    for (itime = 0; itime <= ntimes - 1; itime++) {
        for (iwhit = 0; iwhit <= whitesp - 1; iwhit++) stream << " ";
        stream << "gg  cttttefsdiefsdit";
    }
    stream << std::endl;
    for (iind = 0; iind <= indent - 1; iind++) stream << " ";
    for (itime = 0; itime <= ntimes - 1; itime++) {
        for (iwhit = 0; iwhit <= whitesp - 1; iwhit++) stream << " ";
        stream << "21   21212222211111 ";
    }
    stream << std::endl;
    return stream;
}
//****************************************************************************//
// overload of the << operator
ostream &operator<<(ostream &stream, OutputFromSectorLogic &o) {
    //    stream << "x x ";
    stream.setf(std::ios::dec, std::ios::basefield);
    //  stream << setw(3);
    stream << o.sign2;
    stream << o.sign1;
    stream.width(3);
    stream.fill('x');
    stream << o.bcid;
    stream.width(1);
    stream.fill(' ');
    stream << o.ntrig2;
    stream << o.ntrig1;
    stream << o.pt2;
    stream << o.pt1;
    stream << o.ove2;
    stream << o.ovf2;
    stream << o.r2;
    stream << o.pad2;
    stream << o.roi2;
    stream << o.ove1;
    stream << o.ovf1;
    stream << o.r1;
    stream << o.pad1;
    stream << o.roi1;
    stream << o.ntrig;
    return stream;
}
//****************************************************************************//

// ### InternalRegister ###

// overload of the << operator
ostream &operator<<(ostream &stream, InternalRegister &o) {
    int j = 0;
    for (j = 0; j <= 7; j++) {
        stream << "pad[" << j << "] : ";
        stream << o.pad[j] << std::endl;
    }
    stream.setf(ios::hex, ios::basefield);
    stream << "tile   : " << o.tile << std::endl;
    stream << "sl out : " << o.out;
    stream.setf(ios::hex, ios::basefield);
    return stream;
}
//****************************************************************************//
// ### SectorLogic ###

// standard constructor of the class
SectorLogic::SectorLogic(int run, int event, CMAword debug, ubit16 subsys, ubit16 sect, bool oldSimulation) :
    BaseObject(Hardware, "SectorLogic") {
    // identificazione sector logic di aleandro
    m_run = run;
    m_event = event;
    m_debug = debug;
    m_subsys = subsys;
    m_sect = sect;
    m_oldSimulation = oldSimulation;

    //
    // define m_nBunMax
    //
    m_nBunMax = NOBXS;
    // reset TC and OPL check parameters
    // EnableTCCheck and m_EnableOPLCheck
    m_EnableTCCheckLow = 0x00000000;
    m_EnableTCCheckHigh = 0x00000000;
    m_EnableOPLCheck = 0x00000000;
    // m_SetTCCheck and m_SetOPLCheck
    int j = 0;
    int k = 0;
    for (j = 0; j <= 7; j++) {
        // low pT
        for (k = 0; k <= 2; k++) {
            m_SetTCCheck[j][k] = 0;
            m_SetOPLCheck[j][k] = 0;
        }
        // high pT
        for (k = 3; k <= 5; k++) { m_SetTCCheck[j][k] = 0; }
    }
}
//****************************************************************************//
// destructor of the class
SectorLogic::~SectorLogic(void) {}
//****************************************************************************//

CMAword SectorLogic::outputToMuCTPI(int deltaBC) {
    int bunchID = BCZERO + deltaBC;
    if (bunchID < m_nBunMax && bunchID >= 0) {
        ubit16 bxsafe = (ubit16)bunchID;
        return output(bxsafe);
    } else {
        throw std::out_of_range("SectorLogic::outputToMuCTPI: bunchID out of range: " + std::to_string(bunchID));
    }
}

// OLD VERSION
/*
CMAword SectorLogic::outputToMuCTPI(ubit16 bunchID) {
  ubit16 bxsafe=BCZERO;
  if( bunchID <= m_nBunMax-1 ) {
    bxsafe=bunchID;
  } else {
    cout << "warning : bunchID out of range, set to default value" << bxsafe << std::endl;
  }
  return output(bxsafe);
}
*/

// initializes the array
void SectorLogic::init(void) {}
//****************************************************************************//
// @@@@@ check the internal registers of the Sector Logic
void SectorLogic::check(void) {
    int bcid;

    for (bcid = 0; bcid <= m_nBunMax - 1; bcid++) {
        cout << "LowPtFilter_in BCID is " << m_LowPtFilter_in[bcid].out.bcid << std::endl
             << "LowPtFilter_out BCID is " << m_LowPtFilter_out[bcid].out.bcid << std::endl

             << "TileCalConfirm_in BCID is " << m_TileCalConfirm_in[bcid].out.bcid << std::endl
             << "TileCalConfirm_out BCID is " << m_TileCalConfirm_out[bcid].out.bcid << std::endl

             << "SolveEtaOverlap_in BCID is " << m_SolveEtaOverlap_in[bcid].out.bcid << std::endl
             << "SolveEtaOverlap_out BCID is " << m_SolveEtaOverlap_out[bcid].out.bcid << std::endl

             << "SortHighest_in BCID is " << m_SortHighest_in[bcid].out.bcid << std::endl
             << "SortHighest_out BCID is " << m_SortHighest_out[bcid].out.bcid << std::endl

             << "Sort2ndHighest_in BCID is " << m_Sort2ndHighest_in[bcid].out.bcid << std::endl
             << "Sort2ndHighest_out BCID is " << m_Sort2ndHighest_out[bcid].out.bcid << std::endl;
    }
}
//****************************************************************************//
// tile cal check configuration
void SectorLogic::LoadTCCheck(CMAword EnableTCCheckLow_in, CMAword EnableTCCheckHigh_in, CMAword SetTCCheck_in[8][6]) {
    m_EnableTCCheckLow = EnableTCCheckLow_in;
    m_EnableTCCheckHigh = EnableTCCheckHigh_in;
    int i = 0;
    int j = 0;
    for (i = 0; i <= 7; i++) {
        for (j = 0; j <= 5; j++) { m_SetTCCheck[i][j] = SetTCCheck_in[i][j]; }
    }
}
//****************************************************************************//
// opl check configuration
void SectorLogic::LoadOPLCheck(CMAword EnableOPLCheck_in, ubit16 SetOPLCheck_in[8][3]) {
    m_EnableOPLCheck = EnableOPLCheck_in;
    int i = 0;
    int j = 0;
    for (i = 0; i <= 7; i++) {
        for (j = 0; j <= 2; j++) { m_SetOPLCheck[i][j] = SetOPLCheck_in[i][j]; }
    }
}
//****************************************************************************//
// @@@@@ input, output and clock methods @@@@@
void SectorLogic::dbginput(ubit16 bx, DataFromPad from_pad[8], CMAword from_tile_cal) {
    //  int BCID_now = bx%8;

    int ipad = 0;
    for (ipad = 0; ipad <= 7; ipad++) { m_InFromPad[bx][ipad] = from_pad[ipad]; }
    m_InFromTileCal[bx] = from_tile_cal;
}
//****************************************************************************//
void SectorLogic::load(ubit16 padAdd, ubit16 BX, ubit16 RoIAdd, ubit16 pT, ubit16 OPL, ubit16 overlapPhi, ubit16 overlapEta,
                       ubit16 RoiAmbiguity, ubit16 BCIDcounter) {
    m_InFromPad[BX][padAdd].bcid = BCIDcounter;
    m_InFromPad[BX][padAdd].r = RoiAmbiguity;
    m_InFromPad[BX][padAdd].oveta = overlapEta;
    m_InFromPad[BX][padAdd].ovphi = overlapPhi;
    m_InFromPad[BX][padAdd].opl = OPL;
    m_InFromPad[BX][padAdd].pt = pT;
    m_InFromPad[BX][padAdd].roi = RoIAdd;

    m_InFromTileCal[BX] = 0xff;

    /*
    cout << "input from pad : BC = " << BX << " padAdd = " << padAdd
         << " pT = "<<  pT << " roi = " << RoIAdd << " bcid = " << BCIDcounter << std::endl;
    */
}
//****************************************************************************//
OutputFromSectorLogic SectorLogic::dbgoutput(ubit16 i) { return m_OutFromSectorLogic[i]; }
//****************************************************************************//
CMAword SectorLogic::output(ubit16 i) {
    CMAword fmtout = 0;

    ubit16 slroi1 = 0;
    ubit16 slroi2 = 0;

    if (m_oldSimulation) {
        if (m_OutFromSectorLogic[i].pt1)  // aleandro addendum 6-10-2003
            slroi1 = (m_OutFromSectorLogic[i].pad1) * 4 + (m_OutFromSectorLogic[i].roi1) + 1;
        if (m_OutFromSectorLogic[i].pt2)  // aleandro addendum 6-10-2003
            slroi2 = (m_OutFromSectorLogic[i].pad2) * 4 + (m_OutFromSectorLogic[i].roi2) + 1;
    } else {
        if (m_OutFromSectorLogic[i].pt1)  // aleandro addendum 6-10-2003
            slroi1 = (m_OutFromSectorLogic[i].pad1) * 4 + (m_OutFromSectorLogic[i].roi1);
        if (m_OutFromSectorLogic[i].pt2)  // aleandro addendum 6-10-2003
            slroi2 = (m_OutFromSectorLogic[i].pad2) * 4 + (m_OutFromSectorLogic[i].roi2);
    }

    // MC 2015/7/7 add bc information
    ubit16 bc = 0;
    if (i > BCZERO) {
        bc = i - BCZERO;
    } else if (i < BCZERO) {
        bc = i + 8 - BCZERO;
    }

    if (m_OutFromSectorLogic[i].pt1 == 0) m_OutFromSectorLogic[i].pt1 = 7;
    if (m_OutFromSectorLogic[i].pt2 == 0) m_OutFromSectorLogic[i].pt2 = 7;

    fmtout = fmtout | (((m_OutFromSectorLogic[i].ntrig > 2) & 0x01) << 0);  // >2 candidates in a sector
    fmtout = fmtout | ((slroi1 & 0x1f) << 1);                               // ROI1
    fmtout = fmtout | ((0 & 0x03) << 6);                                    // Reserved1
    fmtout = fmtout | ((m_OutFromSectorLogic[i].ovf1 & 0x01) << 8);         // Phi overlap1
    fmtout = fmtout | ((m_OutFromSectorLogic[i].ove1 & 0x01) << 9);         // Eta overlap1
    fmtout = fmtout | ((slroi2 & 0x1f) << 10);                              // ROI2
    fmtout = fmtout | ((0 & 0x03) << 15);                                   // Reserved2
    fmtout = fmtout | ((m_OutFromSectorLogic[i].ovf2 & 0x01) << 17);        // Phi overlap2
    fmtout = fmtout | ((m_OutFromSectorLogic[i].ove2 & 0x01) << 18);        // Eta overlap2
    fmtout = fmtout | ((m_OutFromSectorLogic[i].pt1 & 0x07) << 19);         // Pt1
    fmtout = fmtout | ((m_OutFromSectorLogic[i].pt2 & 0x07) << 22);         // Pt2
    fmtout = fmtout | ((m_OutFromSectorLogic[i].r1 & 0x01) << 25);          // >1 candidate in ROI1
    fmtout = fmtout | ((m_OutFromSectorLogic[i].r2 & 0x01) << 26);          // >1 candidate in ROI2
    fmtout = fmtout | ((bc & 0x07) << 27);                                  // BCID
    fmtout = fmtout | ((m_OutFromSectorLogic[i].sign1 & 0x01) << 30);       // Candidate1 sign
    fmtout = fmtout | ((m_OutFromSectorLogic[i].sign2 & 0x01) << 31);       // Candidate2 sign

    return fmtout;
}
//****************************************************************************//
void SectorLogic::execute() {
    // executes sector logic algorithm for all input bunches
    int ibx = 0;
    for (ibx = 0; ibx <= m_nBunMax - 1; ibx++) {
        // the content of the SL input registers are copied in the 1st input
        // register
        int ipad = 0;
        for (ipad = 0; ipad <= 7; ipad++) { m_LowPtFilter_in[ibx].pad[ipad] = m_InFromPad[ibx][ipad]; }
        m_LowPtFilter_in[ibx].tile = m_InFromTileCal[ibx];

        // @@@@@ 1st step registers
        // @@@@@ low Pt filter
        // if there is a low Pt track in the pad and the OPL check option is on
        // looks for OPL confirmation in the other pads mapped for OPL check
        m_LowPtFilter_in[ibx].out.bcid = m_LowPtFilter_in[ibx].pad[0].bcid;
        m_LowPtFilter_out[ibx] = m_LowPtFilter_in[ibx];
        int i1 = 0;
        // external cycle on the pads
        int OPLfrompads = 0;
        for (i1 = 0; i1 <= 7; i1++) { OPLfrompads = setbits(OPLfrompads, i1, i1, m_LowPtFilter_in[ibx].pad[i1].opl); }
        // external cycle on the pads
        for (i1 = 0; i1 <= 7; i1++) {
            // check for low Pt track in the pad
            int Pt_reg1 = 0;
            Pt_reg1 = m_LowPtFilter_in[ibx].pad[i1].pt;
            if (0 < Pt_reg1 && Pt_reg1 <= 3) {
                // check for OPL check option for the given pad (depends on track Pt!!)
                if (getbits(m_EnableOPLCheck, i1 * 4 + Pt_reg1, i1 * 4 + Pt_reg1) == 1) {
                    int OPLCheck = 0;
                    // check for OPL confirmation in all mapped pads
                    OPLCheck = OPLfrompads & m_SetOPLCheck[i1][Pt_reg1 - 1];
                    // if the Pt track were not confirmed the track Pt is reset to 0
                    if (OPLCheck == 0) m_LowPtFilter_out[ibx].pad[i1].pt = 0;
                }
            }
        }

        // @@@@@ 2nd step registers
        // @@@@@ tile cal confirmation
        // if there is track in the pad and the tile cal confirmation option is on
        // looks for the tile cal confirmation in the mapped tile cal signals
        m_TileCalConfirm_in[ibx] = m_LowPtFilter_out[ibx];
        m_TileCalConfirm_out[ibx] = m_TileCalConfirm_in[ibx];
        int j1 = 0;
        // external cycle on the pads
        for (j1 = 0; j1 <= 7; j1++) {
            // check for low Pt track in the pad
            int Pt_reg2 = 0;
            Pt_reg2 = m_TileCalConfirm_in[ibx].pad[j1].pt;
            if (0 < Pt_reg2 && Pt_reg2 <= 6) {
                // check for TileCal confirmation option for the given pad
                // (depends on track Pt !!)
                // chose EnableTCCheck depending on pT
                sbit32 EnableTCCheck = 0;
                if (0 < Pt_reg2 && Pt_reg2 <= 3) {
                    EnableTCCheck = m_EnableTCCheckLow;
                } else {
                    EnableTCCheck = m_EnableTCCheckHigh;
                }
                if (getbits(EnableTCCheck, j1 * 4 + Pt_reg2, j1 * 4 + Pt_reg2) == 1) {
                    int TileCalCheck = 0;
                    // check for TileCal confirmation in all the pads mapped
                    TileCalCheck = m_TileCalConfirm_in[ibx].tile & m_SetTCCheck[j1][Pt_reg2 - 1];
                    // if the Pt track has not been confirmed the track Pt is reset to 0
                    if (TileCalCheck == 0) { m_TileCalConfirm_out[ibx].pad[j1].pt = 0; }
                }
            }
        }

        // @@@@@ 3rd step registers
        // @@@@@ solve eta overlaps within sector
        /* ALGORITMO DI SELEZIONE DELLE SOGLIE:
           TRA DUE TRACCE IN OVERLAP PASSA LA TRACCIA CON SOGLIA PIU' ALTA
           A PARITA' DI SOGLIA PASSA QUELLA CON ETA MINORE */
        m_SolveEtaOverlap_in[ibx] = m_TileCalConfirm_out[ibx];
        m_SolveEtaOverlap_out[ibx] = m_SolveEtaOverlap_in[ibx];
        // first check : if the overlap flag is on there must be a valid track
        // in the pad, otherwise send a warning
        int k1 = 0;
        for (k1 = 0; k1 <= 7; k1++) {
            if (m_SolveEtaOverlap_out[ibx].pad[k1].oveta && m_SolveEtaOverlap_out[ibx].pad[k1].pt == 0) {
                if (m_debug) {
                    std::cout << "pad # " << k1 << " bcid # " << ibx << " has eta overlap flag on but no triggered track" << std::endl;
                }
            }
        }
        // run the overlap resolution algorithm on EVEN and then on ODD pads
        int kk = 0;
        for (kk = 0; kk <= 1; kk++) {
            int k3 = 0;
            for (k3 = kk; k3 <= 6 - kk; k3 = k3 + 2) {
                if (m_SolveEtaOverlap_out[ibx].pad[k3].oveta || m_SolveEtaOverlap_out[ibx].pad[k3 + 1].oveta) {
                    if (m_SolveEtaOverlap_out[ibx].pad[k3].oveta && m_SolveEtaOverlap_out[ibx].pad[k3 + 1].oveta) {
                        if ((m_SolveEtaOverlap_out[ibx].pad[k3].roi == 2 && m_SolveEtaOverlap_out[ibx].pad[k3 + 1].roi == 3) ||
                            (m_SolveEtaOverlap_out[ibx].pad[k3].roi == 0 && m_SolveEtaOverlap_out[ibx].pad[k3 + 1].roi == 1) ||
                            (m_SolveEtaOverlap_out[ibx].pad[k3].roi == 3 && m_SolveEtaOverlap_out[ibx].pad[k3 + 1].roi == 2) ||
                            (m_SolveEtaOverlap_out[ibx].pad[k3].roi == 1 && m_SolveEtaOverlap_out[ibx].pad[k3 + 1].roi == 0)) {
                            // set to 0 the lower Pt
                            if (m_SolveEtaOverlap_out[ibx].pad[k3].pt >= m_SolveEtaOverlap_out[ibx].pad[k3 + 1].pt) {
                                m_SolveEtaOverlap_out[ibx].pad[k3 + 1].pt = 0;
                                m_SolveEtaOverlap_out[ibx].pad[k3 + 1].oveta = 0;
                            } else {
                                m_SolveEtaOverlap_out[ibx].pad[k3].pt = 0;
                                m_SolveEtaOverlap_out[ibx].pad[k3].oveta = 0;
                            }
                        } else {
                            if (m_debug) {
                                std::cout << "pads " << k3 << " and " << k3 + 1 << " have eta overlap flags on with wrong RoIs"
                                          << std::endl;
                            }
                        }
                    }
                }
            }
        }

        // @@@@@ 4th step registers
        // @@@@@ sort highest track
        m_SortHighest_in[ibx] = m_SolveEtaOverlap_out[ibx];
        m_SortHighest_out[ibx] = m_SortHighest_in[ibx];
        int Pt_reg4 = 0;
        int ROI_reg4 = 0;
        int pad_reg4 = 0;
        int l1 = 0;
        // external cycle on the pads
        for (l1 = 0; l1 <= 7; l1++) {
            // check for a track with higher Pt
            if (m_SortHighest_in[ibx].pad[l1].pt > Pt_reg4) {
                Pt_reg4 = m_SortHighest_in[ibx].pad[l1].pt;
                ROI_reg4 = m_SortHighest_in[ibx].pad[l1].roi;
                pad_reg4 = l1;
            }
        }
        // if there is a validated track
        if (Pt_reg4 > 0) {
            // put the result in the output part of the register
            m_SortHighest_out[ibx].out.pt1 = Pt_reg4;
            m_SortHighest_out[ibx].out.pad1 = pad_reg4;
            m_SortHighest_out[ibx].out.roi1 = ROI_reg4;
            m_SortHighest_out[ibx].out.ovf1 = m_SortHighest_in[ibx].pad[pad_reg4].ovphi;
            m_SortHighest_out[ibx].out.ove1 = m_SortHighest_in[ibx].pad[pad_reg4].oveta;
            m_SortHighest_out[ibx].out.ntrig1 = m_SortHighest_in[ibx].pad[pad_reg4].ntrig;
            m_SortHighest_out[ibx].out.sign1 = m_SortHighest_in[ibx].pad[pad_reg4].sign;
            m_SortHighest_out[ibx].out.r1 = m_SortHighest_in[ibx].pad[pad_reg4].r;
            // set to 1 the number of valid tracks
            m_SortHighest_out[ibx].out.ntrig = 1;
            // and put the ouput internal register to 0
            // otherwise the track will be counted twice !!!!!
            m_SortHighest_out[ibx].pad[pad_reg4].pt = 0;
            m_SortHighest_out[ibx].pad[pad_reg4].roi = 0;
        }

        // @@@@@ 5th step registers
        // @@@@@ sort 2nd highest track
        m_Sort2ndHighest_in[ibx] = m_SortHighest_out[ibx];
        m_Sort2ndHighest_out[ibx] = m_Sort2ndHighest_in[ibx];
        int Pt_reg5 = 0;
        int ROI_reg5 = 0;
        int pad_reg5 = 0;
        int m1 = 0;
        // external cycle on the pads
        for (m1 = 0; m1 <= 7; m1++) {
            // check for a track with 2nd higer Pt
            if (m_Sort2ndHighest_in[ibx].pad[m1].pt > Pt_reg5 &&
                (m_Sort2ndHighest_in[ibx].pad[m1].pt != Pt_reg4 || m_Sort2ndHighest_in[ibx].pad[m1].roi != ROI_reg4 || m1 != pad_reg4)) {
                Pt_reg5 = m_Sort2ndHighest_in[ibx].pad[m1].pt;
                ROI_reg5 = m_Sort2ndHighest_in[ibx].pad[m1].roi;
                pad_reg5 = m1;
            }
        }

        // if there is a validated track
        if (Pt_reg5 > 0) {
            // put the result in the output part of the register
            m_Sort2ndHighest_out[ibx].out.pt2 = Pt_reg5;
            m_Sort2ndHighest_out[ibx].out.pad2 = pad_reg5;
            m_Sort2ndHighest_out[ibx].out.roi2 = ROI_reg5;
            m_Sort2ndHighest_out[ibx].out.ovf2 = m_Sort2ndHighest_in[ibx].pad[pad_reg5].ovphi;
            m_Sort2ndHighest_out[ibx].out.ove2 = m_Sort2ndHighest_in[ibx].pad[pad_reg5].oveta;
            m_Sort2ndHighest_out[ibx].out.ntrig2 = m_Sort2ndHighest_in[ibx].pad[pad_reg5].ntrig;
            m_Sort2ndHighest_out[ibx].out.sign2 = m_Sort2ndHighest_in[ibx].pad[pad_reg5].sign;
            m_Sort2ndHighest_out[ibx].out.r2 = m_Sort2ndHighest_in[ibx].pad[pad_reg5].r;
            // set to 2 the number of valid tracks
            m_Sort2ndHighest_out[ibx].out.ntrig = 2;
            // and put the ouput internal register to 0
            // othrwise the track will be counted twice !!!!!
            m_Sort2ndHighest_out[ibx].pad[pad_reg5].pt = 0;
            m_Sort2ndHighest_out[ibx].pad[pad_reg5].roi = 0;
        }

        // check for other candidates in sector
        for (m1 = 0; m1 <= 7; m1++) {
            if (m_Sort2ndHighest_out[ibx].pad[m1].pt > 0) (m_Sort2ndHighest_out[ibx].out.ntrig)++;
        }

        // the content of the 5th output internal register is copied in the SL output register
        m_OutFromSectorLogic[ibx] = m_Sort2ndHighest_out[ibx].out;
    }
}
//****************************************************************************//
// overload of the << operator
std::ostream &operator<<(std::ostream &stream, SectorLogic &o) {
    stream << "@@@@@@@@@@ event and sector logic identification @@@@@@@@@@\n\n";

    stream << "run    = " << o.m_run << std::endl;
    stream << "event  = " << o.m_event << std::endl;
    stream << "debug  = " << o.m_debug << std::endl;
    stream << "subsys = " << o.m_subsys << std::endl;
    stream << "sect   = " << o.m_sect << std::endl;

    stream << std::endl;

    // print all the parameters of the sector logic board
    stream << "@@@@@@@@@@ sector logic configuration parameters @@@@@@@@@@\n\n";

    // tccheck
    stream.setf(ios::hex, ios::basefield);
    stream << "EnableTCCheckLow   : ";
    stream.width(8);
    stream.fill('0');
    stream << o.m_EnableTCCheckLow << std::endl;
    stream << "EnableTCCheckHigh  : ";
    stream.width(8);
    stream.fill('0');
    stream << o.m_EnableTCCheckHigh << std::endl;
    int jj = 0;
    int kk = 0;
    for (jj = 0; jj <= 7; jj++) {
        stream << "SetTCCheck pad[";
        stream << jj;
        stream << "]  : ";
        for (kk = 5; kk >= 0; kk--) {
            stream.width(8);
            stream.fill('0');
            stream << o.m_SetTCCheck[jj][kk] << "  ";
        }
        stream << std::endl;
    }
    stream << std::endl;

    // opl check
    stream << "EnableOPLCheck     : ";
    stream.width(8);
    stream.fill('0');
    stream << o.m_EnableOPLCheck << std::endl;
    for (jj = 0; jj <= 7; jj++) {
        stream << "SetOPLCheck pad[";
        stream << jj;
        stream << "] : ";
        for (kk = 2; kk >= 0; kk--) {
            stream.width(2);
            stream.fill('0');
            stream << o.m_SetOPLCheck[jj][kk] << "  ";
        }
        stream << std::endl;
    }
    stream << std::endl;

    stream.width(1);
    stream.fill(' ');

    // internal register (input register)
    InternalRegister *intreginp[5];
    intreginp[0] = o.m_LowPtFilter_in.data();
    intreginp[1] = o.m_TileCalConfirm_in.data();
    intreginp[2] = o.m_SolveEtaOverlap_in.data();
    intreginp[3] = o.m_SortHighest_in.data();
    intreginp[4] = o.m_Sort2ndHighest_in.data();
    // internal register (output register)
    InternalRegister *intregoutp[5];
    intregoutp[0] = o.m_LowPtFilter_out.data();
    intregoutp[1] = o.m_TileCalConfirm_out.data();
    intregoutp[2] = o.m_SolveEtaOverlap_out.data();
    intregoutp[3] = o.m_SortHighest_out.data();
    intregoutp[4] = o.m_Sort2ndHighest_out.data();

    int ibx = 0;

    // print the input registers of the sector logic board
    stream << "@@@@@@@@@@ sector logic input registers @@@@@@@@@@\n\n";

    // print input from pads
    stream.setf(std::ios::dec, std::ios::basefield);
    dfpa(stream, 8, 8, NOBXS);
    int ipad = 0;
    for (ipad = 0; ipad <= 7; ipad++) {
        stream << "pad[" << ipad << "] :";
        for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
            stream << "      ";
            stream << o.m_InFromPad[ibx][ipad];
        }
        stream << std::endl;
    }

    // print input from tilecal
    stream.setf(std::ios::hex, std::ios::basefield);
    stream << "tile   :";
    for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
        stream << "                 ";
        stream.width(8);
        stream.fill('0');
        stream << o.m_InFromTileCal[ibx];
    }
    stream.width(1);
    stream.fill(' ');
    stream << std::endl << std::endl;

    // print all  sector logic internal registers
    stream << "@@@@@@@@@@ sector logic internal registers @@@@@@@@@@\n\n";

    int ireg = 0;
    for (ireg = 0; ireg <= 4; ireg++) {
        // input registers
        stream << "internal registers # " << ireg + 1 << " (input)" << std::endl;
        // DataFromPad
        stream.setf(std::ios::dec, std::ios::basefield);
        dfpa(stream, 8, 8, NOBXS);
        for (ipad = 0; ipad <= 7; ipad++) {
            stream << "pad[" << ipad << "] :";
            for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
                stream << "      ";
                stream << (intreginp[ireg] + ibx)->pad[ipad];
            }
            stream << std::endl;
        }
        // Tile Cal
        stream.setf(ios::hex, ios::basefield);
        stream << "tile   :";
        for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
            stream << "                 ";
            stream.width(8);
            stream.fill('0');
            stream << (intreginp[ireg] + ibx)->tile;
        }
        stream.width(1);
        stream.fill(' ');
        stream << std::endl;
        // Sector Logic Output
        stream.setf(std::ios::dec, std::ios::basefield);
        stream << "sl out :";
        for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
            stream << "     ";
            stream << (intreginp[ireg] + ibx)->out;
        }
        stream << std::endl;
        stream << std::endl;

        // output registers
        stream << "internal registers # " << ireg + 1 << " (output)" << std::endl;
        // DataFromPad
        stream.setf(std::ios::dec, std::ios::basefield);
        dfpa(stream, 8, 8, NOBXS);
        for (ipad = 0; ipad <= 7; ipad++) {
            stream << "pad[" << ipad << "] :";
            for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
                stream << "      ";
                stream << (intregoutp[ireg] + ibx)->pad[ipad];
            }
            stream << std::endl;
        }
        // Tile Cal
        stream.setf(std::ios::hex, std::ios::basefield);
        stream << "tile   :";
        for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
            stream << "                 ";
            stream.width(8);
            stream.fill('0');
            stream << (intregoutp[ireg] + ibx)->tile;
        }
        stream.width(1);
        stream.fill(' ');
        stream << std::endl;
        // Sector Logic Output
        stream.setf(std::ios::dec, std::ios::basefield);
        stream << "sl out :";
        for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
            stream << "     ";
            stream << (intregoutp[ireg] + ibx)->out;
        }
        stream << std::endl << std::endl;
    }

    // print the output registers of the sector logic board
    stream.setf(std::ios::dec, std::ios::basefield);
    stream << "@@@@@@@@@@ sector logic output register @@@@@@@@@@\n\n";
    ofsla(stream, 8, 5, NOBXS);
    stream << "        ";
    for (ibx = 0; ibx <= NOBXS - 1; ibx++) {
        stream << "     ";
        stream << o.m_OutFromSectorLogic[ibx];
    }

    return stream;
}
//****************************************************************************//
