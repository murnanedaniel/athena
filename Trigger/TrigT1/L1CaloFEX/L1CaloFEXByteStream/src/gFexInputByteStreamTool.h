/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//***************************************************************************
//                           gFexInputByteStreamTool  -  This tool decodes Run3 gFEX input data!
//                              -------------------
//     begin                : 10 08 2022
//     email                : cecilia.tosciri@cern.ch
//  ***************************************************************************/

#ifndef GFEXINPUTBYTESTREAMTOOL_H 
#define GFEXINPUTBYTESTREAMTOOL_H

// Trigger includes
#include "TrigT1ResultByteStream/IL1TriggerByteStreamTool.h"
#include "xAODTrigL1Calo/gFexTowerContainer.h"
#include "xAODTrigL1Calo/gFexTowerAuxContainer.h"

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"

// Gaudi includes
#include "Gaudi/Property.h"

#include "gFexPos.h"
namespace gPos = LVL1::gFEXPos;

/** @class gFexInputByteStreamTool
 *  @brief Implementation of a tool for L1 input data conversion from BS to xAOD and from xAOD to BS
 *  (IL1TriggerByteStreamTool interface)
 **/

typedef  std::array<std::array<uint32_t, 7>,  100>        gfiber;
typedef  std::array<std::array<int,      6>,   32>        gEngines;
typedef  std::array<std::array<int,      12>,  32>        gtFPGA;
typedef  std::array<std::array<int,      20>,  100>       gFields;
typedef  std::array<std::array<int,      16>,  100>       gCaloTwr;
typedef  std::array<std::array<char,     20>,  100>       gFieldsChar;

typedef  std::array<std::array<int,      20>,  4>         gType;
typedef  std::array<std::array<char,     20>,  4>         gTypeChar;


class gFexInputByteStreamTool : public extends<AthAlgTool, IL1TriggerByteStreamTool> {
    public:
        gFexInputByteStreamTool(const std::string& type, const std::string& name, const IInterface* parent);
        virtual ~gFexInputByteStreamTool() override = default;

        // ------------------------- IAlgTool methods --------------------------------
        virtual StatusCode initialize() override;

        // ------------------------- IL1TriggerByteStreamTool methods ----------------------
        /// BS->xAOD conversion
        virtual StatusCode convertFromBS(const std::vector<const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment*>& vrobf, const EventContext& eventContext)const override;

        /// xAOD->BS conversion
        virtual StatusCode convertToBS(std::vector<OFFLINE_FRAGMENTS_NAMESPACE_WRITE::ROBFragment*>& vrobf, const EventContext& eventContext) override;

        /// Declare ROB IDs for conversion
        virtual const std::vector<uint32_t>& robIds() const override {
            return m_robIds.value();
        }


    private:
        // ------------------------- Properties --------------------------------------
        // ROBIDs property required by the interface
        Gaudi::Property<std::vector<uint32_t>> m_robIds {this, "ROBIDs", {}, "List of ROB IDs required for conversion to/from xAOD RoI"};

         //Write handle keys for the L1Calo EDMs for BS->xAOD mode of operation
        SG::WriteHandleKey< xAOD::gFexTowerContainer> m_gTowersWriteKey   {this,"gTowersWriteKey"  ,"L1_gTowers","Write gFexEDM Trigger Tower container"};
        
        // Read handle keys for the L1Calo EDMs for xAOD->BS mode of operation
        SG::ReadHandleKey < xAOD::gFexTowerContainer> m_gTowersReadKey    {this,"gTowersReadKey"   ,"L1_gTowers","Read gFexEDM Trigger Tower container"};

        virtual void a_gtrx_map( const gfiber &inputData, gfiber &jf_lar_rx_data) const;

        virtual void b_gtrx_map( const gfiber &inputData, gfiber &jf_lar_rx_data) const;

        virtual void c_gtrx_map( const gfiber &inputData, gfiber &outputData) const;
        
        virtual void gtReconstructABC(  int XFPGA, 
                                        gfiber Xfiber,  int Xin, 
                                        gtFPGA &Xgt, int *BCIDptr,
                                        int do_lconv, 
                                        std::array<int, gPos::MAX_FIBERS> XMPD_NFI,
                                        std::array<int, gPos::MAX_FIBERS>  XCALO_TYPE,
                                        gCaloTwr XMPD_GTRN_ARR,
                                        gType XMPD_DSTRT_ARR,  
                                        gTypeChar XMPD_DTYP_ARR,
                                        std::array<int, gPos::MAX_FIBERS> XMSK) const;

        virtual int crc9d32(std::array<int, 6> inWords,int numWords,int reverse) const;

        virtual int crc9d23(int inword, int in_crc, int  reverse ) const;

        virtual void undoMLE(int &datumPtr ) const;

        virtual void gtRescale(gtFPGA twr, gtFPGA &twrScaled, int scale) const;
        
        virtual void getEtaPhi(float &Eta, float &Phi, int iEta, int iPhi) const;
        
};

#endif // GFEXINPUTBYTESTREAMTOOL_H
