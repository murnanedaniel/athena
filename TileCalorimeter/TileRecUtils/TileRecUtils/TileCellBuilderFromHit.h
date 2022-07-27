/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TILERECUTILS_TILECELLBUILDERFROMHIT_H
#define TILERECUTILS_TILECELLBUILDERFROMHIT_H

/********************************************************************
 *
 * NAME:     TileCellBuilderFromHit
 * PACKAGE:  offline/TileCalorimeter/TileRecUtils
 *
 * AUTHOR :  A. Solodkov
 * CREATED:  10-Oct-2015
 *
 * PURPOSE:  Create Cells from Hits and store them in container
 *
 *  Input: TileHit from TileHitContainer
 *  Output: Container or collection with TileCells
 *  Parameters:
 *    TileHitContainer - Name of input container
 *   
 ********************************************************************/

// Tile includes
#include "TileEvent/TileCellContainer.h"
#include "TileEvent/TileHitContainer.h"
#include "TileIdentifier/TileFragHash.h"
#include "TileIdentifier/TileRawChannelUnit.h"
#include "TileRecUtils/TileCellBuilder.h"
#include "TileConditions/ITileBadChanTool.h"
#include "TileConditions/TileCondToolEmscale.h"
#include "TileConditions/TileCondToolTiming.h"
#include "TileConditions/TileCablingSvc.h"
#include "TileConditions/TileSamplingFraction.h"

// Calo includes
#include "CaloInterface/ICaloCellMakerTool.h"
#include "CaloConditions/CaloAffectedRegionInfo.h"
#include "CaloConditions/CaloNoise.h"

// Atlas includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "Identifier/HWIdentifier.h"
#include "AthenaKernel/IOVSvcDefs.h"
#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteHandleKey.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaKernel/IAthRNGSvc.h"

// Gaudi includes
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"
#include "CLHEP/Random/RandomEngine.h"

// C++ STL includes
#include <string>
#include <vector>
#include <memory>

// forward declarations
class TileID;
class TileTBID;
class TileHWID;
class TileCell;
class MbtsDetDescrManager;
class TileDetDescrManager;
class TileCellCollection;
class CaloCellContainer;
class TileCablingService;


/**
 @class TileCellBuilderFromHit
 @brief This class creates Cells from RawChannels and stores them in a container
 
 */
class TileCellBuilderFromHit
  : public extends<AthAlgTool, ICaloCellMakerTool>
{
  public:
    TileCellBuilderFromHit(const std::string& type, const std::string& name, const IInterface* parent); //!< Contructor

    virtual ~TileCellBuilderFromHit(); //!< Destructor

    virtual StatusCode initialize() override;

    virtual StatusCode finalize() override;

    /// method to process all raw channels and store them in container
    virtual StatusCode process (CaloCellContainer* theCellContainer,
                                const EventContext& ctx) const override;

    //AlgTool InterfaceID
    static const InterfaceID& interfaceID();
    //static const InterfaceID& interfaceID() { return ICaloCellMakerTool; };

  private:
    /// status of every drawer
    typedef TileDrawerEvtStatus TileDrawerEvtStatusArray[5][64];

    // properties
    SG::ReadCondHandleKey<CaloNoise> m_caloNoiseKey{this, "CaloNoise",
                                                    "electronicNoise",
                                                    "CaloNoise object to read"};
    SG::ReadHandleKey<TileHitContainer> m_hitContainerKey{this, "TileHitContainer", 
                                                          "TileHitCnt", 
                                                          "Input Tile hit container key"};

    SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey{this, "EventInfo", 
                                                      "EventInfo", 
                                                      "EventInfo key"};


    SG::WriteHandleKey<TileCellContainer> m_MBTSContainerKey{this, "MBTSContainer", 
                                                             "MBTSContainer", 
                                                             "Output Tile MBTS container key"};

    SG::WriteHandleKey<TileCellContainer> m_E4prContainerKey{this, "E4prContainer", 
                                                             "E4prContainer",
                                                             "Output Tile E4 prime container key"};
    /**
     * @brief Name of TileSamplingFraction in condition store
     */
    SG::ReadCondHandleKey<TileSamplingFraction> m_samplingFractionKey{this,
        "TileSamplingFraction", "TileSamplingFraction", "Input Tile sampling fraction"};


    std::string m_infoName;

    float m_eneForTimeCut;        //!< keep time for channels with energy above cut
    float m_eneForTimeCutMBTS;    //!< similar cut for MBTS in pC
    float m_zeroEnergy;           //!< energy to store in every PMT if both PMT are bad
    int m_qualityCut;           //!< cut on channel quality (set energy to m_zeroEnergy for them)

    float m_maxTime;              //!< maximum time for the PMTs in the cels
    float m_minTime;              //!< minimum time for the PMTs in the cels
    bool m_maskBadChannels;      //!< if true=> bad channels are masked
    float m_noiseSigma;          //!< cell electronic noise if CaloNoise is switched off 

    const TileID* m_tileID{nullptr};   //!< Pointer to TileID
    const TileTBID* m_tileTBID{nullptr}; //!< Pointer to TileTBID
    const TileHWID* m_tileHWID{nullptr}; //!< Pointer to TileHWID
    const TileCablingService* m_cabling{nullptr}; //!< Pointer to TileCabling

    ServiceHandle<IAthRNGSvc> m_rndmSvc  //!< Random number service to use
      { this, "RndmSvc", "AthRNGSvc", "Random Number Service used in TileCellBuildetFromHit" };

    ToolHandle<ITileBadChanTool> m_tileBadChanTool{this,
        "TileBadChanTool", "TileBadChanTool", "Tile bad channel tool"};

    ToolHandle<TileCondToolEmscale> m_tileToolEmscale{this,
        "TileCondToolEmscale", "TileCondToolEmscale", "Tile EM scale calibration tool"};

    /**
     * @brief Name of Tile cabling service
     */
    ServiceHandle<TileCablingSvc> m_cablingSvc{ this,
        "TileCablingSvc", "TileCablingSvc", "The Tile cabling service"};


    const TileDetDescrManager* m_tileMgr{nullptr}; //!< Pointer to TileDetDescrManager
    const MbtsDetDescrManager* m_mbtsMgr{nullptr}; //!< Pointer to MbtsDetDescrManager

    TileFragHash::TYPE m_RChType;        //!< Type of TileRawChannels (Fit, OF2, etc.)
    //unsigned int m_bsflags;              //!< other flags stored in TileRawChannelContainer

    // These were accumulated, but never actually used.
    // They also spoil reentrancy, so leave them commented-out for now.
    // If this information is needed in the future, these can be changed
    // to use atomics.
    ///TileDrawerRunStatus m_drawerRunStatus[5][64]; //!< overall status of drawer in whole run
    //int m_eventErrorCounter[4]; //!< number of events with no errors(0), warnings(1), error(2), total(3)

    std::vector<CaloAffectedRegionInfo> m_affectedRegionInfo_global;
    std::vector<CaloAffectedRegionInfo> m_affectedRegionInfo_current_run;

    //!< method to process raw channels from a given vector and store them in collection
    template<class ITERATOR, class COLLECTION>
    void build(const CaloNoise* caloNoise,
               TileDrawerEvtStatusArray& drawerEvtStatus,
               const ITERATOR & begin,
               const ITERATOR & end,
               COLLECTION * coll,
               TileCellContainer* MBTSCells,
               TileCellContainer* E4prCells,
               const TileSamplingFraction* samplingFraction) const;
               

    /** method to check if channels are good or bad. Puts zero if both channels are bad
     or recovers from single-channel failure. It returns true if cell was changed, false otherwise
     */
    bool maskBadChannel (TileDrawerEvtStatusArray& drawerEvtStatus,
                         TileCell* pCell) const;
    bool maskBadChannels (TileDrawerEvtStatusArray& drawerEvtStatus,
                          TileCell* pCell, bool single_PMT_C10, bool Ecell) const;

    void correctCell(TileCell* pCell, int correction, int pmt, int gain, float ener, float time,
        unsigned char iqual, unsigned char qbit) const; //!< Compute calibrated energy, time, etc. for TileCell and adjust it.

    unsigned char iquality(float qual) const  {//!< method to compute the cell quality
         return std::min(255, abs((int) qual));
    } // keep quality within 8 bits make it "unsigned char"

    unsigned char qbits(TileDrawerEvtStatusArray& drawerEvtStatus,
                        int ros, int drawer, bool count_over, bool good_time, bool good_ener,
        bool overflow, bool underflow, bool good_overflowfit) const; //!< method to compute the cell quality bits

    bool m_RUN2;
    bool m_RUN2plus;
    int m_E1_TOWER;

    static const int NSIDE = 2;
    static const int NPHI = 8;
    static const int NETA = 2;
    static const int NCELLMBTS = NSIDE * NPHI * NETA;
    inline int mbts_index(int side, int phi, int eta) const {
      return (side * NPHI + phi) * NETA + eta;
    }
    static const int E4SIDE = -1;
    static const int E4ETA  = 2;
    static const int E4NPHI = 4;
    static const int NCELLE4PR = E4NPHI;
    inline int e4pr_index(int phi) const {
      return  phi;
    }

};

#endif
