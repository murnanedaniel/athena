/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef _TRIGGER_TRIGT1_TRIGT1CALOCALIBTOOLS_L1CALOCELLS2TRIGGERTOWERS_H_
#define _TRIGGER_TRIGT1_TRIGT1CALOCALIBTOOLS_L1CALOCELLS2TRIGGERTOWERS_H_

// Athena Includes
#include "AsgTools/AsgTool.h"
#include "TrigT1CaloCalibToolInterfaces/IL1CaloMatchCell2Tower.h"
#include "GaudiKernel/ToolHandle.h"
#include "CxxUtils/checker_macros.h"
#include "StoreGate/StoreGateSvc.h"

// Calo includes
#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "CaloIdentifier/TileID.h"
#include "CaloTriggerTool/CaloTriggerTowerService.h"
#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCellContainer.h"

// LAr includes
#include "LArRawEvent/LArDigitContainer.h"
#include "LArRawEvent/LArDigit.h"
#include "LArCabling/LArCablingLegacyService.h"

// Tile includes
#include "TileEvent/TileDigitsContainer.h"
#include "TileEvent/TileDigits.h"
#include "TileConditions/TileCablingService.h"

// PABC interface include
#include "TrigT1CaloCalibToolInterfaces/IL1CaloCells2TriggerTowers.h"

// STL includes
#include <map>
#include <vector>
#include <string>


/**
  * AlgTool that implement several methods
  * used to produced the On/Off and On/Attr
  * trigger tower Id maps.
  *
  *  @author Damien Prieur <damien.prieur@cern.ch>
  *  @author John Morris <john.morris@cern.ch>
  */

namespace LVL1 {

class ATLAS_NOT_THREAD_SAFE   // use of LArCablingBase
L1CaloCells2TriggerTowers final : virtual public IL1CaloCells2TriggerTowers, public asg::AsgTool
{
  /// Create a proper constructor for Athena
  ASG_TOOL_CLASS( L1CaloCells2TriggerTowers , LVL1::IL1CaloCells2TriggerTowers )
public:
  /// constructor
  L1CaloCells2TriggerTowers( const std::string& name );
  /// destructor
  virtual ~L1CaloCells2TriggerTowers() {}
  /// delete the big 4
  L1CaloCells2TriggerTowers() = delete;
  L1CaloCells2TriggerTowers(const L1CaloCells2TriggerTowers& rhs) = delete;
  L1CaloCells2TriggerTowers(L1CaloCells2TriggerTowers&& rhs) = delete;
  L1CaloCells2TriggerTowers& operator=(const L1CaloCells2TriggerTowers& rhs) = delete;

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  // To associate CaloCell to the trigger tower they belong to
  // This method must be called before any other method in order
  // to setup correctly the CaloCell-TT map
  bool initCaloCellsTriggerTowers(const CaloCellContainer& cellContainer);
  bool initLArDigitsTriggerTowers(const LArDigitContainer& larDigitContainer);
  bool initTileDigitsTriggerTowers(const TileDigitsContainer& tileDigitsContainer);

  // John Morris <john.morris@cern.ch> 15/2/2011
  // Moving most useful functions to the top of this file for clarity
  virtual std::vector<const CaloCell*> caloCells(const Identifier& ttId) const override;
  virtual std::vector<std::vector<const CaloCell*> > caloCellsByLayer(const Identifier& ttId) const override;
  virtual std::vector<int> layerNames(const Identifier& ttID) const override;
  virtual int layerNames(const CaloCell* cell) const override;
  virtual float energy(const std::vector<const CaloCell*> &cells) const override;
  virtual float et(const std::vector<const CaloCell*> &cells) const override;
  virtual float energy(const Identifier& ttId) const override;
  virtual float et(const Identifier& ttId) const override;

  // Return a vector digits belonging the to requested trigger tower
  std::vector<double> samples(const Identifier& ttId) const;
  std::vector<const LArDigit*> larDigits(const Identifier& ttId) const;
  std::vector<const TileDigits*> tileDigits(const Identifier& ttId) const;

  void dumpCaloCells() const;
  void dumpCaloCells(const Identifier& ttId) const;

  void dumpDigits(const Identifier& ttId) const;

  void dumpLArDigits(const Identifier& ttId) const; /* ttId must be a LAr TTId */
  void dumpLArDigits() const;

  void dumpTileDigits(const Identifier& ttId) const; /* ttId must be a Tile TTId */
  void dumpTileDigits() const;

protected:
  float calcEnergyOrEt(const std::vector<const CaloCell*> &cells, const unsigned int mode) const;

  void reset();
  void resetCaloCells();
  void resetLArDigits();
  void resetTileDigits();


  void dump(const std::vector<const CaloCell*>& vCells) const;
  void dump(const std::vector<const LArDigit*>& vCells) const;
  void dump(const std::vector<const TileDigits*>& vCells) const;

private:
  // Id helpers
  const CaloLVL1_ID* m_lvl1Helper;
  const TileID* m_tileID;
  const CaloCell_ID* m_caloCellHelper;

  // Services
  LArCablingLegacyService* m_larCablingSvc;
  ToolHandle<LVL1::IL1CaloMatchCell2Tower> m_cellMatch
    { this, "L1CaloMatchCell2Tower",  "LVL1::L1CaloMatchCell2Tower",  "L1CaloMatchCell2Tower" };
  const TileCablingService * m_tileCablingService;
  CaloTriggerTowerService* m_ttSvc;

  bool m_bInitialized;
  bool m_bLArDigitsInitialized;
  bool m_bTileDigitsInitialized;

  // Offline TT Id - vector of CaloCell map
  std::map<unsigned int, std::vector<const CaloCell*> > m_mTTCaloCells;
  std::map<unsigned int, std::vector<const LArDigit*> > m_mTTLArDigits;
  std::map<unsigned int, std::vector<const TileDigits*> > m_mTTTileDigits;

};
} // end of namespace
#endif
