/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_CALOTILERETRIEVER_H
#define JIVEXML_CALOTILERETRIEVER_H

#include <string>
#include <vector>
#include <cstddef>
#include <map>

#include "TileConditions/TileCondToolTiming.h"
#include "TileConditions/TileCondToolEmscale.h"
#include "TileConditions/ITileBadChanTool.h"

#include "CaloIdentifier/CaloCell_ID.h"

#include "JiveXML/IDataRetriever.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

class IToolSvc;

class Identifier;
class CaloCellContainer;

namespace JiveXML{
  
  /**
   * @class CaloTileRetriever
   * @brief Retrieves all @c Tile Calo Cell @c objects 
   *
   *  - @b Properties
   *    - StoreGateKeyTile: default is 'AllCalo'. Don't change.
   *	- CallThreshold: default is 50 MeV
   *	- RetrieveTile: activate retriever, default is true
   *	- DoTileDigits: write Tile digits (ADC), default is false
   *	- DoTileCellDigits: more verbose on cell details
   *	- CellEnergyPrec: output precision, default is 3 digits
   *	- CellTimePrec: output precision, default is 3 digits
   *   
   *  - @b Retrieved @b Data
   *    - location in phi and eta
   *    - identifier and adc counts of each cell 
   *    - various pmt details
   */
  class CaloTileRetriever : virtual public IDataRetriever,
                                   public AthAlgTool {
    
    public:
      
      /// Standard Constructor
      CaloTileRetriever(const std::string& type,const std::string& name,const IInterface* parent);
      
      /// Retrieve all the data
      virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool); 
      const DataMap getCaloTileData(const CaloCellContainer* cellContainer);

      /// Return the name of the data type
      virtual std::string dataTypeName() const { return m_typeName; };
	
      ///Default AthAlgTool methods
      StatusCode initialize();

    private:
      ///The data type that is generated by this retriever
      const std::string m_typeName;

      ToolHandle<TileCondToolTiming> m_tileToolTiming{this,
          "TileCondToolTiming", "TileCondToolTiming", "Tile timing tool"};

      ToolHandle<TileCondToolEmscale> m_tileToolEmscale{this,
          "TileCondToolEmscale", "TileCondToolEmscale", "Tile EM scale calibration tool"};

      ToolHandle<ITileBadChanTool> m_tileBadChanTool{this,
          "TileBadChanTool", "TileBadChanTool", "Tile bad channel tool"};

      void calcTILELayerSub(Identifier&);
      const CaloCell_ID*   m_calocell_id;
    
      std::string m_sgKey; 
      double m_cellThreshold;
      int m_cellEnergyPrec;
      int m_cellTimePrec;
      bool m_tile;
      bool m_doTileDigit;
      bool m_doTileCellDetails;
      bool m_doBadTile;

      DataVect m_sub;
  };
}
#endif
