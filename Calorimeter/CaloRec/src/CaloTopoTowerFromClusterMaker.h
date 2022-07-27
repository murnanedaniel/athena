// -*- c++ -*- 
/* Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration */
#ifndef CALOREC_CALOTOPOTOWERFROMCLUSTERMAKER_H
#define CALOREC_CALOTOPOTOWERFROMCLUSTERMAKER_H

#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteHandleKey.h"

#include "GaudiKernel/ServiceHandle.h"

#include "AthenaBaseComps/AthAlgTool.h" 

#include "CaloRec/CaloTowerCollectionProcessor.h"
#include "CaloProtoCluster.h"

#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloClusterCellLink.h"
#include "CaloEvent/CaloCellClusterWeights.h"
#include "CaloEvent/CaloCellContainer.h"

#include "CaloGeoHelpers/CaloSampling.h"

#include "CaloDetDescr/CaloTowerGeometry.h"

#include "xAODCaloEvent/CaloClusterContainer.h"

#include <string>
#include <vector>
#include <bitset>
#include <map>

#ifndef _CALOTOPOTOWERFROMCLUSTERMAKER_BITSET_SIZE
#define _CALOTOPOTOWERFROMCLUSTERMAKER_BITSET_SIZE 28
#endif

class CaloTopoTowerFromClusterMaker : public AthAlgTool, virtual public CaloTowerCollectionProcessor
{
public:

  ///@brief Tool constructor
  CaloTopoTowerFromClusterMaker(const std::string& type,const std::string& name,const IInterface* pParent);
  ///@name @c AthAlgTool and @c CaloClusterCellProcessor interface implementations
  ///@{
  virtual StatusCode initialize() override;                   ///< Setting up the operational mode and corresponding parameters
  virtual StatusCode execute(const EventContext& ctx,
                             xAOD::CaloClusterContainer* pClusCont,
			                 CaloCellClusterWeights*     cellWeights) const override; ///< Execute the tool and fill the @c xAOD::CaloClusterContainer pointed to by @c pClusCont
  virtual StatusCode finalize() override;                        ///< Finalize the tool (no action)
  ///@}
  

private:

  ///@name Internally used types
  ///@{
  typedef std::vector<CaloProtoCluster> protocont_t; ///< Container for @c CaloProtoCluster objects
  typedef std::size_t                   uint_t;      ///< Unsigned integral type
  ///@}

  /// @name Tool properties
  /// @{
  /** @brief the name of the key of the CaloTowerGeometry object in the 
      ConditonsStore */
  SG::ReadCondHandleKey<CaloTowerGeometry>      m_towerGeoKey{this,"TowerGeometry","CaloTowerGeometry"};
  SG::ReadCondHandleKey<CaloDetDescrManager>    m_caloMgrKey{this,"CaloDetDescrManager", "CaloDetDescrManager"};

  SG::ReadHandleKey<xAOD::CaloClusterContainer> m_clusterContainerKey;                      ///< Topo-cluster container key
  SG::ReadHandleKey<CaloCellContainer>          m_cellContainerKey;                         ///< Calorimeter cell container
  bool                                          m_orderByPt = { false  };                   ///< Orders cluster container by @f$ p_{\text{T}} @f$, default @c true
  bool                                          m_prepareLCW = { false  };                  ///< Prepare LCW calibration, default is @c false
  bool                                          m_useCellsFromClusters = { true  };         ///< Use cells from topo-clusters if @c true, else use all cells, default is @c true
  bool                                          m_applyCellEnergyThreshold = { false  };    ///< Apply cell energy threshold, default is @c false 
  bool                                          m_doCellIndexCheck = { false  };            ///< Check cell hash index consistency if @c true (default @c false)
  bool                                          m_buildCombinedSignal = { false  };         ///< Build topo-clusters within given @f$ y @f$ range, else topo-towers
  double                                        m_energyThreshold;                          ///< Cell energy threshold, default is set in @c m_energyThresholdDef
  double                                        m_clusterRange;                             ///< Range where topo-clusters are used when <tt>m_buildCombinedSignal = true</tt>
  bool                                          m_removeSamplingData = { true };            ///< Remove sampling data for towers
  /// @}

  /// @name Constants and parameters
  /// @{
  //uint_t             m_numberOfCells;          ///< Number of cells (highest cell index + 1)
  //uint_t             m_maxCellHash;            ///< Maximum hash index of cell ( number of cells - 1)
  uint_t             m_numberOfSamplings;      ///< Number of samplings
  //uint_t             m_numberOfTowers;         ///< Number of towers
  static const double      m_energyThresholdDef;     ///< Default energy threshold
  static const double      m_clusterRangeDef;        ///< Default cluster @f$ y @f$ range
  static const uint_t      m_errorValueUINT;         ///< Error value for @c uint_t type values
  /// @}

  ///@name Internally used helpers
  ///@{
  static xAOD::CaloCluster::ClusterSize getClusterSize(uint_t etaBins,uint_t phiBins) ; ///< Returns a cluster size tag from number of eta and phi bins in tower grid
  static xAOD::CaloCluster::ClusterSize getClusterSize(uint_t towerBins) ;              ///< Returns a cluster size tag from number of towers (bins) in tower grid
  int cleanupCells(const CaloTowerGeometry* towerGeo, CaloClusterCellLink* clk,uint_t nclus) const;                      ///< Checks @c CaloClusterCellLink for consistency
  ///@}

  ///@name Tower builders
  ///
  ///@return @c false in case of problems with data access or inconsistent data structures 
  ///
  ///@param pCellCont reference to non-modifiable @c CaloCellContainer
  ///@param pProtoCont reference to @c CaloProtoCluster container filled on output.
  ///@param clusCont reference to non-modifiable @c xAOD::CaloClusterContainer
  ///@param protoCont reference to modifiable proto-cluster container
  ///
  ///@return 
  ///@{
  uint_t buildInclTowers(const CaloTowerGeometry* towerGeo, const CaloCellContainer& pCellCont,protocont_t& pProtoCont)const;            ///< Inclusive towers
  uint_t buildExclTowers(const CaloTowerGeometry* towerGeo, const CaloCellContainer& pCellCont,protocont_t& pProtoCont)const;            ///< Exclusive towers
  uint_t buildEMTopoTowers(const CaloTowerGeometry* towerGeo, const xAOD::CaloClusterContainer& clusCont,protocont_t& protoCont)const;   ///< EM topo-towers
  uint_t buildLCWTopoTowers(const CaloTowerGeometry* towerGeo, const xAOD::CaloClusterContainer& clusCont,protocont_t& protoCont,CaloCellClusterWeights* cellWeights) const;  ///< LCW topo-towers
  ///@}
  /// @brief Adding cells to proto-clusters
  ///
  /// @return @c true if cell successfully added to one (or more) proto-clusters
  ///
  /// @param cptr       pointer ton non-modifiable @c CaloCell object
  /// @param pProtoCont reference to proto-cluster container
  /// @param weight     additional (global) weight of cell (e.g. for geometrical weight for combined EM-scale signals)  
  bool addCellToProtoCluster(const CaloTowerGeometry* towerGeo, const CaloCell* cptr,protocont_t& pProtoCont,double weight=1.) const;
  
  ///@name Helpers
  ///@{
  static bool   filterProtoCluster(const CaloClusterCellLink& clnk)  ; ///< Checks for and removes invalid cell links  
  bool   checkCellIndices(const CaloTowerGeometry* towerGeo, const CaloDetDescrManager* caloDDM,
                          const CaloCellContainer* pCellCont) const; ///< Checks consistency between cell indices and hash identifiers
  bool   isValidIndex(uint_t idx)                             const; ///< Checks if argument is a valid index value 
  uint_t badIndexValue()                                      const; ///< Returns value indicating a bad index
  ///@}

  ///@name Excluded samplings
  ///@{
  std::vector<CaloSampling::CaloSample>                       m_excludedSamplings;         ///< List of excluded samplings (@c CaloSampling::CaloSample enumerators)
  std::vector<std::string>                                    m_excludedSamplingsName;     ///< List of excluded samplings (human-readable names)
  std::bitset< _CALOTOPOTOWERFROMCLUSTERMAKER_BITSET_SIZE >   m_excludedSamplingsPattern;  ///< Bit pattern indicates if sampling is excluded
  ///@}

  ///@name Monitoring
  ///@{
  ///@}
};

inline CaloTopoTowerFromClusterMaker::uint_t CaloTopoTowerFromClusterMaker::badIndexValue()          const { return m_errorValueUINT;       } 
inline bool                                  CaloTopoTowerFromClusterMaker::isValidIndex(uint_t idx) const { return idx != badIndexValue(); }

///@class CaloTopoTowerFromClusterMaker
///
/// @brief A cluster builder tool forming topo-clusters representing calorimeter tower signals on a regular grid in @f$ (\eta,\phi) @f$ space. By default,
///        EM-scale <i>topo-towers</i> are created from cells in topo-clusters.
///
/// This tool fills EM-scale towers and stores them as @c xAOD::CaloCluster. It supports several operational modes, which are
/// controlled by tool properties. It fills a container of type @c xAOD::CaloClusterContainer. The properties controlling its
/// specific behavior are:  
///
/// <table width="90%" align="center" style="border-width:0px;">
/// <tr><td align="center" colspan="4" style="background-color:rgb(245,245,220);color:rgb(165,42,42);"><b>Properties defining tool behavior</b></td></tr>
/// <tr style="color:rgb(131,42,34);">
/// <td align="left" width="25%">Property name</td>
/// <td align="center" width="25%">Property type</td>
/// <td align="center" width="25%">Default value</td>
/// <td align="left" width="25%">Comment</td>
/// </tr> 
/// <tr>
/// <td align="left" valign="top"><tt>OrderClusterByPt</tt></td>
/// <td align="center" valign="top"><tt>bool</tt></td>
/// <td align="center" valign="top"><tt>false</tt></td>
/// <td align="left" valign="top">if @c true, the @c xAOD::CaloClusterContainer is ordered by @f$ p_{\rm T}^{\rm clus} @f$.  See further comments below.</td>
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>PrepareLCW</tt></td>
/// <td align="center" valign="top"><tt>bool</tt></td>
/// <td align="center" valign="top"><tt>false</tt></td>
/// <td align="left" valign="top">if @c true, the tool fills a @c CaloCellClusterWeights object and records it into the event store to be used by @c CaloTopoClusterFromTowerCalibrator</td>
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>UseCellsFromClusters</tt></td>
/// <td align="center" valign="top"><tt>bool</tt></td>
/// <td align="center" valign="top"><tt>true</tt></td>
/// <td align="left" valign="top">if @c true, only cells from topo-clusters are used to fill the towers (<i>topo-towers</i>); else, @a inclusive @a towers are filled with all cells.</td>
/// </tr>
/// <tr><td align="center" colspan="4" style="background-color:rgb(245,245,220);color:rgb(165,42,42);"><b>Properties setting variables for operational modes</b></td></tr>
/// <tr style="color:rgb(131,42,34)">
/// <td align="left" width="25%">Property name</td>
/// <td align="center" width="25%">Property type</td>
/// <td align="center" width="25%">Default value</td>
/// <td align="left" width="25%">Comment</td>
/// </tr> 
/// <tr>
/// <td align="left" valign="top"><tt>CellEnergyThreshold</tt></td>
/// <td align="center" valign="top"><tt>double</tt></td>
/// <td align="center" valign="top"><tt>m_energyThresholdDef</tt></td>
/// <td align="left" valign="top">cell energy threshold used in exclusive mode only. See further comments below.</td>
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>CellContainerKey</tt></td>
/// <td align="center" valign="top"><tt>SG::ReadHandleKey<CaloCellContainer></tt></td>
/// <td align="center" valign="top"><tt>"AllCalo"</tt></td>
/// <td align="left" valign="top">cell container key is needed to pick up @c CaloCellContainer for all operational modes.</td> 
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>ClusterContainerKey</tt></td>
/// <td align="center" valign="top"><tt>SG::ReadHandleKey<xAOD::CaloClusterContainer></tt></td>
/// <td align="center" valign="top"><tt>"CaloTopoClusters"</tt></td>
/// <td align="left" valign="top">cluster container key is needed to pick up @c xAOD::CaloClusterContainer for filtered mode (<tt>UseCellsFromCluster = true</tt>)
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>CellClusterWeightKey</tt></td>
/// <td align="center" valign="top"><tt>SG::WriteHandleKey<CaloCellClusterWeights></tt></td>
/// <td align="center" valign="top"><tt>&minus;N/A&minus;</tt></td>
/// <td align="left" valign="top">key for @c CaloCellClusterWeights object is needed if <tt>PrepareLCW = true</tt>. Default is empty key. 
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>BuildCombinedTopoSignal</tt></td>
/// <td align="center" valign="top"><tt>bool</tt></td>
/// <td align="center" valign="top"><tt>false</tt></td>
/// <td align="left" valign="top">turns on combined topo-cluster/topo-tower output, with topo-clusters used within the rapidity range defined by <tt>TopoClusterRange</tt> and topo-towers elsewhere.</td></tr>
/// </tr>
/// <tr>
/// <td align="left" valign="top"><tt>TopoClusterRange</tt></td>
/// <td align="center" valign="top"><tt>double</tt></td>
/// <td align="center" valign="top"><tt>5.</tt></td>
/// <td align="left" valign="top">sets the range @f$ y_{\rm topo-cluster}^{\rm max} @f$ for using topo-clusters when <tt>BuildCombinedTopoSignal = true</tt>; 
/// topo-clusters with @f$ \left|y_{\rm topo-cluster}\right| < y_{\rm topo-cluster}^{\rm max} @f$ are used. 
/// </tr>
/// </table> 
///                                 
/// The towers can be classified as:
/// -# <b>inclusive cell towers</b>
///    All cells are collected into inclusive towers, independent of their signal. Requires properties <tt>UseCellsFromClusters = false</tt> and <tt>UseCellEnergyThreshold = false</tt>. Only EM
///    towers are possible, as cells not collected into topo-clustersdo not have a calibration applied.
/// -# <b>exclusive cell towers</b> 
///    Cells with @f$ E > E_{\rm min} @f$ are collected into exclusive towers. This behaviour is turned on by  <tt>UseCellsFromClusters = false</tt> and <tt>UseCellEnergyThreshold = true</tt>. A
///    meaningful <tt>CellEnergyThreshold</tt> value needs to be provided in addition.
/// -# <b>filtered mode</b>
///    Cells contributing to standard topo-clusters are collected into topo-towers. This behaviour is triggered by <tt>UseCellsFromClusters = true</tt>. Optionally, LCW calibration can be applied
///    to these towers by setting <tt>PrepareLCW = true</tt> and scheduling a @c CaloTopoClusterFromTowerCalibrator tool after the cluster moment calculators. The values of the <tt>UseEnergyThreshold</tt>
///    and <tt>CellEnergyThreshold</tt> properties are ignored in this mode. A valid event store key needs to be provided in the to pick up the topo-cluster container. Note that building EM 
///    topo-towers requires topo-clusters on EM scale (no LCW applied) to get the correct geometrical cell weights only. LCW topo-towers require LCW scale topo-clusters to get the correct full geometrical 
///    and calibration weights.
/// -# <b>mixed mode</b>
///    Cells contributing to standard topo-clusters are collected into towers if these topo-clusters are outside of a give rapidity range. The rapidity range is defined by the <tt>TopoClusterRange</tt>
///    property. This mode is turned on by setting the property <tt>BuildCombinedTopoSignal = true</tt>. It is turned off by default (<tt>BuildCombinedTopoSignal = false</tt>). 
///    EM scale and LCW scale is possible, as in the filtered mode. 
///
///  Configuration 2 and 3 are exclusive, with 3 overwriting 2. The output topo-clusters represent calorimeter towers on the EM scale. The can be handed to cluster moment
///  tools (needs EM scale) and, if desired, to a dedicated cluster calibration tool of type @c xAOD::CaloTowerClusterFromTowerCalibrator .  
///
///  To avoid multiple retrievals of the same weights by searching for cells in (many) topo-clusters, the overall weight of the cell signal is stored in a random access
///  look-up table stored in a @c CaloCellClusterWeights object in the detector store. This object is created by this tool, if needed. If the tool property 
///  @c CellWeightLookupKey is set, this object will be generated, filled, and recorded. This is essential for calibrated topo-towers!
///
///@note The @c OrderByPt property, which orders the container by descending transverse momentum, is only useful for EM towers. Applyin LCW may lead to a different 
///      order - if a container with LCW towers should be ordered, the corresponding property @c OrderByPt of the @c CaloTopoClusterFromTowerCalibrator tool should
///      be set to @c true. 
///
///@note Many more details on the towers are available on 
///      <a href="https://twiki.cern.ch/twiki/bin/view/AtlasSandboxProtected/CaloTowerPerformance" title="https://twiki.cern.ch/twiki/bin/view/AtlasSandboxProtected/CaloTowerPerformance">this page</a>.
///
/// @author Peter Loch <loch@physics.arizona.edu>
#endif
