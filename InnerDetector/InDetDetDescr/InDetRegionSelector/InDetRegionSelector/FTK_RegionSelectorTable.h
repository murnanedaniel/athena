// emacs: this is -*- c++ -*-
//
//   FTK_RegionSelectorTable.h        
//
//   Create the new RegSelSiLUT for the FTK
// 
//   Copyright (C) 2011 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: FTK_RegionSelectorTable.h, v0.0 Fri 24 Jun 2011 13:10:59 BST $

 
#ifndef InDetRegionSelector_FTK_RegionSelectorTable_h
#define InDetRegionSelector_FTK_RegionSelectorTable_h

#include "RegSelLUT/IRegionFTKLUT_Creator.h"

// #include "GaudiKernel/AlgTool.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"

// class StoreGateSvc;
class RegSelSiLUT;

#include <string>
using std::string;

/////////////////////////////////////////////////////////////////////////////

class FTK_RegionSelectorTable : public AthAlgTool, virtual public IRegionFTKLUT_Creator
{

public:
  FTK_RegionSelectorTable (const std::string&, 
			 const std::string&,
			 const IInterface*);
  virtual ~FTK_RegionSelectorTable();
  StatusCode initialize();
  StatusCode finalize();
  
  virtual RegSelEtaPhiLUT* getLUT() const;

private:
  
  StatusCode createTable();
  
  //  StoreGateSvc*  m_detStore;
  RegSelEtaPhiLUT*   m_regionLUT;

  // Algorithm properties
  std::string m_managerName;
  double m_deltaZ;
  std::string m_roiFileName;
  bool m_printHashId;
  bool m_printTable;

  // cablings

};

#endif // InDetRegionSelector_FTK_RegionSelectorTable_h
