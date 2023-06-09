/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//Dear emacs, this is -*-c++-*-

#ifndef LARCABLING_LARSUPERCELLCABLINGTOOL_H
#define LARCABLING_LARSUPERCELLCABLINGTOOL_H

#include "LArCabling/LArCablingService.h"


/** 
   @class LArSuperCellCablingTool 
   @brief Tool for mapping online and offline identifiers, <br>
   @author Walter Lampl
*/


static const InterfaceID IID_LArSuperCellCablingTool("LArSuperCellCablingTool", 1 , 0); 


class LArSuperCellCablingTool : public LArCablingBase  {

public:
   /** constructor */
  LArSuperCellCablingTool( const std::string& type, const std::string& name, const IInterface* parent ) ;

  ~LArSuperCellCablingTool();
  
  StatusCode initialize( ); 

  StatusCode iovCallBack(IOVSVC_CALLBACK_ARGS_K(keys));

  // Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_LArSuperCellCablingTool; }


  //All other methods inherited from regular LArCablingBase
};

#endif // LARCABLING_LARSUPERCELLCABLINGTOOL_H
