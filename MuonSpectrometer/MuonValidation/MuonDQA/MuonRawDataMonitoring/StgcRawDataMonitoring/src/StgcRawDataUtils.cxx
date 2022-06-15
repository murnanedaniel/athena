/*                                                                             
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration                              
*/

///////////////////////////////////////////////////////////////////////////                          
// Utils for the main sTGCRawDataMonAlg.cxx                                                            
// Part of StgcRawDataMonAlg.h                                                                         
// see StgcRawDataMonAlg.cxx                                
///////////////////////////////////////////////////////////////////////////                      
         
#include "StgcRawDataMonitoring/StgcRawDataMonAlg.h"
#include <string>
#include <stdexcept>

int sTgcRawDataMonAlg::get_sectorPhi_from_stationPhi_stName(int stationPhi, const std::string &stName) const 
{  
  if (stationPhi == 1 && stName == "STS") return 1;
  if (stationPhi == 1 && stName == "STL") return 2;
  if (stationPhi == 2 && stName == "STS") return 3;
  if (stationPhi == 2 && stName == "STL") return 4;
  if (stationPhi == 3 && stName == "STS") return 5;
  if (stationPhi == 3 && stName == "STL") return 6;
  if (stationPhi == 4 && stName == "STS") return 7;
  if (stationPhi == 4 && stName == "STL") return 8;
  if (stationPhi == 5 && stName == "STS") return 9;
  if (stationPhi == 5 && stName == "STL") return 10;
  if (stationPhi == 6 && stName == "STS") return 11;
  if (stationPhi == 6 && stName == "STL") return 12;
  if (stationPhi == 7 && stName == "STS") return 13;
  if (stationPhi == 7 && stName == "STL") return 14;
  if (stationPhi == 8 && stName == "STS") return 15;
  if (stationPhi == 8 && stName == "STL") return 16;

  throw std::invalid_argument("stationPhi and stName are not valid!");
}

