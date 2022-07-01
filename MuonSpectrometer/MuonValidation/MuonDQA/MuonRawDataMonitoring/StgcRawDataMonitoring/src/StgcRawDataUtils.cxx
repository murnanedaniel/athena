/*                                                                             
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration                              
*/

///////////////////////////////////////////////////////////////////////////                          
// Utils for the main sTGCRawDataMonAlg.cxx                                                            
// Part of StgcRawDataMonAlg.h                                                                         
// see StgcRawDataMonAlg.cxx                                
///////////////////////////////////////////////////////////////////////////                      
         
#include "StgcRawDataMonitoring/StgcRawDataMonAlg.h"

int sTgcRawDataMonAlg::get_sectorPhi_from_stationPhi_stName(int stationPhi, const std::string &stName) const {

  
  if (stationPhi==1 && stName=="MML") return 1;
  if (stationPhi==1 && stName=="MMS") return 2;
  if (stationPhi==2 && stName=="MML") return 3;
  if (stationPhi==2 && stName=="MMS") return 4;
  if (stationPhi==3 && stName=="MML") return 5;
  if (stationPhi==3 && stName=="MMS") return 6;
  if (stationPhi==4 && stName=="MML") return 7;
  if (stationPhi==4 && stName=="MMS") return 8;
  if (stationPhi==5 && stName=="MML") return 9;
  if (stationPhi==5 && stName=="MMS") return 10;
  if (stationPhi==6 && stName=="MML") return 11;
  if (stationPhi==6 && stName=="MMS") return 12;
  if (stationPhi==7 && stName=="MML") return 13;
  if (stationPhi==7 && stName=="MMS") return 14;
  if (stationPhi==8 && stName=="MML") return 15;
  if (stationPhi==8 && stName=="MMS") return 16;

  throw std::invalid_argument( "stationPhi and stName are not valid!" );
 
 
}



int sTgcRawDataMonAlg::get_PCB_from_channel(int channel) const {

  if (channel>0 && channel<=1024) return 1; 
  if (channel>1024 && channel<=2048) return 2;
  if (channel>2048 && channel<=3072) return 3;
  if (channel>3072 && channel<=4096) return 4;
  if (channel>4096 && channel<=5120) return 5;

  throw std::invalid_argument( "channel is not valid!" );
}


int sTgcRawDataMonAlg::get_bin_for_occ_CSide_pcb_eta1_hist(int stationEta, int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 5;
  static const int max_gas_gap = 4;
  if (stationEta != -1) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb+ (gas_gap-1)*max_pcb + (PCB-1);

}






