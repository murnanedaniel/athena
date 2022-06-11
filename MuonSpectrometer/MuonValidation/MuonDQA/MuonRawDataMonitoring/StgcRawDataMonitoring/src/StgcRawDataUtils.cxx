/*                                                                                                        													  
Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration                                   
*/
///////////////////////////////////////////////////////////////////////////                                  
//Utils for the main sTGCRawDataMonAlg.cxx                                                                    
        
//Part of StgcRawDataMonAlg.h                                                                                 
         
//Authors                                                                                                   
         
//   see StgcRawDataMonAlg.cxx                                                                                
        
///////////////////////////////////////////////////////////////////////////                                 
         

#include "StgcRawDataMonitoring/StgcRawDataMonAlg.h"
#include <string>
#include <stdexcept>

int StgcRawDataMonAlg::get_PCB_from_channel(int channel) const {

  if (channel>0 && channel<=1024) return 1;
  if (channel>1024 && channel<=2048) return 2;
  if (channel>2048 && channel<=3072) return 3;
  if (channel>3072 && channel<=4096) return 4;
  if (channel>4096 && channel<=5120) return 5;

  throw std::invalid_argument( "channel is not valid!" );
}

int StgcRawDataMonAlg::get_sectorPhi_from_stationPhi_stName(int stationPhi,const std::string & stName) const {
  
  if (stationPhi==1 && stName=="S") return 1;
  if (stationPhi==1 && stName=="L") return 2;
  if (stationPhi==2 && stName=="S") return 3;
  if (stationPhi==2 && stName=="L") return 4;
  if (stationPhi==3 && stName=="S") return 5;
  if (stationPhi==3 && stName=="L") return 6;
  if (stationPhi==4 && stName=="S") return 7;
  if (stationPhi==4 && stName=="L") return 8;
  if (stationPhi==5 && stName=="S") return 9;
  if (stationPhi==5 && stName=="L") return 10;
  if (stationPhi==6 && stName=="S") return 11;
  if (stationPhi==6 && stName=="L") return 12;
  if (stationPhi==7 && stName=="S") return 13;
  if (stationPhi==7 && stName=="L") return 14;
  if (stationPhi==8 && stName=="S") return 15;
  if (stationPhi==8 && stName=="L") return 16;

  throw std::invalid_argument( "stationPhi and stName are not valid!" );

}

int StgcRawDataMonAlg::get_sectorEta_from_stationEta(int stationEta) const {
 
  //  1<-0  0-> 1
  if (std::abs(stationEta)==1) return 0;                                                            
  if (std::abs(stationEta)==2) return 1;

  return -1;

}

int StgcRawDataMonAlg::get_bin_for_occ_CSide_hist(int stationEta, int multiplet, int gas_gap) const {

  static const int max_gas_gap = 4;
  static const int max_multiplet = 2;

  return (stationEta+2)*(max_gas_gap*max_multiplet)+(multiplet-1)*max_gas_gap +(gas_gap-1);

}

int StgcRawDataMonAlg::get_bin_for_occ_ASide_hist(int stationEta, int multiplet, int gas_gap) const {

  static const int max_gas_gap = 4;
  static const int max_multiplet = 2;

  return (stationEta-1)*(max_gas_gap*max_multiplet)+(multiplet-1)*max_gas_gap +(gas_gap-1);

}


int StgcRawDataMonAlg::get_bin_for_occ_CSide_pcb_eta2_hist(int stationEta, int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 3;
  static const int max_gas_gap = 4;
  if (stationEta != -2) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb+ (gas_gap-1)*max_pcb + (PCB-1);

}

int StgcRawDataMonAlg::get_bin_for_occ_CSide_pcb_eta1_hist(int stationEta, int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 5;
  static const int max_gas_gap = 4;
  if (stationEta != -1) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb+ (gas_gap-1)*max_pcb + (PCB-1);

}


int StgcRawDataMonAlg::get_bin_for_occ_ASide_pcb_eta2_hist(int stationEta, int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 3;
  static const int max_gas_gap = 4;
  if (stationEta != 2) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb+ (gas_gap-1)*max_pcb + (PCB-1);

}


int StgcRawDataMonAlg::get_bin_for_occ_ASide_pcb_eta1_hist(int stationEta, int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 5;
  static const int max_gas_gap = 4;
  if (stationEta != 1) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb+ (gas_gap-1)*max_pcb + (PCB-1);
}


int StgcRawDataMonAlg::get_bin_for_occ_lb_CSide_pcb_eta2_hist(int stationEta, int multiplet, int gas_gap, int PCB,int isector) const {

  static const int max_pcb = 3;
  static const int max_gas_gap = 4;
  static const int max_isector = 2;
  if (stationEta != -2) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb*max_isector+ (gas_gap-1)*max_pcb*max_isector + isector*max_pcb+ (PCB-1);

}


int StgcRawDataMonAlg::get_bin_for_occ_lb_CSide_pcb_eta1_hist(int stationEta, int multiplet, int gas_gap, int PCB,int isector) const {

  static const int max_pcb = 5;
  static const int max_gas_gap = 4;
  static const int max_isector = 2;
  if (stationEta != -1) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb*max_isector+ (gas_gap-1)*max_pcb*max_isector + isector*max_pcb + (PCB-1);

}

int StgcRawDataMonAlg::get_bin_for_occ_lb_ASide_pcb_eta1_hist(int stationEta, int multiplet, int gas_gap, int PCB,int isector) const {

  static const int max_pcb = 5;
  static const int max_gas_gap = 4;
  static const int max_isector = 2;
  if (stationEta != 1) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb*max_isector+ (gas_gap-1)*max_pcb*max_isector + isector*max_pcb + (PCB-1);
}


int StgcRawDataMonAlg::get_bin_for_occ_lb_ASide_pcb_eta2_hist(int stationEta, int multiplet, int gas_gap, int PCB,int isector) const {

  static const int max_pcb = 3;
  static const int max_gas_gap = 4;
  static const int max_isector = 2;
  if (stationEta != 2) return -1;

  return  (multiplet-1)*max_gas_gap*max_pcb*max_isector+ (gas_gap-1)*max_pcb*max_isector + isector*max_pcb + (PCB-1);

}

int StgcRawDataMonAlg::get_bin_for_occ_lb_pcb_hist(int multiplet, int gas_gap, int PCB) const {

  static const int max_pcb = 8;
  static const int max_gas_gap = 4;
  
  return (multiplet-1)*max_gas_gap*max_pcb + (gas_gap-1)*max_pcb + (PCB-1);

}
