/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "AtlasDetectorIDHelper.h"
#include "IdDict/IdDictDefs.h"  
#include "AtlasDetDescr/AtlasDetectorID.h"
#include <iostream>

AtlasDetectorIDHelper::AtlasDetectorIDHelper(IMessageSvc* msgSvc) :
    AthMessaging(msgSvc, "AtlasDetectorIDHelper")
{
}

int         
AtlasDetectorIDHelper::initialize_from_dictionary(const IdDictMgr& dict_mgr,
                                                  bool quiet)
{

    if(m_initialized) return(0);
    m_initialized = true;

    AtlasDetectorID atlas_id;
    ExpandedIdentifier id;

    const IdDictDictionary* 	dict = dict_mgr.find_dictionary ("InnerDetector"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access InnerDetector dictionary");
		return 1;
    }
    else {
	// Check if this is High Luminosity LHC layout
    if (dict->m_version=="ITkHGTD" || dict->m_version=="ITkHGTDPLR" || dict->m_version=="P2-RUN4") {
      m_isHighLuminosityLHC = true;
    }

	// Save index to a PIXEL region for unpacking
	id = atlas_id.pixel_exp(); 
	if (dict->find_region(id, m_pixel_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find pixel region index: id, reg "
                        << (std::string)id << " " << m_pixel_region_index);
	}

	// for High Luminosity LHC layout one cannot get the sct region as below, nor
	// is there any trt regions
	if (!m_isHighLuminosityLHC) {
	    
	    // Save index to a SCT region for unpacking
	    id = atlas_id.sct_exp();
	    if (dict->find_region(id, m_sct_region_index)) {
            ATH_MSG_WARNING("initialize_from_dictionary - unable to find sct region index: id, reg "
                            << (std::string)id << " " << m_sct_region_index);
	    }

	    // Save index to a TRT region for unpacking
	    id = atlas_id.trt_exp(); 
	    if (dict->find_region(id, m_trt_region_index)) {
            ATH_MSG_WARNING("initialize_from_dictionary - unable to find trt region index: id, reg "
                            << (std::string)id << " " << m_trt_region_index);
        }
    }

    dict = dict_mgr.find_dictionary ("LArCalorimeter"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access LArCalorimeter dictionary");
		return 1;
    }
    else {
	
	// Save index to a LAR_EM region for unpacking
	id = atlas_id.lar_em_exp(); 
	if (dict->find_region(id, m_lar_em_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find lar_em region index: id, reg "
                        << (std::string)id << " " << m_lar_em_region_index);
	}

	// Save index to a LAR_HEC region for unpacking
	id = atlas_id.lar_hec_exp(); 
	if (dict->find_region(id, m_lar_hec_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find lar_hec region index: id, reg "
                        << (std::string)id << " " << m_lar_hec_region_index);
	}

	// Save index to a LAR_FCAL region for unpacking
	id = atlas_id.lar_fcal_exp(); 
	if (dict->find_region(id, m_lar_fcal_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find lar_fcal region index: id, reg "
                        << (std::string)id << " " << m_lar_fcal_region_index);
    }
    }
    
    // Get Calorimetry dictionary for both LVL1 and Dead material
    dict = dict_mgr.find_dictionary ("Calorimeter"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access Calorimeter dictionary");
		return 1;
    }
    else {
	
	// Save index to a LVL1 region for unpacking
	IdDictRegion* reg = dict->find_region("Lvl1_0");
	if (reg) {
	    m_lvl1_region_index = reg->m_index;
	}
	else {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find lvl1 region");
	}
    
	// Save index to a Dead Material region for unpacking
	reg = dict->find_region("DM_4_1_0_0");
	if (reg) {
	    m_dm_region_index = reg->m_index;
	}
	else {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find dead material region");
    }
    }
    
    dict = dict_mgr.find_dictionary ("TileCalorimeter"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access TileCalorimeter dictionary");
		return 1;
    }
    else {
	
	// Save index to a TILE region for unpacking
	id = atlas_id.tile_exp(); 
	if (dict->find_region(id, m_tile_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find tile region index: id, reg");
	}
    }
    
    dict = dict_mgr.find_dictionary ("MuonSpectrometer"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access MuonSpectrometer dictionary");
		return 1;
    }
    else {
	
	m_station_field = dict->find_field ("stationName");
	if(!m_station_field) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access stationName field");
	}
	else {
	    m_muon_station_index = m_station_field->m_index;
	}
    }
    

	// Save index to a MDT region for unpacking
	IdDictGroup* group = dict->find_group ("mdt");
	if (group) {
	    const std::vector<IdDictRegion*>&    regions = group->regions();
	    if (regions.size() > 0) {
		m_mdt_region_index = regions[0]->m_index;
	    }
	}
    

	if (UNDEFINED == m_mdt_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find mdt region index: group, region size "
                        << group << " " << size);
	}

	// Save index to a CSC region for unpacking
	group = dict->find_group ("csc");
	if (group) {
	    if (group->regions().size() > 0) {
		m_csc_region_index = group->regions()[0]->m_index;
	    }
	}
	if (UNDEFINED == m_csc_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
        ATH_MSG_DEBUG("initialize_from_dictionary - unable to find csc region index: group, region size "
                      << group << " " << size);
    }

	// Save index to a RPC region for unpacking
	group = dict->find_group ("rpc");
	if (group) {
	    if (group->regions().size() > 0) {
		m_rpc_region_index = group->regions()[0]->m_index;
	    }
	}
	if (UNDEFINED == m_rpc_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find rpc region index: group, region size "
                        << group << " " << size);
	}

	// Save index to a TGC region for unpacking
	group = dict->find_group ("tgc");
	if (group) {
	    if (group->regions().size() > 0) {
		m_tgc_region_index = group->regions()[0]->m_index;
	    }
	}
	if (UNDEFINED == m_tgc_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find tgc region index: group, region size "
                        << group << " " << size);
	}

	// Save index to a MM region for unpacking
	group = dict->find_group ("mm");
	if (group) {
	    if (group->regions().size() > 0) {
		m_mm_region_index = group->regions()[0]->m_index;
	    }
	}
	if (UNDEFINED == m_mm_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
            if (!quiet) {
                ATH_MSG_DEBUG("initialize_from_dictionary - unable to find mm region index: group, region size "
                              << group << " " << size);
            }
	}

	// Save index to a sTGC region for unpacking
	group = dict->find_group ("stgc");
	if (group) {
	    if (group->regions().size() > 0) {
		m_stgc_region_index = group->regions()[0]->m_index;
	    }
	}
	if (UNDEFINED == m_stgc_region_index) {
	    int size = 0;
	    if (group) {
	        size = group->regions().size();
	    }
            if (!quiet) {
                ATH_MSG_DEBUG("initialize_from_dictionary - unable to find stgc region index: group, region size "
                              << group << " " << size);
            }
	}
    }


    dict = dict_mgr.find_dictionary ("ForwardDetectors"); 
    if(!dict) {
        ATH_MSG_WARNING("initialize_from_dictionary - cannot access ForwardDetectors dictionary");
		return 1;
    }
    else {
	
	// Save index to a ALFA region for unpacking
	id = atlas_id.alfa_exp(); 
	if (dict->find_region(id, m_alfa_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find alfa region index: id, reg "
                        << (std::string)id << " " << m_alfa_region_index);
	}

	// Save index to a BCM region for unpacking
	id = atlas_id.bcm_exp(); 
	if (dict->find_region(id, m_bcm_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find bcm region index: id, reg "
                        << (std::string)id << " " << m_bcm_region_index);
	}

	// Save index to a LUCID region for unpacking
	id = atlas_id.lucid_exp(); 
	if (dict->find_region(id, m_lucid_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find lucid region index: id, reg "
                        << (std::string)id << " " << m_lucid_region_index);
	}

	// Save index to a ZDC region for unpacking
	id = atlas_id.zdc_exp(); 
	if (dict->find_region(id, m_zdc_region_index)) {
        ATH_MSG_WARNING("initialize_from_dictionary - unable to find zdc region index: id, reg "
                        << (std::string)id << " " << m_zdc_region_index);
	}
    }
    
//      std::cout << "AtlasDetectorIDHelper::initialize_from_dictionary " << std::endl;
//      std::cout << " pixel_region_index     " << m_pixel_region_index      << std::endl;
//      std::cout << " sct_region_index       " << m_sct_region_index        << std::endl;
//      std::cout << " trt_region_index       " << m_trt_region_index        << std::endl;
//      std::cout << " lar_em_region_index    " << m_lar_em_region_index     << std::endl;
//      std::cout << " lar_hec_region_index   " << m_lar_hec_region_index    << std::endl;
//      std::cout << " lar_fcal_region_index  " << m_lar_fcal_region_index   << std::endl;
//      std::cout << " lvl1_region_index      " << m_lvl1_region_index       << std::endl;
//      std::cout << " tile_region_index      " << m_tile_region_index       << std::endl;
//      std::cout << " mdt_region_index       " << m_mdt_region_index        << std::endl;
//      std::cout << " csc_region_index       " << m_csc_region_index        << std::endl;
//      std::cout << " rpc_region_index       " << m_rpc_region_index        << std::endl;
//      std::cout << " tgc_region_index	  " << m_tgc_region_index        << std::endl;
//      std::cout << " muon_station_index	  " << m_muon_station_index      << std::endl;

    return (0);
}
