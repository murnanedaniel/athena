/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "GaudiKernel/MsgStream.h"

#include "HGTD_Identifier/HGTD_ID.h"
#include "Identifier/IdentifierHash.h"
#include "IdDict/IdDictDefs.h"  
#include <set>
#include <algorithm>
#include <iostream>

// Constructor
HGTD_ID::HGTD_ID(void):
    m_hgtd_region_index(0),
    m_INDET_INDEX(0),
    m_HGTD_INDEX(1),
    m_ENDCAP_INDEX(2),
    m_LAYER_INDEX(3),
    m_PHI_MODULE_INDEX(4),
    m_ETA_MODULE_INDEX(5),
    m_PHI_INDEX_INDEX(6),
    m_ETA_INDEX_INDEX(7),
    m_dict(0),
    m_wafer_hash_max(0),
    m_pixel_hash_max(0) {

}

// Destructor
HGTD_ID::~HGTD_ID(void)
{}

void
HGTD_ID::wafer_id_checks ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module ) const
{

    // Check that id is within allowed range

    // Fill expanded id
    ExpandedIdentifier id;
    id << indet_field_value() << hgtd_field_value()
       << endcap << layer << phi_module << eta_module;
    if (!m_full_wafer_range.match(id)) {  // module range check is sufficient
        MsgStream log(m_msgSvc, "HGTD_ID");
        if(m_msgSvc) log << MSG::ERROR << " HGTD_ID::wafer_id result is NOT ok. ID, range " 
                         << (std::string)id << "   " << (std::string)m_full_wafer_range << endmsg;
        else std::cout << " ERROR HGTD_ID::wafer_id result is NOT ok. ID, range " 
                       << (std::string)id << "   " << (std::string)m_full_wafer_range << std::endl;
    }
}

void
HGTD_ID::pixel_id_checks ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module,
                           int phi_index,
                           int eta_index) const
{

    // Check that id is within allowed range

    // Fill expanded id
    ExpandedIdentifier id;
    id << indet_field_value() << hgtd_field_value()
       << endcap << layer << phi_module << eta_module << phi_index << eta_index;

    if (!m_full_pixel_range.match(id)) {
        MsgStream log(m_msgSvc, "HGTD_ID");
        if(m_msgSvc) log << MSG::ERROR << " HGTD_ID::pixel_id result is NOT ok. ID, range " 
                         << (std::string)id << " " << (std::string)m_full_pixel_range << endmsg;
        else std::cout << " ERROR HGTD_ID::pixel_id result is NOT ok. ID, range " 
                       << (std::string)id << " " << (std::string)m_full_pixel_range << std::endl;
    }
}

int
HGTD_ID::layer_max(const Identifier& id) const
{
    // get max from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& layer_field = range[m_LAYER_INDEX];
            if (layer_field.has_maximum()) {
                return (layer_field.get_maximum());
            }
        }
    }
    return (-999);  // default
}

int
HGTD_ID::phi_module_max(const Identifier& id) const
{
    // get max from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& phi_module_field = range[m_PHI_MODULE_INDEX];
            if (phi_module_field.has_maximum()) {
                return (phi_module_field.get_maximum());
            }
        }
    }
    return (-999);  // default
}

int 
HGTD_ID::eta_module_max(const Identifier& id) const
{
    // get max from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& eta_module_field = range[m_ETA_MODULE_INDEX];
            if (eta_module_field.has_maximum()) {
                return (eta_module_field.get_maximum());
            }
        }
    }
    return (-999); // default
}

int
HGTD_ID::eta_module_min(const Identifier& id) const
{
    // get min from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& eta_module_field = range[m_ETA_MODULE_INDEX];
            if (eta_module_field.has_minimum()) {
                return (eta_module_field.get_minimum());
            }
        }
    }
    return (-999); // default
}

int
HGTD_ID::phi_index_max(const Identifier& id) const
{
    // get max from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& phi_field = range[m_PHI_INDEX_INDEX];
            if (phi_field.has_maximum()) {
                return (phi_field.get_maximum());
            }
        }
    }
    return (-999);  // default
}

int
HGTD_ID::eta_index_max(const Identifier& id) const
{
    // get max from dictionary
    ExpandedIdentifier expId;
    IdContext wafer_context1 = wafer_context();
    get_expanded_id(id, expId, &wafer_context1);
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        if (range.match(expId)) {
            const Range::field& eta_field = range[m_ETA_INDEX_INDEX];
            if (eta_field.has_maximum()) {
                return (eta_field.get_maximum());
            }
        }
    }
    return (-999);  // default
}

bool
HGTD_ID::is_phi_module_max(const Identifier& id) const
{
    return (phi_module(id) == phi_module_max(id));
}

int
HGTD_ID::initialize_from_dictionary(const IdDictMgr& dict_mgr)
{

    MsgStream log(m_msgSvc, "HGTD_ID");
    if(m_msgSvc) log << MSG::INFO << "Initialize from dictionary" << endmsg;
    else         std::cout << " INFO Initialize from dictionary" << std::endl;

    // Check whether this helper should be reinitialized
    if (!reinitialize(dict_mgr)) {
        if(m_msgSvc) log << MSG::INFO << "Request to reinitialize not satisfied - tags have not changed" << endmsg;
        else std::cout << " INFO Request to reinitialize not satisfied - tags have not changed" << std::endl;
        return (0);
    }
    else {
        if(m_msgSvc) log << MSG::DEBUG << "(Re)initialize" << endmsg;
        else         std::cout << " DEBUG (Re)initialize" << std::endl;
    }

    // init base object
    if(AtlasDetectorID::initialize_from_dictionary(dict_mgr)) return (1);

    // Register version of InnerDetector dictionary 
    if (register_dict_tag(dict_mgr, "InnerDetector")) return(1); 

    m_dict = dict_mgr.find_dictionary ("InnerDetector"); 
    if(!m_dict) {
        if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::initialize_from_dict - cannot access InnerDetector dictionary "
                         << endmsg;
        else std::cout << " FATAL HGTD_ID::initialize_from_dict - cannot access InnerDetector dictionary "
                       << std::endl;
        return (1);
    }

    AtlasDetectorID::setDictVersion(dict_mgr, "InnerDetector");

    // Initialize the field indices
    if(initLevelsFromDict()) return (1);

    //
    // Build multirange for the valid set of identifiers
    //

    // Find value for the field InnerDetector
    const IdDictDictionary* atlasDict = dict_mgr.find_dictionary ("ATLAS"); 
    int inDetField   = -1;
    if (atlasDict->get_label_value("subdet", "InnerDetector", inDetField)) {
        if(m_msgSvc) log << MSG::FATAL << "Could not get value for label 'InnerDetector' of field 'subdet' in dictionary " 
                         << atlasDict->m_name
                         << endmsg;
        else std::cout << " FATAL Could not get value for label 'InnerDetector' of field 'subdet' in dictionary " 
                       << atlasDict->m_name
                       << std::endl;
        return (1);
    }

    // Find value for the field HGTD
    int hgtdField   = -1;
    if (m_dict->get_label_value("part", "HGTD", hgtdField)) {
        if(m_msgSvc) log << MSG::FATAL << "Could not get value for label 'HGTD' of field 'part' in dictionary " 
                         << m_dict->m_name
                         << endmsg;
        else std::cout << " FATAL Could not get value for label 'HGTD' of field 'part' in dictionary " 
                       << m_dict->m_name
                       << std::endl;
        return (1);
    }
    if(m_msgSvc) log << MSG::DEBUG << " HGTD_ID::initialize_from_dict " 
                     << "Found field values: InDet/HGTD "  
                     << inDetField << "/"
                     << hgtdField
                     << endmsg;
    else std::cout << " DEBUG HGTD_ID::initialize_from_dict " 
                   << "Found field values: InDet/HGTD "  
                   << inDetField << "/"
                   << hgtdField
                   << std::endl;

    // Set up id for region and range prefix
    ExpandedIdentifier region_id;
    region_id.add(inDetField);
    region_id.add(hgtdField);
    Range prefix;
    m_full_wafer_range = m_dict->build_multirange(region_id, prefix, m_moduleInRow);
    m_full_pixel_range = m_dict->build_multirange(region_id, prefix);

    // Setup the hash tables
    if(init_hashes()) return (1);

    // Setup hash tables for finding neighbors
    if(init_neighbors()) return (1);

    if(m_msgSvc) {
        log << MSG::DEBUG << " HGTD_ID::initialize_from_dict " 
            << endmsg;
        log << MSG::DEBUG  
            << "Wafer range -> " << (std::string)m_full_wafer_range
            <<  endmsg;
        log << MSG::DEBUG 
            << "Pixel range -> " << (std::string)m_full_pixel_range
            << endmsg;
    }
    else {
        std::cout << " DEBUG HGTD_ID::initialize_from_dict " 
                  << std::endl;
        std::cout << " DEBUG Wafer range -> " << (std::string)m_full_wafer_range
                  <<  std::endl;
        std::cout << " DEBUG Pixel range -> " << (std::string)m_full_pixel_range
                  << std::endl;
    }
    return 0;
}

int
HGTD_ID::init_hashes(void)
{

    //
    // create a vector(s) to retrieve the hashes for compact ids. For
    // the moment, we implement a hash for wafers but NOT for pixels
    // (too many)
    //

    MsgStream log(m_msgSvc, "HGTD_ID");

    // wafer hash
    m_wafer_hash_max = m_full_wafer_range.cardinality();
    m_wafer_vec.clear();
    m_wafer_vec.resize(m_wafer_hash_max);
    unsigned int nids = 0;
    std::set<Identifier> ids;
    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        Range::const_identifier_factory first = range.factory_begin();
        Range::const_identifier_factory last  = range.factory_end();
        for (; first != last; ++first) {
            const ExpandedIdentifier& exp_id = (*first);
            Identifier id = wafer_id (exp_id[m_ENDCAP_INDEX],
                                      exp_id[m_LAYER_INDEX],
                                      exp_id[m_PHI_MODULE_INDEX],
                                      exp_id[m_ETA_MODULE_INDEX]);
            if(!(ids.insert(id)).second) {
                if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_hashes "
                                 << " Error: duplicated id for wafer id. nid " << nids
                                 << " id " << show_to_string(id)
                                 << " exp id " << (std::string)exp_id 
                                 << " " << (std::string)m_full_wafer_range << endmsg;
                else std::cout << " FATAL HGTD_ID::init_hashes "
                               << " Error: duplicated id for wafer id. nid " << nids
                               << " id " << show_to_string(id)
                               << " exp id " << (std::string)exp_id 
                               << " " << (std::string)m_full_wafer_range << std::endl;
                return (1);
            }
            nids++;
        }
    }
    if(ids.size() != m_wafer_hash_max) {
        if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_hashes "
                         << " Error: set size NOT EQUAL to hash max. size " << ids.size()
                         << " hash max " << m_wafer_hash_max 
                         << endmsg;
        else std::cout << " FATAL HGTD_ID::init_hashes "
                       << " Error: set size NOT EQUAL to hash max. size " << ids.size()
                       << " hash max " << m_wafer_hash_max 
                       << std::endl;
        return (1);
    }

    nids = 0;
    std::set<Identifier>::const_iterator first = ids.begin();
    std::set<Identifier>::const_iterator last  = ids.end();
    for (; first != last && nids < m_wafer_vec.size(); ++first) {
        m_wafer_vec[nids] = (*first);
        nids++;
    }

    // pixel hash - we do not keep a vec for the pixels - too large
    // TODO: is it worthwhile/possible though for HGTD?
    m_pixel_hash_max = m_full_pixel_range.cardinality();

    return 0;

}

int
HGTD_ID::get_prev_in_phi(const IdentifierHash& id, IdentifierHash& prev) const
{
    unsigned short index = id;
    if (index < m_prev_phi_wafer_vec.size()) {
        if (m_prev_phi_wafer_vec[index] == NOT_VALID_HASH) return (1);
        prev =  m_prev_phi_wafer_vec[index];
        return (0);
    }
    return (1);
}

int
HGTD_ID::get_next_in_phi(const IdentifierHash& id, IdentifierHash& next) const
{
    unsigned short index = id;
    if (index < m_next_phi_wafer_vec.size()) {
        if (m_next_phi_wafer_vec[index] == NOT_VALID_HASH) return (1);
        next =  m_next_phi_wafer_vec[index];
        return (0);
    }
    return (1);
}

int
HGTD_ID::get_prev_in_eta(const IdentifierHash& id, IdentifierHash& prev) const
{
    unsigned short index = id;
    if (index < m_prev_eta_wafer_vec.size()) {
        if (m_prev_eta_wafer_vec[index] == NOT_VALID_HASH) return (1);
        prev =  m_prev_eta_wafer_vec[index];
        return (0);
    }
    return (1);
}

int
HGTD_ID::get_next_in_eta(const IdentifierHash& id, IdentifierHash& next) const
{
    unsigned short index = id;
    if (index < m_next_eta_wafer_vec.size()) {
        if (m_next_eta_wafer_vec[index] == NOT_VALID_HASH) return (1);
        next =  m_next_eta_wafer_vec[index];
        return (0);
    }
    return (1);
}

int
HGTD_ID::init_neighbors(void)
{
    //
    // create a vector(s) to retrieve the hashes for compact ids for
    // wafer neighbors.
    //
    MsgStream log(m_msgSvc, "HGTD_ID");

    if(m_msgSvc) log << MSG::DEBUG << "HGTD_ID::init_neighbors " << endmsg;
    else std::cout << " DEBUG HGTD_ID::init_neighbors " << std::endl;

    m_prev_phi_wafer_vec.clear();
    m_next_phi_wafer_vec.clear();
    m_prev_eta_wafer_vec.clear();
    m_next_eta_wafer_vec.clear();
    m_prev_phi_wafer_vec.resize(m_wafer_hash_max, NOT_VALID_HASH);
    m_next_phi_wafer_vec.resize(m_wafer_hash_max, NOT_VALID_HASH);
    m_prev_eta_wafer_vec.resize(m_wafer_hash_max, NOT_VALID_HASH);
    m_next_eta_wafer_vec.resize(m_wafer_hash_max, NOT_VALID_HASH);

    for (unsigned int i = 0; i < m_full_wafer_range.size(); ++i) {
        const Range& range = m_full_wafer_range[i];
        const Range::field& phi_field = range[m_PHI_MODULE_INDEX];
        const Range::field& eta_field = range[m_ETA_MODULE_INDEX];

        Range::const_identifier_factory first = range.factory_begin();
        Range::const_identifier_factory last  = range.factory_end();
        for (; first != last; ++first) {
            const ExpandedIdentifier& exp_id = (*first);
            ExpandedIdentifier::element_type previous_phi;
            ExpandedIdentifier::element_type next_phi;
            ExpandedIdentifier::element_type previous_eta;
            ExpandedIdentifier::element_type next_eta;
            bool pphi = phi_field.get_previous(exp_id[m_PHI_MODULE_INDEX], previous_phi);
            bool nphi = phi_field.get_next    (exp_id[m_PHI_MODULE_INDEX], next_phi);
            bool peta = eta_field.get_previous(exp_id[m_ETA_MODULE_INDEX], previous_eta);
            bool neta = eta_field.get_next    (exp_id[m_ETA_MODULE_INDEX], next_eta);

            IdContext      wcontext = wafer_context();
            
            // First get primary hash id
            IdentifierHash hash_id;
            Identifier id = wafer_id (exp_id[m_ENDCAP_INDEX],
                                      exp_id[m_LAYER_INDEX],
                                      exp_id[m_PHI_MODULE_INDEX],
                                      exp_id[m_ETA_MODULE_INDEX]);
            if (get_hash(id, hash_id, &wcontext)) {
                if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_neighbors - unable to get hash, exp/compact "
                                 << id.getString() << " " << show_to_string(id) << endmsg;
                else std::cout << " FATAL HGTD_ID::init_neighbors - unable to get hash, exp/compact "
                               << id.getString() << " " << show_to_string(id) << std::endl;
                return (1);
            }

            // index for the subsequent arrays
            unsigned short index = hash_id;
            assert (hash_id < m_prev_phi_wafer_vec.size());
            assert (hash_id < m_next_phi_wafer_vec.size());
            assert (hash_id < m_prev_eta_wafer_vec.size());
            assert (hash_id < m_next_eta_wafer_vec.size());
            
            if (pphi) {
                // Get previous phi hash id
                ExpandedIdentifier expId = exp_id;
                expId[m_PHI_MODULE_INDEX] = previous_phi;
                Identifier id = wafer_id (expId[m_ENDCAP_INDEX],
                                          expId[m_LAYER_INDEX],
                                          expId[m_PHI_MODULE_INDEX],
                                          expId[m_ETA_MODULE_INDEX]);
                if (get_hash(id, hash_id, &wcontext)) {
                    if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_neighbors - unable to get previous phi hash, exp/compact "
                                     << id.getString() << " " << show_to_string(id) << endmsg;
                    else std::cout << " FATAL HGTD_ID::init_neighbors - unable to get previous phi hash, exp/compact "
                                   << id.getString() << " " << show_to_string(id) << std::endl;
                    return (1);
                }
                m_prev_phi_wafer_vec[index] = hash_id;
            }
            
            if (nphi) {
                // Get next phi hash id
                ExpandedIdentifier expId = exp_id;
                expId[m_PHI_MODULE_INDEX] = next_phi;
                Identifier id = wafer_id (expId[m_ENDCAP_INDEX],
                                          expId[m_LAYER_INDEX],
                                          expId[m_PHI_MODULE_INDEX],
                                          expId[m_ETA_MODULE_INDEX]);
                if (get_hash(id, hash_id, &wcontext)) {
                    if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_neighbors - unable to get next phi hash, exp/compact " <<
                                     id.getString() << " " << show_to_string(id) << endmsg;
                    else std::cout << " FATAL HGTD_ID::init_neighbors - unable to get next phi hash, exp/compact " <<
                             id.getString() << " " << show_to_string(id) << std::endl;
                    return (1);
                }
                m_next_phi_wafer_vec[index] = hash_id;
            }
            
            if (peta) {
                // Get previous eta hash id
                ExpandedIdentifier expId = exp_id;
                expId[m_ETA_MODULE_INDEX] = previous_eta;
                Identifier id = wafer_id (expId[m_ENDCAP_INDEX],
                                          expId[m_LAYER_INDEX],
                                          expId[m_PHI_MODULE_INDEX],
                                          expId[m_ETA_MODULE_INDEX]);
                if (get_hash(id, hash_id, &wcontext)) {
                    if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_neighbors - unable to get previous eta hash, exp/compact "
                                     << id.getString() << " " << show_to_string(id) << endmsg;
                    else std::cout << " FATAL HGTD_ID::init_neighbors - unable to get previous eta hash, exp/compact "
                                   << id.getString() << " " << show_to_string(id) << std::endl;
                    return (1);
                }
                m_prev_eta_wafer_vec[index] = hash_id;
            }
            
            if (neta) {
                // Get next eta hash id
                ExpandedIdentifier expId = exp_id;
                expId[m_ETA_MODULE_INDEX] = next_eta;
                Identifier id = wafer_id (expId[m_ENDCAP_INDEX],
                                          expId[m_LAYER_INDEX],
                                          expId[m_PHI_MODULE_INDEX],
                                          expId[m_ETA_MODULE_INDEX]);
                if (get_hash(id, hash_id, &wcontext)) {
                    if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::init_neighbors - unable to get next eta hash, exp/compact "
                                     << id.getString() << " " << show_to_string(id) << endmsg;
                    else std::cout << " FATAL HGTD_ID::init_neighbors - unable to get next eta hash, exp/compact "
                                   << id.getString() << " " << show_to_string(id) << std::endl;
                    return (1);
                }
                m_next_eta_wafer_vec[index] = hash_id;
            }
            

//          std::cout << " HGTD_ID::init_neighbors "
//                    << " phi, previous, next " << id[m_PHI_MODULE_INDEX]
//                    << " " << pphi
//                    << " " << previous_phi
//                    << " " << nphi
//                    << " " << next_phi
//                    << " eta, previous, next " << id[m_ETA_MODULE_INDEX]
//                    << " " << peta
//                    << " " << previous_eta
//                    << " " << neta
//                    << " " << next_eta
//                    << " id " << (std::string)(*first)
//                    << std::endl;
        }
    }
    return (0);
}


void HGTD_ID::set_useNewIdentifierScheme(bool switchIntoNewIdentifier)
{
    m_useNewIdentifierScheme = switchIntoNewIdentifier;
}

bool HGTD_ID::get_useNewIdentifierScheme() const
{
    return m_useNewIdentifierScheme;
}

int
HGTD_ID::initLevelsFromDict(void)
{

    MsgStream log(m_msgSvc, "HGTD_ID");
    if(!m_dict) {
        if(m_msgSvc) log << MSG::FATAL << " HGTD_ID::initLevelsFromDict - dictionary NOT initialized " << endmsg;
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - dictionary NOT initialized " << std::endl;
        return (1);
    }

    // Find out which identifier field corresponds to each level. Use
    // names to find each field/level.

    m_INDET_INDEX               = 999;
    m_HGTD_INDEX                = 999;
    m_ENDCAP_INDEX              = 999;
    m_LAYER_INDEX               = 999;
    m_PHI_MODULE_INDEX          = 999;
    m_ETA_MODULE_INDEX          = 999;
    m_PHI_INDEX_INDEX           = 999;
    m_ETA_INDEX_INDEX           = 999;

    // Save index to a PIXEL region for unpacking
    ExpandedIdentifier id;
    id << indet_field_value() << hgtd_field_value();
    if (m_dict->find_region(id, m_hgtd_region_index)) {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find hgtd region index: id, reg "  
                         << (std::string)id << " " << m_hgtd_region_index
                         << endmsg;
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find hgtd region index: id, reg "  
                       << (std::string)id << " " << m_hgtd_region_index
                       << std::endl;
        return (1);
    }
    // Choose the name of the identification scheme based on the dictionary name, 
    // this information is gotten from HGTD_IDDetDescrCnv 
    
    if (m_useNewIdentifierScheme){
        m_endcap_ID               = "endcap";
        m_layer_ID                = "layer";
        m_moduleInLayer_Or_Row    = "moduleInLayer";
        m_moduleInRow             = "dummyVariable";
        m_padInModuleRow          = "padInModuleRow";
        m_padInModuleColumn       = "padInModuleColumn";
    }
    // Get levels
    IdDictField* field = m_dict->find_field("subdet");
    if (field) {
        m_INDET_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'subdet' field "        
                         << endmsg;
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'subdet' field "         
                       << std::endl;
        return (1);
    }

    field = m_dict->find_field("part");
    if (field) {
        m_HGTD_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'part' field " << endmsg;       
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'part' field " << std::endl;     
        return (1);
    }
 
    field = m_dict->find_field(m_endcap_ID);
    if (field) {
        m_ENDCAP_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'endcap' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'endcap' field "  << std::endl;   
            return (1);
    }
    field = m_dict->find_field(m_layer_ID);
    if (field) {
        m_LAYER_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'layer' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'layer' field "  << std::endl;   
        return (1);
    }
    field = m_dict->find_field(m_moduleInLayer_Or_Row);
    if (field) {
        m_PHI_MODULE_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'moduleInLayer' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'moduleInLayer' field "  << std::endl;   
        return (1);
    }
    field = m_dict->find_field(m_moduleInRow);
    if (field) {
        m_ETA_MODULE_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'hgtd_eta_module' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'hgtd_eta_module' field "  << std::endl;   
        return (1);
    }
    field = m_dict->find_field(m_padInModuleRow);
    if (field) {
        m_PHI_INDEX_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'hgtd_phi_index' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'hgtd_phi_index' field "  << std::endl;   
        return (1);
    }
    field = m_dict->find_field(m_padInModuleColumn);
    if (field) {
        m_ETA_INDEX_INDEX = field->m_index;
    }
    else {
        if(m_msgSvc) log << MSG::FATAL << "HGTD_ID::initLevelsFromDict - unable to find 'hgtd_eta_index' field "  << endmsg;     
        else std::cout << " FATAL HGTD_ID::initLevelsFromDict - unable to find 'hgtd_eta_index' field "  << std::endl;   
        return (1);
    }
    // Set the field implementations

    const IdDictRegion& region = *m_dict->m_regions[m_hgtd_region_index];

    m_indet_impl      = region.m_implementation[m_INDET_INDEX];
    m_hgtd_impl       = region.m_implementation[m_HGTD_INDEX];
    m_ec_impl         = region.m_implementation[m_ENDCAP_INDEX];
    m_layer_impl      = region.m_implementation[m_LAYER_INDEX];
    m_phi_mod_impl    = region.m_implementation[m_PHI_MODULE_INDEX];
    m_eta_mod_impl    = region.m_implementation[m_ETA_MODULE_INDEX];
    m_phi_index_impl  = region.m_implementation[m_PHI_INDEX_INDEX];
    m_eta_index_impl  = region.m_implementation[m_ETA_INDEX_INDEX];

    if(m_msgSvc) {
        log << MSG::DEBUG << "decode index and bit fields for each level: " << endmsg;
        log << MSG::DEBUG << "indet     "  << m_indet_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "hgtd      "  << m_hgtd_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "ec        "  << m_ec_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "layer     "  << m_layer_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "phi_mod   "  << m_phi_mod_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "eta_mod   "  << m_eta_mod_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "phi_index "  << m_phi_index_impl.show_to_string() << endmsg;
        log << MSG::DEBUG << "eta_index "  << m_eta_index_impl.show_to_string() << endmsg;
    }
    else{
        std::cout << " DEBUG decode index and bit fields for each level: " << std::endl;
        std::cout << " DEBUG indet     "  << m_indet_impl.show_to_string() << std::endl;
        std::cout << " DEBUG hgtd      "  << m_hgtd_impl.show_to_string() << std::endl;
        std::cout << " DEBUG ec        "  << m_ec_impl.show_to_string() << std::endl;
        std::cout << " DEBUG layer     "  << m_layer_impl.show_to_string() << std::endl;
        std::cout << " DEBUG phi_mod   "  << m_phi_mod_impl.show_to_string() << std::endl;
        std::cout << " DEBUG eta_mod   "  << m_eta_mod_impl.show_to_string() << std::endl;
        std::cout << " DEBUG phi_index "  << m_phi_index_impl.show_to_string() << std::endl;
        std::cout << " DEBUG eta_index "  << m_eta_index_impl.show_to_string() << std::endl;
    }

    std::cout << "indet "  << m_indet_impl.decode_index() << " " 
            << (std::string)m_indet_impl.ored_field() << " " 
            << std::hex << m_indet_impl.mask() << " " 
            << m_indet_impl.zeroing_mask() << " " 
            << std::dec << m_indet_impl.shift() 
            << " " << m_indet_impl.bits() << " " << m_indet_impl.bits_offset() << " ";
    m_indet_impl.ored_field().show();
    std::cout << "hgtd "  << m_hgtd_impl.decode_index() << " " 
            << (std::string)m_hgtd_impl.ored_field() << " " 
            << std::hex << m_hgtd_impl.mask() << " " 
            << m_hgtd_impl.zeroing_mask() << " " 
            << std::dec << m_hgtd_impl.shift() 
            << " " << m_hgtd_impl.bits() << " " << m_hgtd_impl.bits_offset() << " ";
    m_hgtd_impl.ored_field().show();
    std::cout << "ec "  << m_ec_impl.decode_index() << " " 
            << (std::string)m_ec_impl.ored_field() << " " 
            << std::hex << m_ec_impl.mask() << " " 
            << m_ec_impl.zeroing_mask() << " " 
            << std::dec << m_ec_impl.shift() 
            << " " << m_ec_impl.bits() << " " << m_ec_impl.bits_offset() << " ";
    m_ec_impl.ored_field().show();
    std::cout << "layer "  << m_layer_impl.decode_index() << " " 
            << (std::string)m_layer_impl.ored_field() << " " 
            << std::hex << m_layer_impl.mask() << " " 
            << m_layer_impl.zeroing_mask() << " " 
            << std::dec << m_layer_impl.shift() 
            << " " << m_layer_impl.bits() << " " << m_layer_impl.bits_offset() << " ";
    m_layer_impl.ored_field().show();
    std::cout << "phi_mod "  << m_phi_mod_impl.decode_index() << " " 
            << (std::string)m_phi_mod_impl.ored_field() << " " 
            << std::hex << m_phi_mod_impl.mask() << " " 
            << m_phi_mod_impl.zeroing_mask() << " " 
            << std::dec << m_phi_mod_impl.shift() 
            << " " << m_phi_mod_impl.bits() << " " << m_phi_mod_impl.bits_offset() << " ";
    m_phi_mod_impl.ored_field().show();
    std::cout << "eta_mod "  << m_eta_mod_impl.decode_index() << " " 
            << (std::string)m_eta_mod_impl.ored_field() << " " 
            << std::hex << m_eta_mod_impl.mask() << " " 
            << m_eta_mod_impl.zeroing_mask() << " " 
            << std::dec << m_eta_mod_impl.shift() 
            << " " << m_eta_mod_impl.bits() << " " << m_eta_mod_impl.bits_offset() << " ";
    m_eta_mod_impl.ored_field().show();
    std::cout << "phi_index "  << m_phi_index_impl.decode_index() << " " 
            << (std::string)m_phi_index_impl.ored_field() << " " 
            << std::hex << m_phi_index_impl.mask() << " " 
            << m_phi_index_impl.zeroing_mask() << " " 
            << std::dec << m_phi_index_impl.shift() 
            << " " << m_phi_index_impl.bits() << " " << m_phi_index_impl.bits_offset() << " ";
    m_phi_index_impl.ored_field().show();
    std::cout << "eta_index "  << m_eta_index_impl.decode_index() << " " 
            << (std::string)m_eta_index_impl.ored_field() << " " 
            << std::hex << m_eta_index_impl.mask() << " " 
            << m_eta_index_impl.zeroing_mask() << " " 
            << std::dec << m_eta_index_impl.shift() 
            << " " << m_eta_index_impl.bits() << " " << m_eta_index_impl.bits_offset() << " ";
    m_eta_index_impl.ored_field().show();


    std::cout << "HGTD_ID::initLevelsFromDict - found levels " << std::endl;
    std::cout << "subdet        "      << m_INDET_INDEX        << std::endl;
    std::cout << "part          "      << m_HGTD_INDEX         << std::endl;
    std::cout << "endcap        "      << m_ENDCAP_INDEX       << std::endl;
    std::cout << "layer         "      << m_LAYER_INDEX        << std::endl;
    std::cout << "phi_module    "      << m_PHI_MODULE_INDEX   << std::endl;
    std::cout << "eta_module    "      << m_ETA_MODULE_INDEX   << std::endl;
    std::cout << "phi_index     "      << m_PHI_INDEX_INDEX    << std::endl;
    std::cout << "eta_index     "      << m_ETA_INDEX_INDEX    << std::endl;

    return 0;
}

HGTD_ID::size_type
HGTD_ID::wafer_hash_max (void) const
{
    return m_wafer_hash_max;
}

HGTD_ID::size_type
HGTD_ID::pixel_hash_max (void) const
{
    return m_pixel_hash_max;
}

HGTD_ID::const_id_iterator      HGTD_ID::wafer_begin            (void) const
{
    return (m_wafer_vec.begin());
}

HGTD_ID::const_id_iterator      HGTD_ID::wafer_end              (void) const
{
    return (m_wafer_vec.end());
}

HGTD_ID::const_expanded_id_iterator     HGTD_ID::pixel_begin    (void) const
{
    return (m_full_pixel_range.factory_begin());
}

HGTD_ID::const_expanded_id_iterator     HGTD_ID::pixel_end      (void) const
{
    return (m_full_pixel_range.factory_end());
}

void
HGTD_ID::get_expanded_id        (const Identifier& id,
                                 ExpandedIdentifier& exp_id,
                                 const IdContext* context) const
{
    exp_id.clear();
    exp_id << indet_field_value()
           << hgtd_field_value()
           << endcap(id)
           << layer(id)
           << phi_module(id)
           << eta_module(id);
    if(!context || context->end_index() == m_ETA_INDEX_INDEX) {
        exp_id << phi_index(id)
               << eta_index(id);
    }
}

// From hash get Identifier
int
HGTD_ID::get_id (const IdentifierHash& hash_id,
                 Identifier& id,
                 const IdContext* context) const
{
    int result = 1;

    size_t begin = (context) ? context->begin_index(): 0;
    // cannot get hash if end is 0:
    size_t end   = (context) ? context->end_index()  : 0; 
    if (0 == begin) { 
        // No hashes yet for ids with prefixes
        if (m_ETA_MODULE_INDEX == end) {
            if (hash_id < (unsigned int)(m_wafer_vec.end() - m_wafer_vec.begin())) {
                id = m_wafer_vec[hash_id];
                result = 0;
            }
        }
        else if (m_ETA_INDEX_INDEX == end) {
            // Do not know how to calculate pixel id from hash yet!!
            std::cout << "Do not know how to calculate pixel id from hash yet!!" << std::endl;
        }
    }

    return (result);
}

int
HGTD_ID::get_hash       (const Identifier& id,
                         IdentifierHash& hash_id,
                         const IdContext* context) const
{

    // Get the hash code from either a vec (for wafers) or calculate
    // it (pixels). For the former, we convert to compact and call
    // get_hash again. For the latter, we calculate the hash from the
    // Identifier.

    int result = 1;

    hash_id = 0;
    size_t begin = (context) ? context->begin_index(): 0;
    size_t end   = (context) ? context->end_index()  : 0;
    if (0 == begin) {
        // No hashes yet for ids with prefixes
        if (m_ETA_MODULE_INDEX  == end) {
            hash_id = wafer_hash(id);
            if (hash_id.is_valid()) result = 0;
        }
        else if (context && context->end_index() == m_ETA_INDEX_INDEX) {
            // Must calculate for pixel hash
            ExpandedIdentifier new_id;
            get_expanded_id(id, new_id, context);
            hash_id =  m_full_pixel_range.cardinalityUpTo(new_id);
            result = 0;
        }
    }

    return (result);
}

void
HGTD_ID::test_wafer_packing     (void) const
{
    
    MsgStream log(m_msgSvc, "HGTD_ID");
    if (m_dict) {

        int nids = 0;
        IdContext context = wafer_context();
        const_id_iterator first = m_wafer_vec.begin();
        const_id_iterator last  = m_wafer_vec.end();
        for (; first != last; ++first, ++nids) {
            Identifier id = (*first);
            ExpandedIdentifier expId;
            get_expanded_id(id, expId, &context);
            Identifier new_id = wafer_id (expId[m_ENDCAP_INDEX],
                                          expId[m_LAYER_INDEX],
                                          expId[m_PHI_MODULE_INDEX],
                                          expId[m_ETA_MODULE_INDEX]);
            if (id != new_id) {
                if(m_msgSvc) log << MSG::ERROR << "HGTD_ID::test_wafer_packing: new and old compact id not equal. New/old/expanded ids "
                                 << show_to_string(new_id) << " " << show_to_string(id) << " "
                                 << (std::string)expId << endmsg;
                else std::cout << " ERROR HGTD_ID::test_wafer_packing: new and old compact id not equal. New/old/expanded ids " 
                               << show_to_string(new_id) << " " << show_to_string(id) << " "
                               << (std::string)expId << std::endl;
                continue;
            }
        }

        nids = 0;
        context = pixel_context();
        const_expanded_id_iterator      first_pixel = pixel_begin();
        const_expanded_id_iterator      last_pixel  = pixel_end();
        for (; first_pixel != last_pixel && nids < 1000; ++first_pixel, ++nids) {
            const ExpandedIdentifier& exp_id = *first_pixel;
            ExpandedIdentifier new_exp_id;

            Identifier id = wafer_id (exp_id[m_ENDCAP_INDEX],
                                      exp_id[m_LAYER_INDEX], 
                                      exp_id[m_PHI_MODULE_INDEX],
                                      exp_id[m_ETA_MODULE_INDEX]);

            get_expanded_id(id, new_exp_id, &context);
            if (exp_id[0] != new_exp_id[0] ||
                exp_id[1] != new_exp_id[1] ||
                exp_id[2] != new_exp_id[2] ||
                exp_id[3] != new_exp_id[3] ||
                exp_id[4] != new_exp_id[4] ||
                exp_id[5] != new_exp_id[5]) {
                if(m_msgSvc) log << MSG::ERROR << "HGTD_ID::test_wafer_packing: new and old expanded ids not equal. New/old/compact ids "
                                 << (std::string)new_exp_id 
                                 << " " << (std::string)exp_id 
                                 << " " <<  show_to_string(id) << endmsg;
                else std::cout << " ERROR HGTD_ID::test_wafer_packing: new and old expanded ids not equal. New/old/compact ids "
                               << (std::string)new_exp_id 
                               << " " << (std::string)exp_id 
                               << " " <<  show_to_string(id) << std::endl;
            }

            Identifier pid = pixel_id (exp_id[m_ENDCAP_INDEX],
                                       exp_id[m_LAYER_INDEX], 
                                       exp_id[m_PHI_MODULE_INDEX],
                                       exp_id[m_ETA_MODULE_INDEX],
                                       exp_id[m_PHI_INDEX_INDEX],
                                       exp_id[m_ETA_INDEX_INDEX]);
            Identifier pid1 = pixel_id (endcap(pid),
                                        layer(pid),
                                        phi_module(pid),
                                        eta_module(pid),
                                        phi_index(pid),
                                        eta_index(pid));
            if (pid != pid1) {
                if(m_msgSvc) log << MSG::ERROR << "HGTD_ID::test_wafer_packing: new and old pixel ids not equal. New/old ids "
                                 << " " << show_to_string(pid1) << " " 
                                 << show_to_string(pid) << endmsg;
                else std::cout << " ERROR HGTD_ID::test_wafer_packing: new and old pixel ids not equal. New/old ids "
                               << " " << show_to_string(pid1) << " " 
                               << show_to_string(pid) << std::endl;
            }
        }

        if(m_msgSvc) log << MSG::DEBUG << "HGTD_ID::test_wafer_packing: Successful tested " 
                         << nids << " ids. " 
                         << endmsg;
        else std::cout << " DEBUG HGTD_ID::test_wafer_packing: Successful tested " 
                       << nids << " ids. " 
                       << std::endl;
    }
    else {
        if(m_msgSvc) log << MSG::ERROR << "HGTD_ID::test_wafer_packing: Unable to test wafer is packing - no dictionary has been defined. " 
                         << endmsg;
        else std::cout << " ERROR HGTD_ID::test_wafer_packing: Unable to test wafer is packing - no dictionary has been defined. " 
                       << std::endl;
    }
}
