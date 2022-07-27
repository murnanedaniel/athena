/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

//////////////////////////////////////////////////////////
// HGTD_ID.h, (c) ATLAS Detector software
//////////////////////////////////////////////////////////

#ifndef HGTD_IDENTIFIER_HGTD_ID_H
#define HGTD_IDENTIFIER_HGTD_ID_H

#include "AtlasDetDescr/AtlasDetectorID.h"
#include "Identifier/Identifier.h"
#include "Identifier/Identifier32.h"
#include "Identifier/IdentifierHash.h"
#include "Identifier/Range.h"
#include "Identifier/IdHelper.h"
#include "IdDict/IdDictFieldImplementation.h"
#include "AthenaKernel/CLASS_DEF.h"

#include <string>
#include <assert.h>
#include <algorithm>

class IdDictDictionary;

/**
 **  @class HGTD_ID
 **  
 **  @brief This is an Identifier helper class for the HGTD
 **  subdetector. This class is a factory for creating compact
 **  Identifier objects and IdentifierHash or hash ids. And it also
 **  allows decoding of these ids.
 ** 
 **  Definition and the range of values for the levels of the
 **  identifier are:
 ** 
 ** @verbatim
 **    element           range    bits          meaning
 **    -------           -----    ----          -------
 **    TODO: fill in with latest identifier scheme
 ** @endverbatim
 ** 
 **/

class HGTD_ID : public AtlasDetectorID
{
public:

    /// @name public typedefs
    //@{
    typedef Identifier::size_type                       size_type; 
    typedef std::vector<Identifier>::const_iterator     const_id_iterator;
    typedef MultiRange::const_identifier_factory        const_expanded_id_iterator;
    //@}

    /// @name structors
    //@{
    HGTD_ID(void);
    ~HGTD_ID(void);
    //@}

    /// @name Creators for wafer ids and pixel ids
    //@{
    /// For a single crystal
    Identifier  wafer_id ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module ) const;

    /// For a single crystal from a pixel id - DO NOT USE wafer id as
    /// input
    Identifier  wafer_id ( const Identifier& pixel_id ) const;

    /// From hash
    Identifier  wafer_id ( IdentifierHash wafer_hash ) const;

    /// For an individual pixel
    Identifier  pixel_id ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module,
                           int phi_index,
                           int eta_index) const;
    Identifier  pixel_id ( const Identifier& id,
                           int phi_index,
                           int eta_index) const;
    //@}

    /// @name Hash table maximum sizes
    //@{
    size_type   wafer_hash_max          (void) const;
    size_type   pixel_hash_max          (void) const;
    //@}

    /// @name Access to all ids
    //@{
    /// Iterators over full set of ids. Wafer iterator is sorted
    const_id_iterator   wafer_begin     (void) const;
    const_id_iterator   wafer_end       (void) const;
    /// For pixel ids, only expanded id iterators are available. Use
    /// following "pixel_id" method to obtain a compact identifier
    const_expanded_id_iterator  pixel_begin     (void) const;  
    const_expanded_id_iterator  pixel_end       (void) const;
    //@}

    /// wafer hash from id 
    IdentifierHash      wafer_hash      (Identifier wafer_id) const;

    /// Values of different levels (failure returns 0)
    int   endcap          (const Identifier& id) const;
    int   layer           (const Identifier& id) const;
    int   phi_module      (const Identifier& id) const;
    int   eta_module      (const Identifier& id) const;
    int   phi_index       (const Identifier& id) const;
    int   eta_index       (const Identifier& id) const;

    /// Max/Min values for each field (error returns -999)
    int   layer_max       (const Identifier& id) const;
    int   phi_module_max  (const Identifier& id) const;
    int   eta_module_max  (const Identifier& id) const;
    int   eta_module_min  (const Identifier& id) const;
    int   phi_index_max   (const Identifier& id) const;
    int   eta_index_max   (const Identifier& id) const;

    /// @name module eta/phi navigation
    //@{
    /// Previous wafer hash in phi (return == 0 for neighbor found)
    int         get_prev_in_phi (const IdentifierHash& id, IdentifierHash& prev) const;
    /// Next wafer hash in phi (return == 0 for neighbor found)
    int         get_next_in_phi (const IdentifierHash& id, IdentifierHash& next) const;
    /// Previous wafer hash in eta (return == 0 for neighbor found)
    int         get_prev_in_eta (const IdentifierHash& id, IdentifierHash& prev) const;
    /// Next wafer hash in eta (return == 0 for neighbor found)
    int         get_next_in_eta (const IdentifierHash& id, IdentifierHash& next) const;

    /// To check for when phi wrap around may be needed
    bool        is_phi_module_max(const Identifier& id) const;
    //@}

    /// @name contexts to distinguish wafer id from pixel id
    //@{
    IdContext   wafer_context           (void) const;
    IdContext   pixel_context           (void) const;
    //@}

    /// @name methods from abstract interface - slower than opt version
    //@{
    /// Create compact id from hash id (return == 0 for OK)
    virtual int         get_id          (const IdentifierHash& hash_id,
                                         Identifier& id,
                                         const IdContext* context = 0) const;
    
    /// Create hash id from compact id (return == 0 for OK)
    virtual int         get_hash        (const Identifier& id, 
                                         IdentifierHash& hash_id,
                                         const IdContext* context = 0) const;
    //@}

    // TODO: see if these commented out methods ported from PixelID class are needed

    // /// Create a compact id from a value (e.g., from a persistent object).
    // /// This repacks fields in case it's a special pixel channel id.
    // Identifier          pixel_id        (Identifier::value_type val) const;

    // /// Test if this is a valid shortened pixel channel id.
    // bool is_shortened_pixel_id(Identifier32::value_type val) const;
    // bool is_shortened_pixel_id(const Identifier32& id) const;
    // bool is_shortened_pixel_id(const Identifier& id) const;

    // /// Create a compact pixel id from a (fixed format) legacy
    // /// pixel channel id.  If compiled in 32-bit mode, this method
    // /// does nothing.
    // Identifier pixel_id_from_shortened(Identifier32::value_type val) const;

    /// Return the lowest bit position used in the channel id.
    int                 base_bit        (void) const;

    /// Calculate a channel offset between the two identifiers.
    Identifier::diff_type calc_offset(const Identifier& base,
                                      const Identifier& target) const;

    /// Create an identifier with a given base and channel offset
    Identifier pixel_id_offset(const Identifier& base,
                               Identifier::diff_type offset) const;

    /// @name interaction with id dictionary
    //@{
    /// Create pixel Identifier from expanded id, which is returned by the
    /// id_iterators
    Identifier          pixel_id        (const ExpandedIdentifier& pixel_id) const;

    /// Create expanded id from compact id (return == 0 for OK)
    void                get_expanded_id (const Identifier& id,
                                         ExpandedIdentifier& exp_id,
                                         const IdContext* context = 0) const;

    /// Initialization from the identifier dictionary
    virtual int         initialize_from_dictionary(const IdDictMgr& dict_mgr);

    /// Tests of packing
    void        test_wafer_packing      (void) const;

    // switch between identification schemes
    bool m_useNewIdentifierScheme = false;
    void set_useNewIdentifierScheme(bool switchIntoNewIdentifier);
    bool get_useNewIdentifierScheme() const;
    std::string m_endcap_ID               = "hgtd_endcap";
    std::string m_layer_ID                = "hgtd_layer";
    std::string m_moduleInLayer_Or_Row    = "hgtd_phi_module";
    std::string m_moduleInRow             = "hgtd_eta_module";
    std::string m_padInModuleRow          = "hgtd_phi_index";
    std::string m_padInModuleColumn       = "hgtd_eta_index";
    //@}

private:
    enum {NOT_VALID_HASH        = 64000,
          MAX_BIT               = Identifier::MAX_BIT,
          BITS32                = Identifier::ALL_BITS };

    typedef std::vector<Identifier>     id_vec;
    typedef id_vec::const_iterator      id_vec_it;
    typedef std::vector<unsigned short> hash_vec;
    typedef hash_vec::const_iterator    hash_vec_it;

    void wafer_id_checks ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module ) const;
    void pixel_id_checks ( int endcap,
                           int layer,
                           int phi_module,
                           int eta_module,
                           int phi_index,
                           int eta_index) const;

    int         initLevelsFromDict(void);

    int         init_hashes(void);

    int         init_neighbors(void);

    size_type                   m_hgtd_region_index;
    size_type                   m_INDET_INDEX;
    size_type                   m_HGTD_INDEX;
    size_type                   m_ENDCAP_INDEX;
    size_type                   m_LAYER_INDEX;
    size_type                   m_PHI_MODULE_INDEX;
    size_type                   m_ETA_MODULE_INDEX;
    size_type                   m_PHI_INDEX_INDEX;
    size_type                   m_ETA_INDEX_INDEX;

    const IdDictDictionary*     m_dict;
    MultiRange                  m_full_wafer_range;
    MultiRange                  m_full_pixel_range;
    size_type                   m_wafer_hash_max;
    size_type                   m_pixel_hash_max;

    id_vec                      m_wafer_vec;

    hash_vec                    m_prev_phi_wafer_vec;
    hash_vec                    m_next_phi_wafer_vec;
    hash_vec                    m_prev_eta_wafer_vec;
    hash_vec                    m_next_eta_wafer_vec;

    IdDictFieldImplementation   m_indet_impl;
    IdDictFieldImplementation   m_hgtd_impl;
    IdDictFieldImplementation   m_ec_impl;
    IdDictFieldImplementation   m_layer_impl;
    IdDictFieldImplementation   m_phi_mod_impl;
    IdDictFieldImplementation   m_eta_mod_impl;
    IdDictFieldImplementation   m_phi_index_impl;
    IdDictFieldImplementation   m_eta_index_impl;

};

//using the macros below we can assign an identifier (and a version)
//This is required and checked at compile time when you try to record/retrieve
CLASS_DEF(HGTD_ID, 79264207, 1)

///////////////////////////////////////////////////////////////////
// Inline methods:
///////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
inline Identifier
HGTD_ID::wafer_id ( int endcap,
                    int layer,
                    int phi_module,
                    int eta_module ) const
{

    // Build identifier
    Identifier result((Identifier::value_type)0);

    // Pack fields independently
    m_indet_impl.pack    (indet_field_value(), result);
    m_hgtd_impl.pack     (hgtd_field_value(),  result);
    m_ec_impl.pack       (endcap,              result);
    m_layer_impl.pack    (layer,               result);
    m_phi_mod_impl.pack  (phi_module,          result);
    m_eta_mod_impl.pack  (eta_module,          result);

    // Do checks
    if(m_do_checks) {
        wafer_id_checks ( endcap, layer, phi_module, eta_module );
    }

    return result;
}

//----------------------------------------------------------------------------
inline Identifier  
HGTD_ID::wafer_id ( const Identifier& pixel_id ) const
{
    Identifier result(pixel_id);
    m_phi_index_impl.reset      (result);
    m_eta_index_impl.reset      (result);
    return (result);
}

//----------------------------------------------------------------------------
inline Identifier
HGTD_ID::wafer_id ( IdentifierHash wafer_hash ) const
{
    return (m_wafer_vec[wafer_hash]);
}

//----------------------------------------------------------------------------
inline Identifier
HGTD_ID::pixel_id ( int endcap,
                    int layer,
                    int phi_module,
                    int eta_module,
                    int phi_index,
                    int eta_index) const
{

    // Build identifier
    Identifier result((Identifier::value_type)0);
    m_indet_impl.pack          (indet_field_value(), result);
    m_hgtd_impl.pack           (hgtd_field_value(),  result);
    m_ec_impl.pack             (endcap,              result);
    m_layer_impl.pack          (layer,               result);
    m_phi_mod_impl.pack        (phi_module,          result);
    m_eta_mod_impl.pack        (eta_module,          result);
    m_phi_index_impl.pack      (phi_index,           result);
    m_eta_index_impl.pack      (eta_index,           result);

    if(m_do_checks) {

        pixel_id_checks ( endcap,
                          layer,
                          phi_module,
                          eta_module,
                          phi_index,
                          eta_index);
    }
    
    return result;
}

/// Create pixel Identifier from expanded id, which is returned by the
/// id_iterators
//----------------------------------------------------------------------------
inline Identifier
HGTD_ID::pixel_id       (const ExpandedIdentifier& id) const
{
    Identifier result;
    result = pixel_id(id[m_ENDCAP_INDEX],
                      id[m_LAYER_INDEX],
                      id[m_PHI_MODULE_INDEX],
                      id[m_ETA_MODULE_INDEX],
                      id[m_PHI_INDEX_INDEX],
                      id[m_ETA_INDEX_INDEX]);

    if(m_do_checks) {
        pixel_id_checks (id[m_ENDCAP_INDEX],
                         id[m_LAYER_INDEX],
                         id[m_PHI_MODULE_INDEX],
                         id[m_ETA_MODULE_INDEX],
                         id[m_PHI_INDEX_INDEX],
                         id[m_ETA_INDEX_INDEX]);
    }

    return (result);
}

//----------------------------------------------------------------------------
inline Identifier  
HGTD_ID::pixel_id ( const Identifier& wafer_id,
                    int phi_index,
                    int eta_index) const
{
    Identifier result(wafer_id);
    m_phi_index_impl.reset     (result);
    m_eta_index_impl.reset     (result);
    m_phi_index_impl.pack      (phi_index, result);
    m_eta_index_impl.pack      (eta_index, result);
    return (result);
}

//----------------------------------------------------------------------------
inline IdentifierHash      HGTD_ID::wafer_hash      (Identifier wafer_id) const
{
    id_vec_it it = std::lower_bound(m_wafer_vec.begin(),
                                    m_wafer_vec.end(),
                                    wafer_id);
    // Require that wafer_id matches the one in vector
    if (it != m_wafer_vec.end() && wafer_id == (*it)) {
        return (it - m_wafer_vec.begin());
    }
    IdentifierHash result;
    return (result); // return hash in invalid state
}

//----------------------------------------------------------------------------
inline Identifier::diff_type
HGTD_ID::calc_offset(const Identifier& base, const Identifier& target) const
{
  Identifier::diff_type tval = static_cast<Identifier::diff_type>(target.get_compact() >> base_bit());
  Identifier::diff_type bval = static_cast<Identifier::diff_type>(base.get_compact() >> base_bit());
  return (tval - bval);
}

//----------------------------------------------------------------------------
inline Identifier
HGTD_ID::pixel_id_offset(const Identifier& base,
                         Identifier::diff_type offset) const
{
  Identifier::value_type bval = base.get_compact() >> base_bit();
  return Identifier((bval + offset) << base_bit());
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::base_bit ( void ) const
{
    // TODO: determine if anything needs to change here (implementation ported from PixelID)
    // TODO: e.g. m_eta_index_impl is still lowest field base, but where does the 32 come from? 32 bits?
    // TODO: also understand what's happening
    int base = static_cast<int>(m_eta_index_impl.shift()); // lowest field base
    return (base > 32) ? 32 : base;
    // max base is 32 so we can still read old strip id's and differences
    // from non-SLHC releases.
}

//----------------------------------------------------------------------------
inline IdContext
HGTD_ID::wafer_context          (void) const
{
    ExpandedIdentifier id;
    return (IdContext(id, 0, m_ETA_MODULE_INDEX));
}

//----------------------------------------------------------------------------
inline IdContext
HGTD_ID::pixel_context  (void) const
{
    // For pixel only, the prefix is the first two levels
    ExpandedIdentifier id;
    id << indet_field_value() << hgtd_field_value();
    return (IdContext(id, m_ENDCAP_INDEX, m_ETA_INDEX_INDEX));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::endcap       (const Identifier& id) const
{
    return (m_ec_impl.unpack(id));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::layer       (const Identifier& id) const
{
    return (m_layer_impl.unpack(id));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::phi_module       (const Identifier& id) const
{
    return (m_phi_mod_impl.unpack(id));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::eta_module       (const Identifier& id) const
{
    return (m_eta_mod_impl.unpack(id));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::phi_index       (const Identifier& id) const
{
    return (m_phi_index_impl.unpack(id));
}

//----------------------------------------------------------------------------
inline int
HGTD_ID::eta_index       (const Identifier& id) const
{
    return (m_eta_index_impl.unpack(id));
}

#endif // HGTD_IDENTIFIER_HGTD_ID_H
