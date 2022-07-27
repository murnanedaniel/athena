/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRT_RAWDATABYTESTREAM_TRT_RODDECODER_H
#define TRT_RAWDATABYTESTREAM_TRT_RODDECODER_H


/*
 * Base class
 */
#include "AthenaBaseComps/AthAlgTool.h"  

/*
 * Interface class for this Tool
 */
#include "TRT_RawDataByteStreamCnv/ITRT_RodDecoder.h"

/*
 * TRT Tools we use
 */
#include "TRT_Cabling/ITRT_CablingSvc.h"

/*
 * Framework Headers
 *   Service & Tool Handles
 */
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ContextSpecificPtr.h"
#include "GaudiKernel/ThreadLocalContext.h"
#include "GaudiKernel/ICondSvc.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "CxxUtils/CachedUniquePtr.h"

#include "CoralBase/Attribute.h"

/*
 * Identifier
 */
#include "InDetIdentifier/TRT_ID.h"

/*
 * For cache
 */
#include "AthenaKernel/SlotSpecificObj.h"
#include "CxxUtils/checker_macros.h"

/*
 * STL
 */
#include <atomic>
#include <map>
#include <vector>

// the tool to decode a ROB frament

class TRT_RodDecoder : public extends<AthAlgTool, ITRT_RodDecoder>
{
public: 
  //! constructor
  TRT_RodDecoder(const std::string& type, const std::string& name,
                 const IInterface* parent ) ;
  //! destructor 
  virtual ~TRT_RodDecoder(); 
  //! initialize
  virtual StatusCode initialize() override;
  //! finalize
  virtual StatusCode finalize() override;

  //! the method to fill the IDC
  virtual StatusCode fillCollection ( const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment* robFrag,
				      TRT_RDO_Container* rdoIdc,
				      TRT_BSErrContainer* bserr,
				      const std::vector<IdentifierHash>* vecHash = 0) const override;


 private:

   ServiceHandle<ITRT_CablingSvc>   m_CablingSvc;

   /*
    * Do we look for Front-End Errors at all?
    */
   bool m_recordBSErrors;

   /*
    * Do we look for these types of Front-End Errors?
    */
   bool m_lookAtSidErrors;
   bool m_lookAtErrorErrors;
   bool m_lookAtL1idErrors;
   bool m_lookAtBcidErrors;
   bool m_lookAtMissingErrors;

   bool m_loadCompressTableFile;
   bool m_loadCompressTableDB;
   std::vector<int> m_LoadCompressTableVersions;
   const int m_maxCompressionVersion;
   int m_forceRodVersion;

   const TRT_ID*               m_trt_id;   
   IdContext                   m_straw_layer_context;

   BooleanProperty m_TB04_RawData;   // true to create TRT_TB04_RawData RDOs
   BooleanProperty m_LoLumRawData;   // true to create TRT_LoLumRawData RDOs

   bool m_eventTypeIsSim;            // look at first event to decide if RODBlockVersion
                                    // is incorrect.

   uint32_t m_escape_marker;           // Straw word that means escaped literal

#define CTABLE_FC_LENGTH 33
#define CTABLE_LI_LENGTH 33
   typedef struct {
     // The TableVersion is a compression scheme.
     // There are presently 6 TableVersions for the different run periods of ATLAS.
     // Any IoV of /TRT/onl/ROD/Compress has one (and only one) of these as an attribute.
     // The code takes advantage of this and caches locally any TableVersion that is encountered
     int m_TableVersion;
     int m_firstcode[CTABLE_FC_LENGTH];
     int m_lengths_integral[CTABLE_FC_LENGTH];    // ..[i] = Sum(numl[0,i-1])
     std::unique_ptr<unsigned int[]> m_syms;              // Array of symbols (straw data words)
     int m_Nsymbols;
   } t_CompressTable;

   std::vector<CxxUtils::CachedUniquePtr<t_CompressTable> > m_CompressionTables;

   mutable std::atomic<uint32_t> m_Nrdos;              // Number of RDOs created

   mutable std::atomic<int> m_err_count_fillCollection{0};
   mutable std::atomic<int> m_err_count_int_fillMinimalCompress{0};
   mutable std::atomic<int> m_err_count_int_fillFullCompress{0};

   // This replaces the IOVCALLBACK
   SG::ReadCondHandleKey<CondAttrListCollection> m_CompressKey{this,"keyName","/TRT/Onl/ROD/Compress","in-key"};
   StatusCode update() const;

   //! private methods
private:

   StatusCode int_fillExpanded ( const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment* robFrag,
				TRT_RDO_Container* rodIdc,
				const std::vector<IdentifierHash>* vecHash = 0) const;

   StatusCode int_fillMinimalCompress ( const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment* robFrag,
				       TRT_RDO_Container* rdoIdo,
				       const std::vector<IdentifierHash>* vecHash = 0) const;

   StatusCode int_fillFullCompress ( const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment* robFrag,
				     TRT_RDO_Container* rdoIdo,
				     const t_CompressTable& Ctable,
				     const std::vector<IdentifierHash>* vecHash = 0) const;

   StatusCode ReadCompressTableFile( const std::string& TableFilename );
   StatusCode ReadCompressTableDB( std::string Tag );

   mutable SG::SlotSpecificObj<std::atomic<EventContext::ContextEvt_t> > m_lastPrint ATLAS_THREAD_SAFE;

   mutable std::atomic<unsigned int> m_skip{};
   mutable std::atomic<unsigned int> m_accept{};
};

#endif
