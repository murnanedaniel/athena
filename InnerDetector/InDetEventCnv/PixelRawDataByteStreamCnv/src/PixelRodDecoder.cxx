/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "PixelRodDecoder.h"
#include "CxxUtils/AthUnlikelyMacros.h"
#include "InDetIdentifier/PixelID.h"
#include "PixelReadoutGeometry/PixelDetectorManager.h"
#include "ExtractCondensedIBLhits.h"
#include "PixelByteStreamModuleMask.h"
#include "ByteStreamData/RawEvent.h"
#include "eformat/SourceIdentifier.h"
#include "xAODEventInfo/EventInfo.h"

#include <fstream>
#include <iostream>
#include <unordered_set>
#include <string>

//#define PIXEL_DEBUG ;
namespace {
  constexpr unsigned long long operator"" _BIT(unsigned long long bitPosition){
    return 1ull << bitPosition;
  }
}

#define generalwarning(x)						\
  if (this->m_numGenWarnings < this->m_maxNumGenWarnings)		\
    {ATH_MSG_WARNING(x); ++this->m_numGenWarnings;}			\
  else if (this->m_numGenWarnings == this->m_maxNumGenWarnings)		\
    {ATH_MSG_INFO("PixelRodDecoder: suppressing further general warnings"); ++this->m_numGenWarnings;} \
  else {++this->m_numGenWarnings; /* No warning */}

#define lvl1id_bcid_warning(x)						\
  if (this->m_numBCIDWarnings < this->m_maxNumBCIDWarnings)		\
    {ATH_MSG_WARNING(x); ++this->m_numBCIDWarnings;}			\
  else if (this->m_numBCIDWarnings == this->m_maxNumBCIDWarnings)	\
    {ATH_MSG_INFO("PixelRodDecoder: suppressing further BCID/LVL1ID warnings"); ++this->m_numBCIDWarnings;} \
  else {++this->m_numBCIDWarnings; /* No warning */}


using OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment;

inline bool isIBL( uint32_t robId ) { return ((robId>>16) & 0xFF)==0x14; }
inline bool isDBM( uint32_t robId ) { return ((robId>>16) & 0xFF)==0x15; }

//--------------------------------------------------------------------------- constructor
PixelRodDecoder::PixelRodDecoder
( const std::string& type, const std::string& name,const IInterface* parent )
  :  AthAlgTool(type,name,parent),
     m_is_ibl_present(false){
  declareInterface< IPixelRodDecoder  >( this );
}

//--------------------------------------------------------------------------- destructor
//destructor
PixelRodDecoder::~PixelRodDecoder() {

}

//--------------------------------------------------------------------------- initialize


StatusCode PixelRodDecoder::initialize() {
  ATH_MSG_INFO("PixelRodDecoder::initialize()");
  ATH_CHECK(m_pixelReadout.retrieve());
  ATH_CHECK(detStore()->retrieve(m_pixel_id, "PixelID"));
  const InDetDD::PixelDetectorManager* pixelManager;
  ATH_CHECK(detStore()->retrieve(pixelManager, "Pixel"));
  // check if the ibl layer is present or not (this is necessary for backward compatibility with non-IBL geometries)
  m_is_ibl_present = false;
  const InDetDD::SiNumerology& pixSiNum = pixelManager->numerology();
  m_is_ibl_present = (pixSiNum.numLayers() == 4);
  ATH_MSG_DEBUG( "m_is_ibl_present = " << m_is_ibl_present );
  // Retrieve Pixel Errors Service
  ATH_CHECK(m_condCablingKey.initialize());
  ATH_CHECK(m_condHitDiscCnfgKey.initialize(m_is_ibl_present));
  return StatusCode::SUCCESS;
}


//--------------------------------------------------------------------------- finalize
StatusCode PixelRodDecoder::finalize() {

#ifdef PIXEL_DEBUG
  ATH_MSG_VERBOSE( "in PixelRodDecoder::finalize" );
  ATH_MSG_DEBUG(m_masked_errors << " times BCID and LVL1ID error masked" );
#endif

  ATH_MSG_INFO("Total number of warnings (output limit)");
  ATH_MSG_INFO("General (corruption / illegal value)\t" << m_numGenWarnings << " (" << m_maxNumGenWarnings << ")");
  ATH_MSG_INFO("Unexpected BCID / LVL1ID            \t" << m_numBCIDWarnings << " (" << m_maxNumBCIDWarnings << ")");
  if (m_checkDuplicatedPixel) ATH_MSG_INFO("Duplicated pixels                     \t" << m_numDuplicatedPixels);
  ATH_MSG_INFO(" --- Summary from PixelRodDecoder --- ");
  ATH_MSG_INFO(" - Number of invalid Identifiers Errors: " << m_numInvalidIdentifiers);
  ATH_MSG_INFO(" - Number of Preamble Errors: " << m_numPreambleErrors);
  ATH_MSG_INFO(" - Number of TimeOut Errors: " << m_numTimeOutErrors);
  ATH_MSG_INFO(" - Number of LVL1ID Errors: " << m_numLVL1IDErrors);
  ATH_MSG_INFO(" - Number of BCID Errors: " << m_numBCIDErrors);
  ATH_MSG_INFO(" - Number of Flagged Errors: " << m_numFlaggedErrors);
  ATH_MSG_INFO(" - Number of Trailer Errors: " << m_numTrailerErrors);
  ATH_MSG_INFO(" - Number of Disabled FE Errors: " << m_numDisabledFEErrors);
  ATH_MSG_INFO(" - Number of ROD Errors: " << m_numRODErrors);
  ATH_MSG_INFO(" - Number of links masked by PPC: " << m_numLinkMaskedByPPC);
  ATH_MSG_INFO(" - Number of header/trailer limit errors: " << m_numLimitError);
  ATH_MSG_INFO(" - Number of Unknown word Errors: " << m_numDecodingErrors);
  ATH_MSG_INFO(" --------------------------------------------- ");
  return StatusCode::SUCCESS;
}

void PixelRodDecoder::propagateROBErrorsToModules(const PixelCablingCondData *pixCabling,
                                                  uint32_t robId,
                                                  std::array<uint64_t, PixelRodDecoder::ERROR_CONTAINER_MAX> &bsErrWord,
                                                  IDCInDetBSErrContainer& decodingErrors,
                                                  PixelByteStreamErrors::PixelErrorsEnum error_code,
                                                  const char *error_description) const {
   assert( pixCabling);
   const std::deque<Identifier> offlineIdList = pixCabling->find_entry_offlineList(robId);
   for (const Identifier& id: offlineIdList) {
      IdentifierHash idHash = m_pixel_id->wafer_hash(id);
      PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(idHash)],error_code);
   }
   ATH_MSG_DEBUG("ROB status word for robid 0x"<< std::hex << robId << std::dec <<" indicates " << error_description << ".");
   assert( bsErrWord.size() <= decodingErrors.maxSize() );
   for (size_t i=0; i<static_cast<size_t>(bsErrWord.size()); i++) {
      if (bsErrWord[i]>0) {
         decodingErrors.setOrDrop(i,bsErrWord[i]);
      }
   }
}

//---------------------------------------------------------------------------------------------------- fillCixollection
StatusCode PixelRodDecoder::fillCollection( const ROBFragment *robFrag, IPixelRDO_Container* rdoIdc, IDCInDetBSErrContainer& decodingErrors, std::vector<IdentifierHash>* vecHash) const {
#ifdef PIXEL_DEBUG
  ATH_MSG_VERBOSE( "-------------------------------------------------------------------------------------------------------------");
  ATH_MSG_DEBUG("Entering PixelRodDecoder");
#endif

  const Identifier invalidPixelId = Identifier(); // used by Cabling for an invalid entry
  Identifier pixelId;
  uint64_t onlineId(0);
  bool isIBLModule = false;
  bool isDBMModule = false;
  uint32_t rodId = robFrag->rod_source_id(); // get source ID => method of ROBFragment, Returns the source identifier of the ROD fragment
  uint32_t robId = robFrag->rob_source_id(); // get source ID => returns the Identifier of the ROB fragment. More correct to use w.r.t. rodId.
  uint32_t robBCID = robFrag->rod_bc_id();
  if ( isIBL( robId )  && robId != rodId ) { // isIBL(robId)
    generalwarning("Discrepancy in IBL SourceId: ROBID 0x" << std::hex << robId << " unequal to RODID 0x" << rodId);
  }
  SG::ReadCondHandle<PixelCablingCondData> pixCabling(m_condCablingKey);
  std::unique_ptr<SG::ReadCondHandle<PixelHitDiscCnfgData> > pixHitDiscCnfg;
  unsigned int errorRecoverable = 0;
  uint32_t serviceCodeCounter = 0;  // frequency of the serviceCode (with the exceptions of serviceCode = 14,15 or 16)
  uint32_t serviceCode = 0;         // is a code. the corresponding meaning is listed in the table in the FE-I4 manual, pag. 105

  // @TODO find better solution for the error counter to avoid complex index computations and hard coded maximum size.
  // The index array is defined in PixelRawDataProviderTool::SizeOfIDCInDetBSErrContainer()
  std::array<uint64_t, PixelRodDecoder::ERROR_CONTAINER_MAX> bsErrWord{};
  std::fill(bsErrWord.begin(),bsErrWord.end(),0);
  // Check ROD status
  if (robFrag->nstatus()!=0) {
    const uint32_t* rob_status;
    robFrag->status(rob_status);
    // check the ROD status for truncation
    if (ATH_UNLIKELY((*rob_status)!=0)) {
      ATH_MSG_DEBUG( "ROB status word for robid 0x"<< std::hex << robId << " is non-zero 0x" << (*rob_status) << std::dec);
      /*
        Definition of the status words in a ROB fragment header is found in
           https://twiki.cern.ch/twiki/bin/view/Atlas/ROBINFragmentErrors#Definition_of_the_first_status_e
      */
      if ((*rob_status) & (0x1u << 27)) {
         propagateROBErrorsToModules(pixCabling.cptr(),robId,bsErrWord,decodingErrors,PixelByteStreamErrors::TruncatedROB, "data truncation");
         return StatusCode::RECOVERABLE;
      }
      if ((*rob_status) & (0x1u << 31)) {
         propagateROBErrorsToModules(pixCabling.cptr(),robId,bsErrWord,decodingErrors,PixelByteStreamErrors::MaskedROB, "resource was masked off");
         return StatusCode::RECOVERABLE;
      }
      // in case the ROB fragment has a seemingly invalid size check the fragment and reject it if the check is not passed.
      // Note: there are usable ROB fragments (i.e. ROB fragments which contribute pixel hits to tracks) which do not pass the check, so
      // rejecting all fragments which do not pass the test would reject also seemingly "good" data.
      if (robFrag->rod_ndata() > robFrag->payload_size_word() && !robFrag->check_rod_noex(robFrag->rod_version() >> 16)) {
         propagateROBErrorsToModules(pixCabling.cptr(),robId,bsErrWord,decodingErrors,PixelByteStreamErrors::TruncatedROB,
                                     " invalid ROD fragment, invalid payload size");
         return StatusCode::RECOVERABLE;
      }
    }
  }
  StatusCode sc = StatusCode::SUCCESS;

  // Fill the Collections for this ROB.
  unsigned int mBCID           = 0x3FFF;   // create module BCID variable
  unsigned int mLVL1ID         = 0x3FFF;   // create module LVL1ID variable
  unsigned int mLink           = 0x3FFF;   // create module Link variable
  int mLVL1A                   = 0x3FFF;   // create LVL1A variable
  int offsetBCID_ROB_FE        = 0x3FFF;   // offset between BCID in ROD header and FE header
  // take the Header "nnnnn" bits as a reference for the other words (hit, feFlag, trailer words)
  uint32_t linkNum_IBLheader = 0; // create variable containing 5 bits, corresponding to the "nnnnn" bits in the header IBL/DBM word
  uint32_t fe_IBLheader = 0; // create fe variable, corresponding to the IBL FE number w.r.t. the sLink (bits 24-26 "nnn" in Header)
  uint32_t linkNum_IBLword = 0; // create variable containing 5 bits, corresponding to the "nnnnn" bits in the hit/feFlag/trailer IBL/DBM word // At least temporarily not used, because IBL data format is not clear (Franconi, 17.06.2014)
  uint32_t sLinkSourceId = 0; // create sLink variable, corresponding to the number of the sLink of IBL/DBM, runs from 1 to 4. Retrieved from the robId
  uint32_t sLinkHeader = 0; // create sLink variable, retrieved from the header words
  // Various bookkeeping variables
  bool link_start = false;   // link (module) header found flag - reset by link (module) trailer
  IdentifierHash offlineIdHash; // offlineIdHash
  IdentifierHash previous_offlineIdHash;   // previous module offline identifier
  int mBCID_offset = -1;   // BCID offset to calculate L1A
  unsigned int mBCID_max_pix = 0x100;   // cycle of BCID values
  unsigned int mBCID_max_IBL = 0x400;
  bool receivedCondensedWords = false;    // to cross-check the condensed mode bit ('c')
  int countHitCondensedWords = 0; // counter to the "condensed hits" of the IBL
  uint32_t condensedWords[nCondensedWords] = {}; // this vector saves the 4 condensed hit words that are found
  bool are_4condensed_words = false;
  unsigned int prevBCID           = 0x3FFF;   // initialize to out-of-range value
  unsigned int prevLVL1ID         = 0x3FFF;
  uint32_t prevLinkNum            = 99;
  int prevOffsetBCID_ROB_FE       = 0x3FFF;
  unsigned int prevStartingBCID   = 0x3FFF;
  unsigned int nFragmentsPerFE[8] = {0};    // number of header-trailer pairs per IBL FE
  IdentifierHash skipHash = 0xffffffff, lastHash = 0xffffffff; // used for decoding to not have to search the vector all the time
  PixelRawCollection* coll = nullptr;
  // get the data from the word
  OFFLINE_FRAGMENTS_NAMESPACE::PointerType vint;
  robFrag->rod_data(vint);
  int hitDiscCnfg = 3;    // FE-I4 hit discrimination setting
  // Do a check on IBL Slink ID
  eformat::helper::SourceIdentifier sid_rob(robId);
  sLinkSourceId = (sid_rob.module_id()) & 0x000F; // retrieve S-Link number from the source Identifier (0xRRRL, L = Slink number)
  if ( isIBL( robId ) && sLinkSourceId>0x3) { // Check if SLink number for the IBL is correct!
    generalwarning("In ROB 0x" << std::hex << robId <<  ": IBL/DBM SLink number not in correct range (0-3): SLink = " << std::dec << sLinkSourceId);
  }
  // Store length of module fragments, for monitoring
  uint32_t nwords_in_module_fragment = 0;
  bool corruptionError = false;
  // To check duplicated pixels
  std::unordered_set<Identifier> foundPixels;
  // ============== ============== ============== ============== ============== ============== ============== ============== //
  // loop over the data in the fragment
  const int nwords = robFrag->rod_ndata();
  for (int dataword_it = 0; dataword_it < nwords; ++dataword_it) {
    const uint32_t rawDataWord = vint[dataword_it];
    corruptionError = corruptionError || checkDataWordsCorruption( rawDataWord );
    uint32_t word_type = getDataType(rawDataWord, link_start);   // get type of data word
    ++nwords_in_module_fragment;
    serviceCodeCounter = 0;
    serviceCode = 0;

    int indexModule = -1;
    int indexFE     = -1;
//    if (word_type!=PRB_LINKHEADER) {
      if (offlineIdHash!=0xffffffff) { // now write the error word to the service
        isIBLModule = false;
        isDBMModule = false;
        if (m_is_ibl_present) {
          isIBLModule = isIBL(robId);
          isDBMModule = isDBM(robId);
        }

        int chFE = 0;
        if (isIBLModule || isDBMModule) {   // get FE channel id for IBL
          linkNum_IBLheader = decodeModule_IBL(rawDataWord);
          chFE = (extractFefromLinkNum(linkNum_IBLheader) & 0x1);
          if (m_pixelReadout->getModuleType(m_pixel_id->wafer_id(offlineIdHash))==InDetDD::PixelModuleType::IBL_3D) { chFE=0; }
        }
        else {                              // for PIXEL
          chFE = decodeFE(rawDataWord);
        }
        indexModule = static_cast<int>(offlineIdHash); // [0,2047]
        indexFE     = (1+chFE)*m_pixel_id->wafer_hash_max()+static_cast<int>(offlineIdHash);  // index for IDCInDetBSErrContainer
      }
//    }

    bool headerErr_decoding = false;
    bool headerErr_flagged  = false;
    bool headerErr_preamble = false;
    bool headerErr_timeout  = false;
    bool headerErr_lvl1id   = false;
    bool headerErr_bcid     = false;
    switch (word_type) { // depending on type of word call the right decoding method
      //-------------------------------------------------------------------------------------------- HEADER WORD FOUND
      case PRB_LINKHEADER:   // link (module) header found

        ATH_MSG_DEBUG("Header word found");
        errorRecoverable = 0; // reset errorcode

        if (link_start) {    // if header found before last header was closed by a trailer -> error
          generalwarning("In ROB 0x" << std::hex << robId << ": Unexpected link header found: 0x" << std::hex << rawDataWord << ", data corruption" << std::dec);
          m_numDecodingErrors++;
          headerErr_decoding = true;
        }
        else {
          ATH_MSG_DEBUG( "Header decoding starts" );
        }
        link_start = true;   // setting link (module) header found flag
        are_4condensed_words = false;
        receivedCondensedWords = false;
        isIBLModule = false;
        isDBMModule = false;
        countHitCondensedWords = 0;
        nwords_in_module_fragment = 1;
        if (m_is_ibl_present) {
          isIBLModule = isIBL(robId);
          isDBMModule = isDBM(robId);
        }
        // Keep track of IDs for previous fragment before decoding the new values
        prevLVL1ID = mLVL1ID;
        prevBCID = mBCID;
        prevLinkNum = mLink;
        prevOffsetBCID_ROB_FE = offsetBCID_ROB_FE;

        if (isIBLModule || isDBMModule) { // decode IBL/DBM header word.
#ifdef PIXEL_DEBUG
          ATH_MSG_VERBOSE( "Decoding the IBL/DBM header word: 0x" << std::hex << rawDataWord << std::dec );
#endif

          // IBL header data format: 001nnnnnFLLLLLLLLLLLLLBBBBBBBBBB
          mBCID = decodeBCID_IBL(rawDataWord);   // decode IBL BCID: B
          mLVL1ID = decodeL1ID_IBL(rawDataWord);   // decode IBL LVL1ID: L
          linkNum_IBLheader = decodeModule_IBL(rawDataWord);   // decode IBL FE number on the S-Link (range [0,7], as 8 FEs are connected to 1 S-Link in IBL) and the S-Link number itself: n
          if (decodeFeI4Bflag_IBL(rawDataWord)) {
            m_numFlaggedErrors++;
            headerErr_flagged = true;
          }
          fe_IBLheader = extractFefromLinkNum(linkNum_IBLheader);
          mLink = fe_IBLheader; // this is used to retrieve the onlineId. It contains only the 3 LSBs of the nnnnn, indicating the number of FE w.r.t. the SLink
          sLinkHeader = extractSLinkfromLinkNum(linkNum_IBLheader); // this is used to check that the (redundant) info is correctly transmitted
          if (sLinkHeader != sLinkSourceId) {
            generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": SLink discrepancy: Slink number from SourceId = 0x" << std::hex << sLinkSourceId << ", number from link header = 0x" << sLinkHeader << std::dec);
          }
          // If decoding fragment from same FE as previous one, do LVL1ID and BCID checks
          if (mLink == prevLinkNum) {
            // Check that L1ID is the same for all fragments
            if (mLVL1ID != prevLVL1ID && prevLVL1ID != 0x3FFF) {
              lvl1id_bcid_warning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": frame header LVL1ID differs from previous one (current frame LVL1ID = 0x" << mLVL1ID << ", previous frame LVL1ID = 0x" << prevLVL1ID << ")" << std::dec);
            }
            // Check that BCIDs are consecutive
            if ((mBCID != prevBCID + 1) && prevBCID != 0x3FFF && prevBCID != mBCID_max_IBL) {
              lvl1id_bcid_warning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": frame header with non-consecutive BCID (current BCID = 0x" << mBCID << ", previous BCID = 0x" << prevBCID << ")" << std::dec);
            }
          } 
          else {  // If decoding new FE, check BCID offset
            offsetBCID_ROB_FE = static_cast<int>(mBCID) - robBCID;
            if (offsetBCID_ROB_FE != prevOffsetBCID_ROB_FE && (offsetBCID_ROB_FE != 0x3FFF && prevOffsetBCID_ROB_FE != 0x3FFF)) {
              lvl1id_bcid_warning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << std::dec << ": FE header with unexpected BCID offset" << " wrt to ROB header (offset = " << offsetBCID_ROB_FE << ", expected " << prevOffsetBCID_ROB_FE << " from ROB)");
            }
            // Check that first fragment from each FE starts at the same BCID
            if (mBCID != prevStartingBCID && prevStartingBCID != 0x3FFF) {
              lvl1id_bcid_warning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": FE BCID starts at different value than in previous FE (current BCID = 0x" << mBCID << ", prev starting BCID = 0x" << prevStartingBCID << ")" << std::dec);
            }
            prevStartingBCID = mBCID;
          }
          // Count number of headers received
          if (mLink < 0x8) { ++nFragmentsPerFE[mLink]; }
        } 
        else { // this is Pixels detector decode header word. Data format: 001PtlbxdnnnnnnnMMMMLLLLBBBBBBBB (NOT IBL OR PIXEL)
          ATH_MSG_VERBOSE( "Decoding Pixel header word: 0x" << std::hex << rawDataWord << std::dec );
          mBCID = decodeBCID(rawDataWord);   // decode Pixel BCID: B
          mLVL1ID = decodeL1ID(rawDataWord);   // decode Pixel LVL1ID: L
          mLink = decodeModule(rawDataWord);   // decode Pixel link (module): n
          const uint32_t headerError = decodeHeaderErrors(rawDataWord);   // get link (module) header errors
          if ( headerError != 0 )  {
            sc = StatusCode::RECOVERABLE;
            errorRecoverable = errorRecoverable | (headerError << 20); //encode error as HHHHMMMMMMMMFFFFFFFFTTTT for header, flagword, trailer errors

            { // only treatment for header errors now, FIXME
              if (headerError & (1 << 3)) {
                m_numPreambleErrors++;
                headerErr_preamble = true;
              }
              if (headerError & (1 << 2)) {
                m_numTimeOutErrors++;
                headerErr_timeout  = true;
              }
              if (headerError & (1 << 1)) {
                m_numLVL1IDErrors++;
                headerErr_lvl1id   = true;
              }
              if (headerError & (1 << 0)) {
                m_numBCIDErrors++;
                headerErr_bcid     = true;
              }
            }

          }
        } // end of detector type decoding
        // Get onlineId
        onlineId = pixCabling->getOnlineIdFromRobId(robId, mLink);
        if (onlineId == 0) {
          generalwarning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": Got invalid onlineId (= 0) in FE header - dataword = 0x" << rawDataWord);
        }

#ifdef PIXEL_DEBUG
        ATH_MSG_VERBOSE("In decoder: got onlineId 0x" << std::hex << onlineId );
#endif

        offlineIdHash = m_pixel_id->wafer_hash(pixCabling->find_entry_onoff(onlineId));
        if (offlineIdHash != previous_offlineIdHash) {
          mBCID_offset = mBCID;   // set BCID offset if this is first LVL1A
          // The index array is defined in PixelRawDataProviderTool::SizeOfIDCInDetBSErrContainer()
          // Here, 52736 is a separator beween error flags and isActive flags.
          bsErrWord[static_cast<int>(offlineIdHash)+52736] = 1;
        }
        mLVL1A = mBCID - mBCID_offset;   // calculate the LVL1A
        if (mLVL1A < 0) {
          if (isIBLModule || isDBMModule) { mLVL1A=mLVL1A+mBCID_max_IBL; }
          else                            { mLVL1A=mLVL1A+mBCID_max_pix; }
        }

        // Fill header error word per IBL FE
        if (isIBLModule || isDBMModule) {
          if (indexFE>-1) {
            if (headerErr_decoding) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
            if (headerErr_flagged)  { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Flagged);  }
            if (headerErr_preamble) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Preamble); }
            if (headerErr_timeout)  { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::TimeOut);  }
            if (headerErr_lvl1id)   { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::LVL1ID);   }
            if (headerErr_bcid)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::BCID);     }
          }
        }
        else {
          // Fill header error word per module
          if (headerErr_decoding) { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::Decoding); }
          if (headerErr_flagged)  { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::Flagged);  }
          if (headerErr_preamble) { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::Preamble); }
          if (headerErr_timeout)  { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::TimeOut);  }
          if (headerErr_lvl1id)   { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::LVL1ID);   }
          if (headerErr_bcid)     { PixelByteStreamErrors::addError(bsErrWord[static_cast<int>(offlineIdHash)],PixelByteStreamErrors::BCID);     }
        }


#ifdef PIXEL_DEBUG
        ATH_MSG_VERBOSE( "link header with BCID: " << mBCID << " LVL1ID: " << mLVL1ID << " LVL1A: " << mLVL1A << " link: " << mLink << " found");
        ATH_MSG_VERBOSE( "got OfflineIDHash: 0x" << std::hex << offlineIdHash << " from OnlineID: 0x" <<  onlineId << std::dec << " link: " << mLink);
#endif

        if (offlineIdHash == 0xffffffff) {   // if link (module) online identifier (ROBID and link number) not found by mapping
          generalwarning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": Unknown OnlineId identifier in FE header - not found by mapping" << std::dec);
          m_numDecodingErrors++;
          link_start = false;   // resetting link (module) header found flag
          continue;   // skip this word and process next one
        }
        break;
        //-------------------------------------------------------------------------------------------- HIT DATA WORD FOUND
      case PRB_DATAWORD: // hit data found
        if (link_start) {   // if header found before hit -> expected
          ATH_MSG_DEBUG("Hit word found");
          unsigned int mFE = 0;
          unsigned int mRow = 0;
          unsigned int mColumn = 0;
          int mToT = 0;
          uint32_t row[nHits], col[nHits];
          int tot[nHits];
          if (isIBLModule || isDBMModule) { // check all the different types of IBL/DBM hit words (condensed and not condensed)

#ifdef PIXEL_DEBUG
            ATH_MSG_DEBUG ("Decoding IBL/DBM hit word: 0x" << std::hex << rawDataWord << std::dec);
#endif

            if (((rawDataWord & PRB_DATAMASK) == PRB_FIRSTHITCONDENSEDWORD) && (countHitCondensedWords == 0)) { // it is the first of the 4 hit condensed words (IBL/DBM)
              countHitCondensedWords = 1;
              are_4condensed_words = false; // All 4 words have to be found before turning it on
              condensedWords[0] = rawDataWord;
              continue;
            } 
            else if (countHitCondensedWords == 1) { // it is the second of the 4 hit condensed words (IBL)
              ++countHitCondensedWords;
              condensedWords[1] = rawDataWord;
              continue;
            } 
            else if (countHitCondensedWords == 2) { // it is the third of the 4 hit condensed words (IBL)
              ++countHitCondensedWords;
              condensedWords[2] = rawDataWord;
              continue;
            } 
            else if (((rawDataWord & PRB_DATAMASK) == PRB_LASTHITCONDENSEDWORD) && (countHitCondensedWords == 3)) { // it is the fourth of the 4 hit condensed words (IBL)
              are_4condensed_words = true;
              receivedCondensedWords = true;
              condensedWords[3] = rawDataWord;
              uint32_t hitwords[nHits];
              //Condensed words
              //  101RRRRRTTTTTTTTCCCCCCCRRRRRRRRR
              //  1CCCRRRRRRRRRTTTTTTTTCCCCCCCRRRR
              //  1TTTCCCCCCCRRRRRRRRRTTTTTTTTCCCC
              //  111TTTTTTTTCCCCCCCRRRRRRRRRTTTTT
              // They become 5 hitwords with format : TTTTTTTTCCCCCCCRRRRRRRRR
              // first hitword: simple mask of the first condensed word
              hitwords[0] = (condensedWords[0] & mask24);
              // second hitword: takes the "RRRRR" from the first condensed word + the "TTTTTTTTCCCCCCCRRRR" from the second (shifted left 5 bits to accomodate the "RRRRR" as LSB)
              hitwords[1] = ((condensedWords[0] >> (skipRow + skipCol+ skipTOT) ) & mask5) + ((condensedWords[1] & mask19) << skip5);
              // third hitword: takes remaining "CCCRRRRRRRRR" from the second condensed word + the "TTTTTTTTCCCC" from the third (shifted left 12 bits to accomodate the LSB coming from the second condensed word
              hitwords[2] = ((condensedWords[1] >> (skip4 + skipCol + skipTOT)) & mask12) + ((condensedWords[2] & mask12) << (skipRow + skip3));
              // fourth hitword: takes remaning "TTTCCCCCCCRRRRRRRRR" from the third condensed word +  the "TTTTT" from the fourth (shifted left 19 bits to accomodate the LSB coming from the third condensed word
              hitwords[3] = ((condensedWords[2] >> (skipTOT + skip4)) & mask19) + ((condensedWords[3] & mask5) << (skip3 + skipRow + skipCol));
              // fifth hitword: simply shift 5 right to eliminate the "TTTTT" and mask24
              hitwords[4] = (condensedWords[3] >> skip5) & mask24;
              
#ifdef PIXEL_DEBUG
              ATH_MSG_VERBOSE( "4 consecutive IBL hit words found. Condensed hits are being decoded");
#endif

              for (unsigned int i(0); i < nHits; ++i) {
                row[i] = divideHits (hitwords[i], 0, 8);
                col[i] = divideHits (hitwords[i], 9, 15);
                tot[i] = divideHits (hitwords[i], 16, 23);
                ATH_MSG_VERBOSE( "hitword[" << i << "] = 0x" << std::hex << hitwords[i] << ",  row: 0x" << row[i] << ",  col: 0x" << col[i] << ",  8-bit ToT: 0x" << tot[i] << std::dec);
              }
              countHitCondensedWords = 0;
              linkNum_IBLword = linkNum_IBLheader;
            } // end of the condensed hit words (IBL)
            else if ((rawDataWord & PRB_DATAMASK) == PRB_DATAWORD) { // it's a IBL not-condensed hit word
              if (countHitCondensedWords != 0) { // I received some IBL words, but less than 4, there was an error in the rod transmission
                generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": Interruption of IBL condensed words - hit(s) ignored (current dataword: 0x" << std::hex << rawDataWord << std::dec << ")");
                m_numDecodingErrors++;
                if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Decoding); }
                if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
                countHitCondensedWords = 0;
              }
              // IBL non-condensed Hit Words decoding:
              are_4condensed_words = false;
              mRow = decodeRow(rawDataWord, isIBLModule || isDBMModule );
              mColumn = decodeColumn(rawDataWord, isIBLModule || isDBMModule);
              mToT = decodeToT(rawDataWord, isIBLModule || isDBMModule); // this needs, later in the code, to be divided into Tot(0) and Tot(1), because the IBL TOT information can give information of the ToTs from up to two adjacent pixels
              linkNum_IBLword = decodeLinkNumHit_IBL(rawDataWord); // this corresponds to the "nnnnn" bits of the bytestream convention.
              //They represent: 3 LSB = number of the FE linked by the S-Link sending the block (range [0,7]); 2 MSB = S-Link number (range [1,4])
              if (linkNum_IBLword != linkNum_IBLheader) {
                generalwarning("In ROB 0x" << std::hex << robId << ": Link number mismatch - nnnnn (hit word) = 0x" << linkNum_IBLword << ", nnnnn (header) = 0x" << linkNum_IBLheader << std::dec);
              }

#ifdef PIXEL_DEBUG
              ATH_MSG_DEBUG("hit dataword found for module offlineIDHash: " << offlineIdHash << " Row: 0x" << std::hex << mRow << " Col: 0x" <<  mColumn  << " Tot: 0x" << mToT << std::dec);
#endif

            } 
            else { // it is an IBL/DBM hit word, but it hasn't been recognised
              generalwarning("In ROB 0x" << std::hex << robId << ", FE: 0x" << mLink << ": IBL/DBM hit word 0x" << rawDataWord << " not recognised" << std::dec);
              m_numDecodingErrors++;
              if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
              continue;
            }
            // computing the FE number on the silicon wafer for IBL ( mFE is = 0 for IBL 3D and the lower-eta FE on the planar sensor,
            // while is = 1 for the higher-eta FE on the planar sensor)
            mFE = getLocalFEI4(fe_IBLheader, onlineId);
            if (m_pixelReadout->getModuleType(pixCabling->find_entry_onoff(onlineId))==InDetDD::PixelModuleType::IBL_3D) {
              mFE = 0;
            }

          } // end of the if (isIBLModule || isDBMModule)
          else { // Pixel Hit Case

#ifdef PIXEL_DEBUG
            ATH_MSG_VERBOSE("Decoding Pixel hit word: 0x" << std::hex << rawDataWord << std::dec);
#endif

            if (countHitCondensedWords != 0){ // I received some IBL words, but less than 4, there was an error in the rod transmission
              generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink
                  << ": Interruption of IBL condensed words - hit(s) ignored (current dataword: 0x"
                  << std::hex << rawDataWord << std::dec << ")");
              m_numDecodingErrors++;
              if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Decoding); }
              if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
              countHitCondensedWords = 0;
            }
            //Pixel Hit Words decoding:
            are_4condensed_words = false;
            mFE = decodeFE(rawDataWord);
            mRow = decodeRow(rawDataWord, isIBLModule || isDBMModule );
            mColumn = decodeColumn(rawDataWord, isIBLModule || isDBMModule);
            mToT = decodeToT(rawDataWord, isIBLModule || isDBMModule);
            
#ifdef PIXEL_DEBUG
            ATH_MSG_DEBUG( "hit dataword found for module offlineIDHash: " << offlineIdHash << " FE: " << mFE << " Row: " << mRow << " Col: " <<  mColumn );   // hit found debug message
#endif

          } // end of the condition "!isIBLModule"
          // ------------------------
          // Create pixel collections

          // Decode only modules in provided list of hashes (i.e. not doing full scan)
          // Disabled for now.
          if (offlineIdHash == skipHash) continue;
          if (offlineIdHash != lastHash) {
            lastHash = offlineIdHash;
            if (vecHash!=nullptr) {
              std::vector<IdentifierHash>::iterator p = std::find(vecHash->begin(),vecHash->end(),offlineIdHash);
              if (p == vecHash->end()) { // is the Hash to be skipped (not in the request list)?
                skipHash = offlineIdHash;
                continue;
              }
            }
            ATH_CHECK( rdoIdc->naughtyRetrieve(offlineIdHash, coll) );//Returns null if not present
            if( coll == nullptr ) {
              coll = new PixelRawCollection (offlineIdHash);
              // get identifier from the hash, this is not nice
              Identifier ident = m_pixel_id->wafer_id(offlineIdHash);
              // set the Identifier to be nice to downstream clients
              coll->setIdentifier(ident);

              StatusCode sc = rdoIdc->addCollection(coll, offlineIdHash);
              ATH_MSG_DEBUG("Adding " << offlineIdHash);
              if (sc.isFailure()){
                ATH_MSG_ERROR("failed to add Pixel RDO collection to container" );
              }
            }
          }
          // ------------------------
          // Fill collections
          if (isIBLModule || isDBMModule) {
            unsigned int num_cycles_toFillCollection(0);
            if (are_4condensed_words) {
              num_cycles_toFillCollection = nHits;
            } 
            else {
              num_cycles_toFillCollection = 1;
              row[0] = mRow;
              col[0] = mColumn;
              tot[0] = mToT;
            }
            int IBLtot[2];
            for (unsigned int i(0); i < num_cycles_toFillCollection; ++i) {
            
#ifdef PIXEL_DEBUG
              ATH_MSG_VERBOSE( "ROW[" << i << "] = 0x"  << std::hex << row[i] << std::dec << ",  COL[" << i << "] = 0x" << std::hex << col[i] << std::dec << ",  8-bit TOT[" << i << "] = 0x" << std::hex << tot[i] << std::dec );
#endif

              // ToT1 equal 0 means no hit, regardless of HitDiscCnfg
              if ((tot[i] & 0xF0) == 0x00) {
                generalwarning("In ROB 0x" << std::hex << robId << ", FE 0x" << mLink << ": First IBL ToT field is 0 - hit word 0x" << rawDataWord << " decodes to ToT1 = 0" << (tot[i] >> 4) << ", ToT2 = 0x" << (tot[i] & 0xF) << ", row = " << std::dec << row[i] << " col = " << col[i] << std::dec);
                continue;
              } 
              else {
                if (((row[i] != 0) && (col[i] != 0) && (row[i] <= 336) && (col[i] <= 80)) || isDBMModule) { // FIXME : Hardcoded numbers. Should be ok as they are features of the FE-I4, but worth checking!
                  pixelId = m_pixelReadout ->getPixelIdfromHash (offlineIdHash, mFE, row[i], col[i]);

#ifdef PIXEL_DEBUG
                  ATH_MSG_VERBOSE( "PixelId: " << pixelId );
                  int eta_i = m_pixel_id->eta_index(pixelId);
                  int phi_i = m_pixel_id->phi_index(pixelId);
                  int eta_m = m_pixel_id->eta_module(pixelId);
                  int phi_m = m_pixel_id->phi_module(pixelId);
                  ATH_MSG_VERBOSE( "          eta_i: " << eta_i << ", phi_i: " << phi_i << ",  eta_m: " <<  eta_m << ", phi_m: " << phi_m );
#endif
                  if (pixelId == invalidPixelId) {
                    generalwarning("In ROB 0x" << std::hex << robId << ", FE 0x" << mLink << ": Illegal pixelId - hit word 0x" << rawDataWord << " decodes to ToT1 = 0" << (tot[i] >> 4) << ", ToT2 = 0x" << (tot[i] & 0xF) << ", row = " << std::dec << row[i] << " col = " << col[i] << std::dec);
                    m_numInvalidIdentifiers++;
                    if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
                    continue;
                  }
                  IBLtot[0] = divideHits(tot[i], 4, 7); // corresponds to (col, row)
                  IBLtot[1] = divideHits(tot[i], 0, 3); // corresponds to (col, row+1)
#ifdef PIXEL_DEBUG
                  ATH_MSG_VERBOSE("Starting from tot = 0x" << std::hex << tot[i] << " IBLtot[0] = 0x" << std::hex << IBLtot[0] << "   IBLtot[1] = 0x" << IBLtot[1] << std::dec );
#endif

                  if (!pixHitDiscCnfg) {
                    pixHitDiscCnfg = std::make_unique<SG::ReadCondHandle<PixelHitDiscCnfgData> > (m_condHitDiscCnfgKey);
                  }
                  // Get the hit discrimination configuration setting for this FE
                  if (m_pixelReadout->getModuleType(pixelId) == InDetDD::PixelModuleType::IBL_PLANAR || m_pixelReadout->getModuleType(pixelId) == InDetDD::PixelModuleType::DBM) {
                    hitDiscCnfg = (*pixHitDiscCnfg)->getHitDiscCnfgPL();
                  } 
                  else if (m_pixelReadout->getModuleType(pixelId) == InDetDD::PixelModuleType::IBL_3D) {
                    hitDiscCnfg = (*pixHitDiscCnfg)->getHitDiscCnfg3D();
                  }
                  // Now do some interpreting of the ToT values
                  if (hitDiscCnfg == 2 && IBLtot[0] == 2) IBLtot[0] = 16;
                  if (hitDiscCnfg == 2 && IBLtot[1] == 2) IBLtot[1] = 16;
                  if (not m_checkDuplicatedPixel or thisRdoIsUnique(pixelId, foundPixels)) {
                    // Insert the first part of the ToT info in the collection
                    coll->push_back(new RDO(pixelId, IBLtot[0], mBCID, mLVL1ID, mLVL1A));
                  } 
                  else {
                    m_numDuplicatedPixels++;
                  }

#ifdef PIXEL_DEBUG
                  ATH_MSG_VERBOSE( "Collection filled with pixelId: " << pixelId << " TOT = 0x" << std::hex << IBLtot[0] << std::dec << " mBCID = " << mBCID << " mLVL1ID = " << mLVL1ID << "  mLVL1A = " << mLVL1A );
#endif

                  // Second ToT field:
                  // For HitDiscCnfg 0-2, consider 0 to indicate no hit
                  // For HitDiscCng = 3, consider also F to indicate no hit
                  if (IBLtot[1] != 0x0 && (IBLtot[1] != 0xF || hitDiscCnfg != 3)) {
                    if ((row[i] + 1) > 336) { // FIXME: hardcoded number - but it should still be ok, because it's a feature of the FE-I4!
                      // this should never happen. If row[i] == 336, (row[i]+1) == 337. This row does not exist, so the TOT(337) should always be 0 (== no hit)
                      generalwarning("In ROB = 0x" << std::hex << robId << ", link 0x" << mLink << ": Illegal IBL row number for second ToT field, hit word 0x" << rawDataWord << " decodes to row = " << std::dec << row[i]+1 << " col = " << col[i] << " (ToT1 = 0x" << std::hex << IBLtot[0] << " ToT2 = 0x" << IBLtot[1] << ")");
                      if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
                      m_numInvalidIdentifiers++;
                      continue;
                    } 
                    else {
                      pixelId = m_pixelReadout->getPixelIdfromHash (offlineIdHash, mFE, row[i] + 1, col[i]);

#ifdef PIXEL_DEBUG
                      int eta_i = m_pixel_id->eta_index(pixelId);
                      int phi_i = m_pixel_id->phi_index(pixelId);
                      int eta_m = m_pixel_id->eta_module(pixelId);
                      int phi_m = m_pixel_id->phi_module(pixelId);
                      int phi_i_max = m_pixel_id->phi_index_max(pixelId);
                      int eta_i_max = m_pixel_id->eta_index_max(pixelId);
                      ATH_MSG_VERBOSE( "PixelId = " << pixelId);
                      ATH_MSG_VERBOSE( "          eta_i: " << eta_i << ", phi_i: " << phi_i << ",  eta_m: " <<  eta_m << ", phi_m: " << phi_m  << ", eta_i_max: " << eta_i_max << ", phi_i_max: " << phi_i_max );
#endif

                      if (pixelId == invalidPixelId) {
                        generalwarning("In ROB 0x" << std::hex << robId << ", FE 0x" << mLink << ": Illegal pixelId - hit word 0x" << rawDataWord << " decodes to ToT1 = 0" << (tot[i] >> 4) << ", ToT2 = 0x" << (tot[i] & 0xF) << ", row = " << std::dec << row[i] << " col = " << col[i] << std::dec);
                        m_numInvalidIdentifiers++;
                        if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
                        continue;
                      }
                      if (not m_checkDuplicatedPixel or thisRdoIsUnique(pixelId, foundPixels)) {
                        coll->push_back(new RDO(pixelId, IBLtot[1], mBCID, mLVL1ID, mLVL1A));
                      } 
                      else {
                        m_numDuplicatedPixels++;
                      }
                      
#ifdef PIXEL_DEBUG
                      ATH_MSG_VERBOSE( "Collection filled with pixelId: " << pixelId << " TOT = 0x" << std::hex << IBLtot[1] << std::dec << " mBCID = " << mBCID << " mLVL1ID = " << mLVL1ID << "  mLVL1A = " << mLVL1A );
#endif

                    }
                  }
                } // end of the if that checks that Row and Column of the IBL have a value > 0 and smaller than the maximum of the FE-I4.
                else {
                  generalwarning("In ROB = 0x" << std::hex << robId << ", link 0x" << mLink << ": Illegal IBL row/col number, hit word 0x" << rawDataWord << " decodes to row = " << std::dec << row[i] << " col = " << col[i] << " (ToT1 = 0x" << std::hex << (tot[i] >> 4) << " ToT2 = 0x" << (tot[i] & 0xF) << ")");
                  m_numInvalidIdentifiers++;
                  if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
                  continue;
                }
              } // end of the else from the condition "(tot[i] & 0xF0) == 0xF", corresponding to an invalid tot sent by the ROD
            } // end of the loop over the hit information coming from the IBL words (condensed or not)
          } // end IBL/DBM words pushed into Collection

          else { // it is a Pixel hit word to be saved into the Collection
            pixelId = m_pixelReadout->getPixelIdfromHash(offlineIdHash, mFE, mRow, mColumn);
            if (pixelId == invalidPixelId) {
              ATH_MSG_DEBUG("In ROB = 0x" << std::hex << robId << ", link 0x" << mLink << ": Illegal pixelId - row = " << std::dec << mRow << ", col = " << mColumn << ", dataword = 0x" << std::hex << rawDataWord << std::dec);
              m_numInvalidIdentifiers++;
              if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Invalid); }
              if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
              continue;
            }
            // Now the Collection is there for sure. Create RDO and push it into Collection.
            if (not m_checkDuplicatedPixel or thisRdoIsUnique(pixelId, foundPixels)) {
              coll->push_back(new RDO(pixelId, mToT, mBCID, mLVL1ID, mLVL1A));
            } 
            else {
              m_numDuplicatedPixels++;
            }

#ifdef PIXEL_DEBUG
            ATH_MSG_VERBOSE( "Collection filled with pixelId: " << pixelId << " TOT = 0x" << std::hex << mToT << std::dec << " mBCID = " << mBCID << " mLVL1ID = " << mLVL1ID << "  mLVL1A = " << mLVL1A );
#endif

          } // end Pixel words pushed into Collection
        } // end of the if (link_start)
        else { // no header found before hit -> error
          generalwarning("In ROB = 0x" << std::hex << robId << ", link 0x" << mLink << ": Unexpected hit dataword: " << rawDataWord << " - hit ignored");
          m_numDecodingErrors++;
          if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Decoding); }
          if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
          continue;
        }
        break;
      //-------------------------------------------------------------------------------------------- TRAILER WORD FOUND
      case PRB_LINKTRAILER:  // link (module) trailer found
        if (link_start) {   // if header found before trailer -> expected
          ATH_MSG_DEBUG( "link trailer found" );   // trailer found debug message
        } 
        else {   // no header found before trailer -> error
          generalwarning("In ROB = 0x" << std::hex << robId << ", link 0x" << mLink << ": Unexpected trailer found: 0x" << std::hex << rawDataWord << ", data corruption");
          m_numDecodingErrors++;
          if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Decoding); }
          if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }
          continue;
        }
        previous_offlineIdHash = offlineIdHash;  // save offlineId for next header;
        link_start = false;   // resetting link (module) header found flag
        are_4condensed_words = false;
        // Trailer decoding and error handling
        if (isIBLModule || isDBMModule) { // decode IBL/DBM Trailer word: 010nnnnnECPplbzhvMMMMMMMMMMxxxxx
          ATH_MSG_VERBOSE( "Decoding IBL/DBM trailer word: 0x" << std::hex << rawDataWord << std::dec );
          const uint32_t trailerError = decodeTrailerErrors_IBL(rawDataWord); // => E cPpl bzhv // taking all errors together.

          if (trailerError & (1 << 8)) { // time out error bit => E
            m_numTimeOutErrors++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::TimeOut); }
          }

          if (trailerError & (1 << 6)) {// link masked by PPC => P
            m_numLinkMaskedByPPC++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::LinkMaskedByPPC); }
          }
          if (trailerError & (1 << 5)) { // preamble error bit => p
            m_numPreambleErrors++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Preamble); }
          }
          if (trailerError & (1 << 4)) { // LVL1 error bit => l
            m_numLVL1IDErrors++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::LVL1ID); }
          }
          if (trailerError & (1 << 3)) {// BCID error bit => b
            m_numBCIDErrors++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::BCID); }
          }
          if (trailerError & (1 << 2)) { // trailer error bit => z
            m_numTrailerErrors++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Trailer); }
          }
          if (trailerError & (1 << 1)) { // header/trailer limit error=> h
            m_numLimitError++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Limit); }
          }
          if (trailerError & (1 << 0)) { // data overflow error=> v
            m_numInvalidIdentifiers++;
            if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Invalid); }
          }

          if (trailerError & (1 << 7)) { // condensed mode bit => W
            if (!receivedCondensedWords) { generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": condensed mode bit is set, but no condensed words received" << std::dec); }	
          }
          //At least temporarily removed because the data format is not clear (Franconi, 17.06.2014)
          linkNum_IBLword = decodeLinkNumTrailer_IBL(rawDataWord); // => n
          // Do checks on info in trailer vs header
          if (linkNum_IBLword != linkNum_IBLheader) {
            generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink <<  ": Link number mismatch - nnnnn (trailer) = 0x" << linkNum_IBLword << ", nnnnn (header) = 0x" << linkNum_IBLheader << std::dec);
          }
          if (decodeBcidTrailer_IBL(rawDataWord) != (mBCID & PRB_BCIDSKIPTRAILERmask_IBL)) {
            generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": Trailer BCID does not match header (trailer BCID = 0x" << decodeBcidTrailer_IBL(rawDataWord) << ", 5LSB of header BCID = 0x" << (mBCID & PRB_BCIDSKIPTRAILERmask_IBL) << ")" << std::dec);
          }
        } 
        else { // decode Pixel trailer word
          ATH_MSG_VERBOSE( "Decoding Pixel trailer word: 0x" << std::hex << rawDataWord << std::dec );
          const uint32_t trailerError = decodeTrailerErrors(rawDataWord);   // creating link (module) trailer error variable
          if (trailerError != 0) {
            sc = StatusCode ::RECOVERABLE;
            errorRecoverable = errorRecoverable | trailerError; //encode error as HHHHMMMMMMMMFFFFFFFFTTTT for header, flagword, trailer errors

            const std::array<unsigned long long, 4> bitPosition{0_BIT, 1_BIT, 2_BIT, 3_BIT};
            for (const auto bit:bitPosition){
              if (trailerError & bit) ++m_numTrailerErrors;
            }

            if ( trailerError & 0xF ) {
              if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Trailer); }
              if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Trailer); }
            }

            // Extra for pixel monitoring
            if (trailerError != 0) {
              if (indexModule>-1) {
                if (trailerError & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::RODTrailerBitError); }
                if (trailerError & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::RODHeaderLimit); }
                if (trailerError & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::RODDataOVerflow); }
              }
              if (indexFE>-1) {
                if (trailerError & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::RODTrailerBitError); }
                if (trailerError & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::RODHeaderLimit); }
                if (trailerError & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::RODDataOVerflow); }
              }
            }

          }
        }
        break;
        //-------------------------------------------------------------------------------------------- FE ERROR TYPE 2 WORD FOUND
      case PRB_FEERROR2 :    // type 2 flag found

        are_4condensed_words = false;
        if (isIBLModule || isDBMModule) { // IBL flag found
          ATH_MSG_VERBOSE( "Decoding IBL/DBM FEflag word: 0x" << std::hex << rawDataWord << std::dec);
          linkNum_IBLword = decodeLinkNumFEFlag_IBL(rawDataWord);
          if (linkNum_IBLword != linkNum_IBLheader) {
            generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": Link number mismatch - nnnnn (error word) = 0x" << linkNum_IBLword << ", nnnnn (header) = 0x" << linkNum_IBLheader << std::dec);
          }
          serviceCodeCounter = decodeServiceCodeCounter_IBL(rawDataWord)+1; // frequency of the serviceCode (with the exceptions of serviceCode = 14,15 or 16)
          serviceCode = decodeServiceCode_IBL(rawDataWord); // is a code. the corresponding meaning is listed in the table in the FE-I4 manual, pag. 105
          // Treat the service code only if its meaning is valid (i.e. value < 31)
          if (serviceCode>31) {
            generalwarning("In ROB 0x" << std::hex << robId << ", link 0x" << mLink << ": Got out-of-bounds service code: " << std::dec << serviceCode << " (counter: " << serviceCodeCounter << "), ignored");
          }

          int chFE = (extractFefromLinkNum(linkNum_IBLheader) & 0x1);
          if (m_pixelReadout->getModuleType(m_pixel_id->wafer_id(offlineIdHash))==InDetDD::PixelModuleType::IBL_3D) { chFE=0; }

          if (serviceCodeCounter>0 && serviceCode<32) {
            if (serviceCode!=14) {
              // Monitor service record for IBL (not DBM)
              if (static_cast<int>(offlineIdHash)>155 && static_cast<int>(offlineIdHash)<436) {
                // The index array is defined in PixelRawDataProviderTool::SizeOfIDCInDetBSErrContainer()
                int indexOffset = 17*m_pixel_id->wafer_hash_max();
                int indexSvcCounter = indexOffset+serviceCode*280*2+2*(static_cast<int>(offlineIdHash)-156)+chFE;

                bsErrWord[indexSvcCounter] = serviceCodeCounter;
              }
            }
          }

        } 
        else { // Pixel type2 flag found
          ATH_MSG_VERBOSE( "Decoding Pixel FEflag word: 0x" << std::hex << rawDataWord << std::dec );
          uint32_t FEFlags = decodeFEFlags2(rawDataWord);   // get FE flags
          uint32_t MCCFlags = decodeMCCFlags(rawDataWord);   // get MCC flags
          FEFlags = FEFlags & 0xff;
          if ((MCCFlags | FEFlags) != 0) {
            sc = StatusCode::RECOVERABLE;
            errorRecoverable = errorRecoverable | (MCCFlags << 12) | (FEFlags << 4); //encode error as HHHHMMMMMMMMFFFFFFFFTTTT for header, flagword, trailer errors
            //for now just sum all flagged errors_
            const std::array<unsigned long long, 8> bitPosition {0_BIT, 1_BIT, 2_BIT, 3_BIT, 4_BIT, 5_BIT, 6_BIT, 7_BIT};
            for (const auto thisBit:bitPosition){
              if (MCCFlags & thisBit) ++m_numFlaggedErrors;
              if (FEFlags & 0xf3 & thisBit) ++m_numFlaggedErrors;
            }
           
            if ( MCCFlags & 0xff or FEFlags & 0xf3 ) {
              if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Flagged); }
              if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Flagged); }
            }

            if (indexModule>-1) {
              if (MCCFlags & (1 << 7)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 6)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 5)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 4)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCLVL1IDEoECheck); }
              if (MCCFlags & (1 << 3)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCBCIDEoECheck); }
              if (MCCFlags & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCLVL1IDCheck); }
              if (MCCFlags & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCEoEOverflow); }
              if (MCCFlags & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::MCCHitOverflow); }

              if (FEFlags & (1 << 4)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::FEWarning); }
              if (FEFlags & (1 << 3)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::FEHitParity); }
              if (FEFlags & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::FERegisterParity); }
              if (FEFlags & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::FEHammingCode); }
              if (FEFlags & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::FEEoCOverflow); }
            }

            if (indexFE>-1) {
              if (MCCFlags & (1 << 7)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 6)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 5)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCUndefined); }
              if (MCCFlags & (1 << 4)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCLVL1IDEoECheck); }
              if (MCCFlags & (1 << 3)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCBCIDEoECheck); }
              if (MCCFlags & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCLVL1IDCheck); }
              if (MCCFlags & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCEoEOverflow); }
              if (MCCFlags & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::MCCHitOverflow); }

              if (FEFlags & (1 << 4)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::FEWarning); }
              if (FEFlags & (1 << 3)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::FEHitParity); }
              if (FEFlags & (1 << 2)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::FERegisterParity); }
              if (FEFlags & (1 << 1)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::FEHammingCode); }
              if (FEFlags & (1 << 0)) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::FEEoCOverflow); }
            }
          } 
          else {
            m_numDisabledFEErrors++;
            if (indexFE>-1) { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::BadFE); }
            ATH_MSG_DEBUG( "Disabled Pixel chip " << ((rawDataWord & 0x0F000000) >> 24) );
          }
        }

        break;

        //-------------------------------------------------------------------------------------------- RAWDATA WORD FOUND
      case PRB_RAWDATA : // add treatment for raw data word here!

        are_4condensed_words = false;

        break;

        //-------------------------------------------------------------------------------------------- UNKNOWN WORD
      default:
        generalwarning("In ROB 0x" << std::hex << robId << ", FE 0x" << mLink << ": Unknown word type found, 0x" << std::hex << rawDataWord << std::dec << ", ignoring");
        m_numDecodingErrors++;
        if (indexModule>-1) { PixelByteStreamErrors::addError(bsErrWord[indexModule],PixelByteStreamErrors::Decoding); }
        if (indexFE>-1)     { PixelByteStreamErrors::addError(bsErrWord[indexFE],PixelByteStreamErrors::Decoding); }

    } // end of switch
  }   // end of loop over ROD

  ATH_CHECK( updateEventInfoIfEventCorrupted( corruptionError ) );

  // Verify that all active IBL FEs sent the same number of headers
  if (isIBLModule || isDBMModule) {
    checkUnequalNumberOfHeaders( nFragmentsPerFE, robId );
  }

  // Fill error information in IDC
  for (size_t i=0; i<static_cast<size_t>(bsErrWord.size()); i++) {
    if (bsErrWord[i]>0) {
      decodingErrors.setOrDrop(i,bsErrWord[i]);
    }
  }

  if (sc == StatusCode::RECOVERABLE) {

    if (errorRecoverable == (3 << 20) ){  // Fix for M8, this error always occurs, masked out REMOVE FIXME !!
      m_masked_errors++;
      return StatusCode::SUCCESS;
    }
    ATH_MSG_DEBUG( "Recoverable errors found in PixelRodDecoder, errorcode: " << errorRecoverable );
  }
  return sc;
}

// ****************************************************************************************************** DECODING FUNCTIONS

////////////////////////
// extract the FE (w.r.t. SLink) from the Link Number (nnnnn), present in the IBL header, long hit, fe flag error and trailer words
////////////////////////
uint32_t PixelRodDecoder::extractFefromLinkNum(const uint32_t linkNum) const
{
  return linkNum & 0x7;
}

////////////////////////
// extract the FE (w.r.t. SLink) from the Link Number (nnnnn), present in the IBL header, long hit, fe flag error and trailer words
////////////////////////
uint32_t PixelRodDecoder::extractSLinkfromLinkNum(const uint32_t linkNum) const
{
  return (linkNum >> 3) & 0x3;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ HEADER BITS

////////////////////////
// decode IBL BCID (bits 0-9 of IBL module header word)
////////////////////////
uint32_t PixelRodDecoder::decodeBCID_IBL(const uint32_t word) const
{
  return ((word >> PRB_BCIDskip_IBL) & PRB_BCIDmask_IBL);
}

////////////////////////
// decode IBL LVL1ID (bits 10-22 of IBL module header word)
////////////////////////
uint32_t PixelRodDecoder::decodeL1ID_IBL(const uint32_t word) const
{
  return ((word >> PRB_L1IDskip_IBL) & PRB_L1IDmask_IBL);
}

////////////////////////
// decode IBL FeI4B flag bit (bit 23 of IBL module header word)
////////////////////////
uint32_t PixelRodDecoder::decodeFeI4Bflag_IBL(const uint32_t word) const
{
  return ((word >> PRB_FeI4BFLAGskip_IBL) & PRB_FeI4BFLAGmask_IBL);
}

////////////////////////
// decode IBL module link number  (bits 24-28 of IBL module header word)
////////////////////////
uint32_t PixelRodDecoder::decodeModule_IBL(const uint32_t word) const
{
  return ((word >> PRB_MODULEskip_IBL) & PRB_MODULEmask_IBL);
}

////////////////////////
// decode BCID (bits 0-7 of module Pixel header word)
////////////////////////
uint32_t PixelRodDecoder::decodeBCID(const uint32_t word) const
{
  return ((word >> PRB_BCIDskip) & PRB_BCIDmask);
}

////////////////////////
// decode LVL1ID (bits 8-11 of module Pixel header word)
////////////////////////
uint32_t PixelRodDecoder::decodeL1ID(const uint32_t word) const
{
  return ((word >> PRB_L1IDskip) & PRB_L1IDmask);
}

////////////////////////
// decode LVL1ID skipped (bits 12-15 of module Pixel header word)
////////////////////////
uint32_t PixelRodDecoder::decodeL1IDskip(const uint32_t word) const
{
  return ((word >> PRB_L1IDSKIPskip) & PRB_L1IDSKIPmask);
}

////////////////////////
// decode module link number (bits 16-22 of module Pixel header word)
////////////////////////
uint32_t PixelRodDecoder::decodeModule(const uint32_t word) const
{
  return ((word >> PRB_MODULEskip) & PRB_MODULEmask);
}

////////////////////////
// decode header errors (bits 25-28 of module Pixel header word)
////////////////////////
uint32_t PixelRodDecoder::decodeHeaderErrors(const uint32_t word) const
{
  return ((word >> PRB_HEADERERRORSskip) & PRB_HEADERERRORSmask);
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ HIT WORD BITS

////////////////////////
// decode hit row number (bits 0-7 of Pixel hit word, bits 0-8 of IBL hit word)
////////////////////////
uint32_t PixelRodDecoder::decodeRow(const uint32_t word, bool isIBLorDBM) const {
  if( isIBLorDBM ) {
    return ((word >> PRB_ROWskip_IBL) & PRB_ROWmask_IBL);
  }
  else {
    return ((word >> PRB_ROWskip) & PRB_ROWmask);
  }
}

////////////////////////
// decode hit column number (bits 8-12 of Pixel hit word, bits 9-15 of IBL hit word)
////////////////////////
uint32_t PixelRodDecoder::decodeColumn(const uint32_t word, bool isIBLorDBM) const {
  if( isIBLorDBM ) {
    return ((word >> PRB_COLUMNskip_IBL) & PRB_COLUMNmask_IBL);
  }
  else {
    return ((word >> PRB_COLUMNskip) & PRB_COLUMNmask);
  }
}

////////////////////////
// decode hit TimeOverThreshold value (bits 16-23 of Pixel hit word, bits 16-23 of IBL hit word)
////////////////////////
uint32_t PixelRodDecoder::decodeToT(const uint32_t word, bool isIBLorDBM) const {
  if ( isIBLorDBM ) { // IBL Tot => 8 bits to be divided into 2 Tot info, each of which is 4 bits
    uint32_t Tot8Bits = ((word >> PRB_TOTskip_IBL) & PRB_TOTmask_IBL); // this Tot information returns a 8-bit word. In the case of the IBL, it will be splitted into Tot(1) and Tot(2), corresponding to pixels (col,row) and (col, row+1)
    return Tot8Bits;
  }
  else { // Pixel Tot
    return ((word >> PRB_TOTskip) & PRB_TOTmask);
  }
}

////////////////////////
// decode Link number in the IBL not-condensed hit word (bits 24-26 of IBL hit word)
////////////////////////
uint32_t PixelRodDecoder::decodeLinkNumHit_IBL(const uint32_t word) const {
  return ((word >> PRB_LINKNUMHITskip_IBL) & PRB_LINKNUMHITmask_IBL);
}

////////////////////////
// decode module FE number (bits 24-27 of Pixel hit word)
////////////////////////
uint32_t PixelRodDecoder::decodeFE(const uint32_t word) const {
  return ((word >> PRB_FEskip) & PRB_FEmask);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ TRAILER BITS

////////////////////////
// decode IBL "M" in the trailer (bits 5-14 of IBL module trailer word)
////////////////////////
uint32_t PixelRodDecoder::decodeSkippedTrigTrailer_IBL(const uint32_t word) const {
  return ((word >> PRB_SKIPPEDTRIGGERTRAILERskip_IBL) & PRB_SKIPPEDTRIGGERTRAILERmask_IBL);
}

////////////////////////
// decode IBL trailer errors (bits 15-23 of IBL module trailer word)
////////////////////////
uint32_t PixelRodDecoder::decodeTrailerErrors_IBL(const uint32_t word) const {
  return ((word >> PRB_TRAILERERRORSskip_IBL) & PRB_TRAILERERRORSmask_IBL);
}

////////////////////////
// decode Link Number in the IBL Trailer (bits 24-26 of IBL module trailer word)
////////////////////////
uint32_t PixelRodDecoder::decodeLinkNumTrailer_IBL(const uint32_t word) const {
  return ((word >> PRB_LINKNUMTRAILERskip_IBL) & PRB_LINKNUMTRAILERmask_IBL);
}


////////////////////////
// decode BCID in trailer word
////////////////////////
uint32_t PixelRodDecoder::decodeBcidTrailer_IBL(const uint32_t word) const {
  return ((word >> PRB_BCIDSKIPTRAILERskip_IBL) & PRB_BCIDSKIPTRAILERmask_IBL);
}


////////////////////////
// decode Condensed Mode Bit in the IBL trailer (bit 22 of IBL module trailer word)
////////////////////////
uint32_t PixelRodDecoder::decodeCondensedModeBit_IBL(const uint32_t word) const {
  return ((word >> PRB_CONDENSEDMODEskip_IBL) & PRB_CONDENSEDMODEmask_IBL);
}

////////////////////////
// decode trailer errors (bits 26-28 of Pixel module trailer word)
////////////////////////
//template <class ROBData>
uint32_t PixelRodDecoder::decodeTrailerErrors(const uint32_t word) const {
  return ((word >> PRB_TRAILERERRORSskip) & PRB_TRAILERERRORSmask);
}

////////////////////////
// look for bitflips in trailer word (bits 0-25 of Pixel module trailer word)
////////////////////////
uint32_t PixelRodDecoder::decodeTrailerBitflips(const uint32_t word) const {
  uint32_t flipword = ((word >> PRB_TRAILERBITFLIPskip) & PRB_TRAILERBITFLIPmask);
  uint32_t errorcount = 0;
  if (flipword != 0) {
    for (int i=0; i<=25; i++){    // this calculates the checksum of the flipword = the number of bitflips
      if (flipword & 1) errorcount++;
      flipword >>= 1;
    }
  }
  return errorcount;     // returns the number of errors
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ FLAG WORDS BITS

////////////////////////
// decodeServiceCodeCounter (bits 0-9 of IBL flag word)
////////////////////////
uint32_t PixelRodDecoder::decodeServiceCodeCounter_IBL(const uint32_t word) const {
  return((word >> PRB_SERVICECODECOUNTERskip_IBL) & PRB_SERVICECODECOUNTERmask_IBL);
}

////////////////////////
// decodeServiceCode (bits 17-22 of IBL flag word)
////////////////////////
uint32_t PixelRodDecoder::decodeServiceCode_IBL(const uint32_t word) const {
  return((word >> PRB_SERVICECODEskip_IBL) & PRB_SERVICECODEmask_IBL);
}

////////////////////////
// decodeLinkNumFEFlag (bits 25-27 of IBL flag word)
////////////////////////
uint32_t PixelRodDecoder::decodeLinkNumFEFlag_IBL(const uint32_t word) const {
  return((word >> PRB_LINKNUMBERFEFLAGskip_IBL) & PRB_LINKNUMBERFEFLAGmask_IBL);
}


////////////////////////
// decode FE flag errors (bits 0-7 of Pixel flag type 2 word)
////////////////////////
uint32_t PixelRodDecoder::decodeFEFlags2(const uint32_t word) const {
  return ((word >> PRB_FEFlagskip2) & PRB_FEFlagmask2);
}

////////////////////////
// decode MCC flag errors (bits 8-12 of Pixel flag type 2 word)
////////////////////////
uint32_t PixelRodDecoder::decodeMCCFlags(const uint32_t word) const {
  return ((word >> PRB_MCCFlagskip) & PRB_MCCFlagmask);
}

////////////////////////
// determine module word type
////////////////////////
uint32_t PixelRodDecoder::getDataType(unsigned int rawDataWord, bool link_start) const {
  if (link_start) { // there was a header, so if a timeout is found it's really a timeout and not an empty header, same for raw data
    if ((rawDataWord & PRB_DATAMASK) == PRB_RAWDATA ) return PRB_RAWDATA;   // module flag word found
    if (rawDataWord == PRB_TIMEOUT) return PRB_TIMEOUT;   // time out word found
  }
  if ((rawDataWord & PRB_DATAMASK) == PRB_LINKHEADER) return PRB_LINKHEADER;   // module header word found
  if (((rawDataWord & PRB_FEERRORMASK) == PRB_FEERROR2) && ((rawDataWord & PRB_FEERROR2CHECK) == PRB_FEERROR2CHECK)) return PRB_FEERROR2;   // type 2 word found Pixel
  if (((rawDataWord & PRB_DATAMASK) == PRB_FEERROR_IBL)) return PRB_FEERROR2;   // type word found IBL
  if ((rawDataWord & PRB_HITMASK) == PRB_DATAWORD) return PRB_DATAWORD;   // module hit word found
  if ((rawDataWord & PRB_DATAMASK) == PRB_LINKTRAILER) return PRB_LINKTRAILER;   // module trailer word found
  return PRB_UNKNOWNWORD;   // unknown word found
}


bool PixelRodDecoder::checkDataWordsCorruption( uint32_t rawDataWord) const {
    if ( ATH_UNLIKELY( rawDataWord==0xaa1234aa )) {
      generalwarning("Evt marker encountered during loop on ROD datawords");
      return true;
    } 
    else if ( ATH_UNLIKELY( rawDataWord==0xdd1234dd ) ){
      generalwarning("ROB marker encountered during loop on ROD datawords");
      return true;
    } 
    else if ( ATH_UNLIKELY( rawDataWord==0xee1234ee) ){
      generalwarning("ROD marker encountered during loop on ROD datawords");
      return true;
    }
    return false; // no corruption
}

uint32_t PixelRodDecoder::treatmentFEFlagInfo(unsigned int serviceCode, unsigned int serviceCodeCounter) const {
  unsigned int etc = 0, l1req = 0;
  // Log the code count
  if (serviceCode == 16) {                        // I'm not like those other codes
    etc = (serviceCodeCounter >> 4) & 0x1F;
    l1req = serviceCodeCounter & 0x7;
  }
  // Return a 19-bit code to be used for monitoring
  // Service codes of type 'inf' are excluded, except for codes 15, 16 and 25,
  // see FE-I4B manual table 39
  // The bit position of each monitored service code is sequential starting from the
  // least significant bit, see below
  uint32_t code = 0;
  switch (serviceCode) {
  case 0: // BCID counter error, from EOCHL
    ATH_MSG_DEBUG("BCID counter error (retrieved in IBL FE Flag word)");
    code = code | (1 << 0);
    break;

  case 1: // Hamming code error in word 0, from EOCHL
    ATH_MSG_DEBUG("Hamming code error in word 0 (retrieved in IBL FE Flag word)");
    code = code | (1 << 1);
    break;

  case 2: // Hamming code error in word 1, from EOCHL
    ATH_MSG_DEBUG("Hamming code error in word 1 (retrieved in IBL FE Flag word)");
    code = code | (1 << 2);
    break;

  case 3: // Hamming code error in word 2, from EOCHL
    ATH_MSG_DEBUG("Hamming code error in word 2 (retrieved in IBL FE Flag word)");
    code = code | (1 << 3);
    break;

  case 4: // L1_in counter error, from EOCHL
    ATH_MSG_DEBUG("L1_in counter error (retrieved in IBL FE Flag word)");
    code = code | (1 << 4);
    break;

  case 5: // L1 request counter error from EOCHL
    ATH_MSG_DEBUG("L1 request counter error (retrieved in IBL FE Flag word)");
    code = code | (1 << 5);
    break;

  case 6: // L1 register error from EOCHL
    ATH_MSG_DEBUG("L1 register error (retrieved in IBL FE Flag word)");
    code = code | (1 << 6);
    break;

  case 7: // L1 Trigger ID error from EOCHL
    ATH_MSG_DEBUG("L1 Trigger ID error (retrieved in IBL FE Flag word)");
    code = code | (1 << 7);
    break;

  case 8: // readout processor error from EOCHL
    ATH_MSG_DEBUG("Readout processor error (retrieved in IBL FE Flag word)");
    code = code | (1 << 8);
    break;

  case 9: // Fifo_Full flag pulsed from EOCHL
    ATH_MSG_DEBUG("Fifo_Full flag pulsed (retrieved in IBL FE Flag word)");
    break;

  case 10: // HitOr bus pulsed from PixelArray
    ATH_MSG_DEBUG("HitOr bus pulsed (retrieved in IBL FE Flag word)");
    break;

    // case 11 to 13 not used
  case 11:
  case 12:
  case 13:
    ATH_MSG_DEBUG("Received invalid service code: " << serviceCode << " (payload = " << serviceCodeCounter << ")");
    break;

  case 14: // 3 MSBs of bunch counter and 7 MSBs of L1A counter from EOCHL
    ATH_MSG_DEBUG("Received service code 14 - ignored");
    break;

  case 15: // Skipped trigger counter
    ATH_MSG_DEBUG("Skipped trigger counter (retrieved in IBL FE Flag word). There are " << serviceCodeCounter << " skipped triggers.");
    code = code | (1 << 9);
    break;

  case 16: // Truncated event flag and counter from EOCHL
    ATH_MSG_DEBUG("Truncated event flag: ETC = " << etc << ", L1Req = " << l1req);
    code = code | (1 << 10);
    break;

    // case 17 to 20 not used
  case 17:
  case 18:
  case 19:
  case 20:
    ATH_MSG_DEBUG("Received invalid service code: " << serviceCode << " (payload = " << serviceCodeCounter << ")");
    break;

  case 21: // Reset bar RA2b pulsedd from Pad, PRD
    ATH_MSG_DEBUG("Reset bar RA2b pulsed (retrieved in IBL FE Flag word)");
    break;

  case 22: // PLL generated clock phase faster than reference from CLKGEN
    ATH_MSG_DEBUG("PLL generated clock phase faster than reference (retrieved in IBL FE Flag word)");
    break;

  case 23: // Reference clock phase faster than PLL, from CLKGEN
    ATH_MSG_DEBUG("Reference clock phase faster than PLL (retrieved in IBL FE Flag word)");
    break;

  case 24: // Triple redundant mismatchfrom CNFGMEM
    ATH_MSG_DEBUG("Triple redundant mismatch (retrieved in IBL FE Flag word)");
    code = code | (1 << 11);
    break;

  case 25: // Write register data error from CMD
    ATH_MSG_DEBUG("Write register data error (retrieved in IBL FE Flag word)");
    code = code | (1 << 12);
    break;

  case 26: // Address error from CMD
    ATH_MSG_DEBUG("Address error (retrieved in IBL FE Flag word)");
    code = code | (1 << 13);
    break;

  case 27: // Other command decoder error- see CMD section from CMD
    ATH_MSG_DEBUG("Other command decoder error- see CMD section (retrieved in IBL FE Flag word)");
    code = code | (1 << 14);
    break;

  case 28: // Bit flip detected in command decoder input stream from CMD
    ATH_MSG_DEBUG("Bit flip detected in command decoder input stream (retrieved in IBL FE Flag word)");
    code = code | (1 << 15);
    break;

  case 29: // SEU upset detected in command decoder (triple redundant mismatch) from CMD
    ATH_MSG_DEBUG("SEU upset detected in command decoder (triple redundant mismatch) (retrieved in IBL FE Flag word)");
    code = code | (1 << 16);
    break;

  case 30: // Data bus address error from CMD
    ATH_MSG_DEBUG("Data bus address error (retrieved in IBL FE Flag word)");
    code = code | (1 << 17);
    break;

  case 31: // Triple redundant mismatch from CMD
    ATH_MSG_DEBUG("Triple redundant mismatch (retrieved in IBL FE Flag word)");
    code = code | (1 << 18);
    break;

  default:
    ATH_MSG_DEBUG("ServiceCode not used at the moment");

  }
  return code;
}

StatusCode PixelRodDecoder::updateEventInfoIfEventCorrupted( bool isCorrupted ) const {
  if ( not isCorrupted )
    return StatusCode::SUCCESS;
  //Set EventInfo error
  const xAOD::EventInfo* eventInfo=nullptr;
  ATH_CHECK(evtStore()->retrieve(eventInfo));
  if (!eventInfo->updateErrorState(xAOD::EventInfo::Pixel,xAOD::EventInfo::Error)) {
    ATH_MSG_WARNING(" cannot set EventInfo error state for Pixel " );
  }
  if (!eventInfo->updateEventFlagBit(xAOD::EventInfo::Pixel,0x1)) { //FIXME an enum at some appropriate place to indicating 0x1 as
    ATH_MSG_WARNING(" cannot set flag bit for Pixel " );
  }
  return StatusCode::SUCCESS;
}


void PixelRodDecoder::checkUnequalNumberOfHeaders( const unsigned int nFragmentsPerFE[8], uint32_t robId ) const {
  unsigned int nFrags = 0;
  bool foundIssue = false;
  for (unsigned int i = 0; i < 8; ++i) {
    if (nFrags == 0) {
      if (nFragmentsPerFE[i] != 0) nFrags = nFragmentsPerFE[i];    // set nFrags on first non-zero occurence
    } 
    else {
      if (nFragmentsPerFE[i] != 0 && nFragmentsPerFE[i] != nFrags) {
	      foundIssue = true;
      }
    }
  }
  if ( not foundIssue ) return;

  // Got unequal number of headers per FE for same ROD, this means trouble.
  // Print value for each FE
  std::string errmsg;
  for (unsigned int j = 0; j < 8; ++j) {
    if (nFragmentsPerFE[j] != 0) {
      errmsg += "FE "+ std::to_string(j) + " " +std::to_string(nFragmentsPerFE[j]) + "  ";
    }
  }
  generalwarning("In ROB 0x" << std::hex << robId << ": got unequal number of headers per FE" << std::dec);
  generalwarning("[FE number] : [# headers] - " << errmsg);
}


unsigned int PixelRodDecoder::getLocalFEI4(const uint32_t fe, const uint64_t onlineId) const
{
  unsigned int linknum40 = (onlineId>>24) & 0xFF;
  unsigned int linknum80 = (onlineId>>32) & 0xFF;

  if (fe == linknum40) { 
    return 0;
  } else if (fe == linknum80) {
    return 1;
  } else {
    ATH_MSG_ERROR("Error in retrieving local FE-I4 number: linknumber " << fe << " not found in onlineID " << std::hex << onlineId);
  }
  return 0xF;
}
