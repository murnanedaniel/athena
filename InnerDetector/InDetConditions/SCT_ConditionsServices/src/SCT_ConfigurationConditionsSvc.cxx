/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_ConfigurationConditionsSvc.h"

// STL includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "boost/lexical_cast.hpp"

// Gaudi includes
#include "GaudiKernel/StatusCode.h"

// Athena includes
#include "StoreGate/StoreGateSvc.h"
#include "Identifier/IdentifierHash.h"

#include "InDetIdentifier/SCT_ID.h"
#include "InDetReadoutGeometry/SCT_DetectorManager.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"

//#include "AthenaKernel/IIOVDbSvc.h"

// Local includes
#include "SCT_ReadoutTool.h"
#include "SCT_Chip.h"

// Static folder names 
static const std::string coolChannelFolderName("/SCT/DAQ/Configuration/Chip");
static const std::string coolModuleFolderName("/SCT/DAQ/Configuration/Module");
static const std::string coolMurFolderName("/SCT/DAQ/Configuration/MUR");

// Static value to represent all 128 channels being good (signed int)
static const int noMask(-1);
// in case the chip number cannot be retrieved, this is the invalid value
static const int invalidChipNumber(-1);
// Constructor
SCT_ConfigurationConditionsSvc::SCT_ConfigurationConditionsSvc( const std::string& name, ISvcLocator* pSvcLocator ) : 
  AthService(name, pSvcLocator),
  m_badChannelIds(0),
  m_badModuleIds(0),
  m_badWaferIds(0),
  m_badLinks(0),
  m_badChips(0),
  m_filled(false),
  m_detStore("DetectorStore", name),
  m_IOVSvc("IOVSvc", name),
  m_IOVDbSvc("IOVDbSvc", name),
  m_pHelper(0),
  m_cablingSvc("SCT_CablingSvc", name),
  m_readoutTool("SCT_ReadoutTool"),
  m_pManager(0),
  m_checkStripsInsideModules(true) 
{ 
  declareProperty("checkStripsInsideModule" , m_checkStripsInsideModules);
}

// Initialize
StatusCode SCT_ConfigurationConditionsSvc::initialize(){
  msg(MSG:: INFO)<< "Initializing configuration" << endreq;

  // Retrieve cabling service
  if (m_cablingSvc.retrieve().isFailure())                return msg(MSG:: ERROR)<< "Can't get the cabling service." << endreq, StatusCode::FAILURE;

  // Retrieve detector store
  if (m_detStore.retrieve().isFailure())                  return msg(MSG:: FATAL)<< "Detector service  not found !" << endreq, StatusCode::FAILURE;

  // Retrieve SCT Detector Manager 
  if (m_detStore->retrieve(m_pManager,"SCT").isFailure()) return msg(MSG:: ERROR)<< "SCT mgr failed to retrieve" << endreq, StatusCode::FAILURE;
  
  // Retrieve SCT ID helper
  if (m_detStore->retrieve(m_pHelper, "SCT_ID").isFailure()) return msg(MSG::FATAL) << "Could not get SCT ID helper" << endreq, StatusCode::FAILURE;

  // Retrieve readout tools
  if (m_readoutTool.retrieve().isFailure())               return msg(MSG:: ERROR)<< "Could not retrieve SCT_ReadoutTool" << endreq, StatusCode::FAILURE;

  // Retrieve IOV service
  if (m_IOVSvc.retrieve().isFailure())                    return msg(MSG:: ERROR)<< "Failed to retrieve IOVSvc " << endreq, StatusCode::FAILURE;

  // Retrieve IOVDb service
  if (m_IOVDbSvc.retrieve().isFailure())                  return msg(MSG:: ERROR)<< "Failed to retrieve IOVDbSvc " << endreq, StatusCode::FAILURE;

  // Assign memory for structres
  m_badChannelIds =  new std::set<Identifier>;
  m_badModuleIds  =  new std::set<Identifier>;
  m_badWaferIds   =  new std::set<Identifier>;
  m_badLinks      =  new std::map<Identifier, std::pair<bool, bool> >;
  m_badChips      =  new std::map<Identifier, unsigned int>;

  // Register callbacks for folders 
  if (m_detStore->regFcn(&SCT_ConfigurationConditionsSvc::fillData,this, m_dataChannel,coolChannelFolderName).isFailure() or
      m_detStore->regFcn(&SCT_ConfigurationConditionsSvc::fillData, this, m_dataModule, coolModuleFolderName).isFailure() or
      m_detStore->regFcn(&SCT_ConfigurationConditionsSvc::fillData, this, m_dataMur, coolMurFolderName).isFailure())
    return msg(MSG:: ERROR)<< "Failed to register callback" << endreq, StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

// Finalize
StatusCode SCT_ConfigurationConditionsSvc::finalize(){
  msg(MSG:: INFO)<< "Configuration finalize" << endreq;

  if (m_badChannelIds) delete m_badChannelIds; 
  if (m_badModuleIds)  delete m_badModuleIds; 
  if (m_badWaferIds  ) delete m_badWaferIds;   
  if (m_badLinks     ) delete m_badLinks;      
  if (m_badChips     ) delete m_badChips;      

  return StatusCode::SUCCESS;
}

// Query interfaces.
StatusCode SCT_ConfigurationConditionsSvc::queryInterface(const InterfaceID& riid, void** ppvInterface) {

  if ( ISCT_ConfigurationConditionsSvc::interfaceID().versionMatch(riid) ) {
    *ppvInterface = this;
  } else if ( ISCT_ConditionsSvc::interfaceID().versionMatch(riid) ) {
    *ppvInterface = dynamic_cast<ISCT_ConditionsSvc*>(this);
  } else {
    // Interface is not directly available : try out a base class
    return AthService::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}

// What level of element can this service report about
bool SCT_ConfigurationConditionsSvc::canReportAbout(InDetConditions::Hierarchy h){
  return (h == InDetConditions::SCT_STRIP or h == InDetConditions::SCT_MODULE or 
	  h == InDetConditions::SCT_SIDE or h == InDetConditions::DEFAULT); 
}

// Is an element with this Identifier and hierachy good?
bool SCT_ConfigurationConditionsSvc::isGood(const Identifier & elementId, InDetConditions::Hierarchy  h){
  if (not canReportAbout(h)) return true;

  bool result(true);
  if (h == InDetConditions::SCT_STRIP) {
    result = (m_badChannelIds->find(elementId) == m_badChannelIds->end());
    // If strip itself is not bad, check if it's in a bad module
    if (result and m_checkStripsInsideModules) result = !isStripInBadModule(elementId);
  } else if (h == InDetConditions::SCT_MODULE) {
    result = (m_badModuleIds->find(elementId) == m_badModuleIds->end());
  } else if (h == InDetConditions::SCT_SIDE or h == InDetConditions::DEFAULT) {
    result = (m_badWaferIds->find(elementId) == m_badWaferIds->end());
  }
  return result;
}

// Is a wafer with this IdentifierHash good?
bool SCT_ConfigurationConditionsSvc::isGood(const IdentifierHash & hashId){
  Identifier elementId(m_pHelper->wafer_id(hashId));
  return isGood(elementId);
}

  // Callback funtion (with arguments required by IovDbSvc) to fill channel, module and link info
StatusCode SCT_ConfigurationConditionsSvc::fillData(int& /*i*/ , std::list<std::string>& l){  
  std::list<std::string>::iterator itr(l.begin());
  std::list<std::string>::iterator itrEnd(l.end());

  // Fill module data if Module folder has changed
  if (find(itr, itrEnd, coolModuleFolderName) != itrEnd) {
    if (fillModuleData().isFailure()) return StatusCode::FAILURE;
  }

  // Fill strip, chip and link info if Chip or MUR folders change
  if (find(itr, itrEnd, coolChannelFolderName) != itrEnd or
      find(itr, itrEnd, coolMurFolderName)     != itrEnd) {
    if(fillChannelData().isFailure()) return StatusCode::FAILURE;
  }

  // The bad channel list contains all the information
  m_filled = (m_badChannelIds->size() != 0);

  return StatusCode::SUCCESS;
}

// Fill bad strip, chip and link info
StatusCode SCT_ConfigurationConditionsSvc::fillChannelData(){

  // Clear previous information at callback
  m_badChannelIds->clear();
  m_badChips->clear();

  // Fill link status
  if (fillLinkStatus().isFailure()) return StatusCode::FAILURE;

  // Get Chip folder
  if (retrieveFolder(m_dataChannel, coolChannelFolderName).isFailure()) {
    return msg(MSG:: ERROR)<< "Could not fill channel configuration data" << endreq, StatusCode::FAILURE;
  } else {
    msg(MSG:: INFO)<< "fillChannelData: IOV callback resulted in a Chip CondAttrListVec of size " << m_dataChannel->size() << endreq;
  }

  // Loop over modules (i.e groups of 12 chips) in DB folder 
  const unsigned int nChips(12);
  CondAttrListVec::const_iterator modItr(m_dataChannel->begin());
  CondAttrListVec::const_iterator modEnd(m_dataChannel->end());

  for (;modItr != modEnd; modItr += nChips) {

    // Get SN and identifiers (the channel number is serial number+1)
    const unsigned int truncatedSerialNumber(modItr->first - 1);
    const IdentifierHash& hash = m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber);
    if (not hash.is_valid()) continue;
    Identifier  waferId(m_pHelper->wafer_id(hash));
    Identifier  moduleId(m_pHelper->module_id(waferId));

    // Don't need to bother checking chips if the module is already bad
    // Commented out until fully tested
    //if (m_badModuleIds->find(moduleId) == m_badModuleIds->end()) continue;

    // Get link status 
    // Can maybe be smarter if both links are bad (but the module will probably be bad then)
    bool link0ok(true), link1ok(true);
    std::map<Identifier, std::pair<bool, bool> >::const_iterator linkItr =  m_badLinks->find(moduleId);
    if (linkItr != m_badLinks->end()) {
      link0ok = (*linkItr).second.first;
      link1ok = (*linkItr).second.second;
    }

    // Loop over chips within module
    CondAttrListVec::const_iterator channelItr(modItr);
    CondAttrListVec::const_iterator channelEnd(modItr + nChips);
    std::vector<SCT_Chip*> chipsInMod;
    chipsInMod.reserve(12);

    for(; channelItr != channelEnd; ++channelItr){
      // Get chip id, config and masks and store as SCT_Chip object
      // Can get AttributeList from second (see http://lcgapp.cern.ch/doxygen/CORAL/CORAL_1_9_3/doxygen/html/classcoral_1_1_attribute_list.html)
      short id      = channelItr->second[2].data<short>();
      short config  = channelItr->second[5].data<short>();
      int mask0     = channelItr->second[6].data<int>();  // These look like ints in the DB, I think they should be unsigned int
      int mask1     = channelItr->second[7].data<int>();  // since a mask of 0xFFFFFFFF (= none masked) shows as -1 
      int mask2     = channelItr->second[8].data<int>();  // (=noMask, declared as static int at top of this file)
      int mask3     = channelItr->second[9].data<int>();
      chipsInMod.push_back(new SCT_Chip(id, config, mask0, mask1, mask2, mask3));
    }

    // Check the module readout to look for bypassed chips, disabled links etc
    if (m_readoutTool->determineReadout(moduleId, chipsInMod, link0ok, link1ok).isFailure()) return StatusCode::FAILURE; 

    // Loop over chips again now know whether they're in the readout
    std::vector<SCT_Chip*>::const_iterator chipItr(chipsInMod.begin());
    std::vector<SCT_Chip*>::const_iterator chipEnd(chipsInMod.end());
    std::vector<int> badStripsVec;

    unsigned int chipStatusWord(0);
    for(; chipItr != chipEnd; ++chipItr){
      
      // Bad strips (only need to do this if at least one bad channel)
      if ((*chipItr)->numberOfMaskedChannels() != 0){
	// Add bad stips to vector
        badStripsVec.clear();
        (*chipItr)->appendBadStripsToVector(badStripsVec);

	// Loop over bad stips and insert strip ID into set
        std::vector<int>::const_iterator stripItr(badStripsVec.begin());
        std::vector<int>::const_iterator stripEnd(badStripsVec.end());

        for(;stripItr != stripEnd; ++stripItr){
          Identifier stripId(getStripId(truncatedSerialNumber, (*chipItr)->id(), *stripItr));
	  // If in rough order, may be better to call with itr of previous insertion as a suggestion	  
          if (stripId.is_valid()) m_badChannelIds->insert(stripId);
	}
      }

      // Bad chips (= all strips bad) bitpacked
      // Should only do this for modules with at least one chip bad?
      if ((*chipItr)->numberOfMaskedChannels() == stripsPerChip) chipStatusWord |= (1<<(*chipItr)->id());
    }

    // Store chip status if not all good (==0)
    if (chipStatusWord != 0) (*m_badChips)[moduleId] = chipStatusWord;

    // Clear up memory associated with chips    
    for (chipItr = chipsInMod.begin(); chipItr != chipEnd; ++chipItr){
      delete *chipItr;
    }
  }

  // No longer need the conditions folder as stored locally
  m_IOVDbSvc->dropObject(coolChannelFolderName,true); 

  const unsigned int totalBad(m_badChannelIds->size());
  msg(MSG:: INFO)<< "Total number of bad strip identifiers is " << totalBad << endreq;
  return StatusCode::SUCCESS;
}

// Fill bad module info
StatusCode SCT_ConfigurationConditionsSvc::fillModuleData(){
  unsigned int totalNumberOfModules(0);
  unsigned int totalNumberOfValidModules(0);

  // Clear previous information at callback
  m_badModuleIds->clear();
  m_badWaferIds->clear();

  // Get Module folder
  if (retrieveFolder(m_dataModule, coolModuleFolderName).isFailure()) {
    return msg(MSG:: ERROR)<< "Could not fill module configuration data" << endreq, StatusCode::FAILURE;
  } else {
    msg(MSG:: INFO)<< "fillModuleData: IOV callback resulted in a CondAttrListVec of size " << m_dataModule->size() << endreq;
  }

  // Loop over modules in DB folder
  CondAttrListVec::const_iterator pModule(m_dataModule->begin());
  CondAttrListVec::const_iterator pLastModule(m_dataModule->end());

  for (;pModule != pLastModule; ++pModule){
    // Get SN and Identifiers (the module number is serial number+1)
    const unsigned int truncatedSerialNumber(pModule->first - 1);
    const IdentifierHash &hash=m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber);
    ++totalNumberOfModules;
    if (not hash.is_valid()) continue;
    Identifier  waferId(m_pHelper->wafer_id(hash));
    ++totalNumberOfValidModules;
    IdentifierHash oppWaferHash;
    m_pHelper->get_other_side(m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber) , oppWaferHash);
    Identifier     oppWaferId(m_pHelper->wafer_id(oppWaferHash));
    Identifier     moduleId(m_pHelper->module_id(waferId));

    // Get AttributeList from second (see http://lcgapp.cern.ch/doxygen/CORAL/CORAL_1_9_3/doxygen/html/classcoral_1_1_attribute_list.html)
    // and get module info from this.  Bad module has a -ve group .
    short group(pModule->second[3].data<short>());

    if (group < 0) { 
      // Insert module/wafer ID into set of bad modules/wafers IDs
      m_badModuleIds->insert(moduleId);
      m_badWaferIds->insert(waferId);
      m_badWaferIds->insert(oppWaferId);
    }
  }  

  // No longer need the conditions folder as stored locally
  m_IOVDbSvc->dropObject(coolModuleFolderName,true); 

  const unsigned int totalBad(m_badModuleIds->size());
  msg(MSG:: INFO)<< "Total number of module identifiers is " << totalNumberOfModules << endreq;
  msg(MSG:: INFO)<< "Total number of modules also found in the cabling is " << totalNumberOfValidModules << endreq;
  msg(MSG:: INFO)<< "Total number of bad module identifiers is " << totalBad << endreq;    
  return StatusCode::SUCCESS;  
}

// Is the information filled?
bool SCT_ConfigurationConditionsSvc::filled() const{
  return m_filled;
}

// Get a DB folder
StatusCode SCT_ConfigurationConditionsSvc::retrieveFolder(const DataHandle<CondAttrListVec> &pDataVec, const std::string & folderName){
  if (not m_detStore) return (msg(MSG:: FATAL) << "The detector store pointer is NULL" << endreq), StatusCode::FAILURE;

  if (m_detStore->retrieve(pDataVec, folderName).isFailure()) 
    return (msg(MSG:: FATAL) << "Could not retrieve AttrListVec for " << folderName << endreq), StatusCode::FAILURE;

  if (0 == pDataVec->size()) return (msg(MSG:: FATAL) << "This folder's data set appears to be empty: " << folderName << endreq), StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

// Construct the strip ID from the module SN, chip number and strip number
Identifier SCT_ConfigurationConditionsSvc::getStripId(const unsigned int truncatedSerialNumber, const unsigned int chipNumber, const unsigned int stripNumber) const{

  Identifier     waferId;
  unsigned int   strip(0);
  IdentifierHash waferHash;

  if (not m_cablingSvc) msg(MSG:: FATAL)<< "The cabling tool pointer is zero." << endreq;

  // If the chip is 0-5 we are in side 0, otherwise in side 1. 'getHash' only 
  // returns the side 0 hash, so we use the helper for side 1

  if (chipNumber<6){
    waferHash = m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber);
    strip     = chipNumber * stripsPerChip + stripNumber;
    if (waferHash.is_valid()) waferId   = m_pHelper->wafer_id(waferHash);
  } else {
    m_pHelper->get_other_side( m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber) , waferHash);
    strip   = (chipNumber - 6) * stripsPerChip + stripNumber;
    if (waferHash.is_valid()) waferId = m_pHelper->wafer_id(waferHash);
  }

  if (not waferId.is_valid()) return waferId;

  const InDetDD::SiDetectorElement* pElement = (m_pManager->getDetectorElement(waferHash));
  if (! pElement) msg(MSG:: FATAL)<< "Element pointer is NULL in 'getStripId' method" << endreq;
  strip = (pElement->swapPhiReadoutDirection()) ? lastStrip - strip: strip;

  return m_pHelper->strip_id(waferId, strip);
}

// Fill link info
StatusCode SCT_ConfigurationConditionsSvc::fillLinkStatus() {
  // Clear previous information at call back
  m_badLinks->clear();

  // Get MUR folders for link info 
  if (retrieveFolder(m_dataMur, coolMurFolderName).isFailure()) {
    return msg(MSG:: ERROR)<< "Could not fill MUR configuration data" << endreq, StatusCode::FAILURE;
  } else {
    msg(MSG:: INFO)<< "fillLinkStatus: IOV callback resulted in a MUR CondAttrListColl of size " << m_dataMur->size() << endreq;
  }

  // loop over MUR folder
  CondAttrListVec::const_iterator pMur(m_dataMur->begin());
  CondAttrListVec::const_iterator pLastMur(m_dataMur->end());

  for (; pMur != pLastMur; ++pMur) {
    // Check for null values
    if (pMur->second[4].isNull()) continue;
    long long serialNumber     = pMur->second[4].data<long long>();
    int truncatedSerialNumber  = truncateSerialNumber(serialNumber);    
    const IdentifierHash& hash = m_cablingSvc->getHashFromSerialNumber(truncatedSerialNumber);
    if (not hash.is_valid()) continue;

    Identifier  waferId(m_pHelper->wafer_id(hash));
    Identifier  moduleId(m_pHelper->module_id(waferId));

    int link0 = pMur->second[6].data<int>();
    int link1 = pMur->second[7].data<int>();

    // Store the modules with bad links, represented by badLink (enum in header) = 255 = 0xFF 
    if (link0 == badLink or link1 == badLink) {
      (*m_badLinks)[moduleId] = std::make_pair((link0!=badLink), (link1!=badLink));
    }

  }

  // No longer need the conditions folder as stored locally
  m_IOVDbSvc->dropObject(coolMurFolderName,true); 
  return StatusCode::SUCCESS;
}

// Truncate a serial number
int SCT_ConfigurationConditionsSvc::truncateSerialNumber(long long serialNumber) {
  std::string snString = boost::lexical_cast<std::string>(serialNumber);
  return boost::lexical_cast<int>(snString.substr(5));
}

// Check if a strip is within a bad module
bool SCT_ConfigurationConditionsSvc::isStripInBadModule(const Identifier& stripId){
  Identifier moduleId(m_pHelper->module_id(m_pHelper->wafer_id(stripId)));
  return (m_badModuleIds->find(moduleId) != m_badModuleIds->end());
}

// Check if a wafer is within a bad module
bool SCT_ConfigurationConditionsSvc::isWaferInBadModule(const Identifier& waferId){
  Identifier moduleId(m_pHelper->module_id(waferId));
  return (m_badModuleIds->find(moduleId) != m_badModuleIds->end());
}

// Find the chip number containing a particular strip Identifier
int SCT_ConfigurationConditionsSvc::getChip(Identifier stripId) {

  // Find side and strip number
  int side(m_pHelper->side(stripId));
  int strip(m_pHelper->strip(stripId));

  // Check for swapped readout direction
  IdentifierHash waferHash = m_pHelper->wafer_hash(m_pHelper->wafer_id(stripId));
  const InDetDD::SiDetectorElement* pElement = m_pManager->getDetectorElement(waferHash);
  if (! pElement){
     msg(MSG:: FATAL)<< "Element pointer is NULL in 'badStrips' method" << endreq;
     return invalidChipNumber;
  }
  strip = (pElement->swapPhiReadoutDirection()) ? lastStrip - strip: strip;

  // Find chip number
  return (side==0 ? strip/stripsPerChip : strip/stripsPerChip + 6);
}

void SCT_ConfigurationConditionsSvc::badStrips(Identifier moduleId,  std::set<Identifier>& strips, bool ignoreBadModules, bool ignoreBadChips) {
  // Bad strips for a given module

  if (ignoreBadModules) {
    // Ignore strips in bad modules
    if (m_badModuleIds->find(moduleId) != m_badModuleIds->end()) return;    
  }

  std::set<Identifier>::const_iterator chanItr(m_badChannelIds->begin());
  std::set<Identifier>::const_iterator chanEnd(m_badChannelIds->end());
  for (; chanItr != chanEnd; ++chanItr) {
    if (ignoreBadChips) {
      // Ignore strips in bad chips
      int chip = getChip(*chanItr);
      if (chip!=invalidChipNumber){
        std::map<Identifier, unsigned int>::const_iterator chipItr(m_badChips->find(moduleId));
        if ((chipItr != m_badChips->end()) && ((*chipItr).second & (1 << chip)) != 0) continue;
      }   
    }
    if (m_pHelper->module_id(m_pHelper->wafer_id((*chanItr))) == moduleId) strips.insert(*chanItr);
  }
}
       
std::pair<bool, bool> SCT_ConfigurationConditionsSvc::badLinks(Identifier moduleId) {
  // Bad links for a given module
  std::map<Identifier, std::pair<bool, bool> >::const_iterator linkItr(m_badLinks->find(moduleId));
  return ((linkItr != m_badLinks->end()) ? (*linkItr).second : std::make_pair(true,true));
}

unsigned int SCT_ConfigurationConditionsSvc::badChips(Identifier moduleId) {
  // Bad chips for a given module
  std::map<Identifier, unsigned int>::const_iterator chipItr(m_badChips->find(moduleId));  
  return ((chipItr != m_badChips->end()) ? (*chipItr).second : static_cast<unsigned int>(0));
}

void SCT_ConfigurationConditionsSvc::badStrips(std::set<Identifier>& strips, bool ignoreBadModules, bool ignoreBadChips) {
  if (ignoreBadModules == false && ignoreBadChips == false) {
    //return m_badChannelIds;
    std::copy(m_badChannelIds->begin(), m_badChannelIds->end(), std::inserter(strips,strips.begin()));
    return;
  }

  // How deal with memory management - THIS NEEDS SORTING OUT (can't have user repsonsible in one case and service in another)
  //std::set<Identifier>* strips = new std::set<Identifier>;

  std::set<Identifier>::const_iterator chanItr(m_badChannelIds->begin());
  std::set<Identifier>::const_iterator chanEnd(m_badChannelIds->end());
  for (; chanItr != chanEnd; ++chanItr) {
    Identifier moduleId(m_pHelper->module_id(m_pHelper->wafer_id(*chanItr)));
    // Ignore strips in bad modules
    if (ignoreBadModules) {
      if (m_badModuleIds->find(moduleId) != m_badModuleIds->end()) continue;    
    }
    // Ignore strips in bad chips
    if (ignoreBadChips) {    
      int chip = getChip(*chanItr);
      if (chip !=invalidChipNumber){
        std::map<Identifier, unsigned int>::const_iterator chipItr(m_badChips->find(moduleId));
        if ((chipItr != m_badChips->end()) && ((*chipItr).second & (1 << chip)) != 0) continue;
      }  
    }

    strips.insert(*chanItr);
  }

  //return strips;
  //return m_badChannelIds;
}
