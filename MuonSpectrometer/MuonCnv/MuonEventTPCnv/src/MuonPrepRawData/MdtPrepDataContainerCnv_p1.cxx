/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonPrepRawData/MdtPrepData.h"
#include "MuonPrepRawData/MdtPrepDataContainer.h"
#include "MuonEventTPCnv/MuonPrepRawData/MdtPrepData_p1.h"
#include "MuonEventTPCnv/MuonPrepRawData/MuonPRD_Container_p1.h"

#include "MuonIdHelpers/MdtIdHelper.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonEventTPCnv/MuonPrepRawData/MdtPrepDataCnv_p1.h"
#include "MuonEventTPCnv/MuonPrepRawData/MdtPrepDataContainerCnv_p1.h"

// Gaudi
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IIncidentSvc.h"

// Athena
#include "StoreGate/StoreGateSvc.h"



#include "DataModel/DataPool.h"

StatusCode Muon::MdtPrepDataContainerCnv_p1::initialize(MsgStream &log) {
   // Do not initialize again:
    m_isInitialized=true;

   // Get Storegate, ID helpers, and so on
    ISvcLocator* svcLocator = Gaudi::svcLocator();
   // get StoreGate service
    StatusCode sc = svcLocator->service("StoreGateSvc", m_storeGate);
    if (sc.isFailure()) {
        log << MSG::FATAL << "StoreGate service not found !" << endreq;
        return StatusCode::FAILURE;
    }

   // get DetectorStore service
    StoreGateSvc *detStore;
    sc = svcLocator->service("DetectorStore", detStore);
    if (sc.isFailure()) {
        log << MSG::FATAL << "DetectorStore service not found !" << endreq;
        return StatusCode::FAILURE;
    } else {
        if (log.level() <= MSG::DEBUG) log << MSG::DEBUG << "Found DetectorStore." << endreq;
    }

   // Get the pixel helper from the detector store
    sc = detStore->retrieve(m_MdtId);
    if (sc.isFailure()) {
        log << MSG::FATAL << "Could not get Mdt ID helper !" << endreq;
        return StatusCode::FAILURE;
    } else {
        if (log.level() <= MSG::DEBUG) log << MSG::DEBUG << "Found the Mdt ID helper." << endreq;
    }

    sc = detStore->retrieve(m_muonDetMgr);
    if (sc.isFailure()) {
        log << MSG::FATAL << "Could not get PixelDetectorDescription" << endreq;
        return sc;
    }

    if (log.level() <= MSG::DEBUG) log << MSG::DEBUG << "Converter initialized." << endreq;
    return StatusCode::SUCCESS;
}

void Muon::MdtPrepDataContainerCnv_p1::transToPers(const Muon::MdtPrepDataContainer* transCont,  Muon::MuonPRD_Container_p1* persCont, MsgStream &log) 
{

    // The transient model has a container holding collections and the
    // collections hold channels.
    //
    // The persistent model flattens this so that the persistent
    // container has two vectors:
    //   1) all collections, and
    //   2) all PRDs
    //
    // The persistent collections, then only maintain indexes into the
    // container's vector of all channels. 
    //
    // So here we loop over all collection and add their channels
    // to the container's vector, saving the indexes in the
    // collection. 

    typedef Muon::MdtPrepDataContainer TRANS;
    typedef ITPConverterFor<Trk::PrepRawData> CONV;

    MdtPrepDataCnv_p1  chanCnv;
    TRANS::const_iterator it_Coll     = transCont->begin();
    TRANS::const_iterator it_CollEnd  = transCont->end();
    unsigned int collIndex;
    unsigned int chanBegin = 0;
    unsigned int chanEnd = 0;
    int numColl = transCont->numberOfCollections();
    // if(numColl == transCont->hashFunc().max() ) { // let's count how many collections we have:
    //  numColl = 0;
    //  for ( ; it_Coll != it_CollEnd; it_Coll++)
    //     numColl++;
    //  it_Coll     = transCont->begin(); // reset the iterator, we used it!
    // }
    persCont->m_collections.resize(numColl);
    if (log.level() <= MSG::DEBUG) log << MSG::DEBUG  << " Preparing " << persCont->m_collections.size() << "Collections" << endreq;

    for (collIndex = 0; it_Coll != it_CollEnd; ++collIndex, it_Coll++)  {
        // Add in new collection
        if (log.level() <= MSG::DEBUG) log << MSG::DEBUG  << " New collection" << endreq;
        const Muon::MdtPrepDataCollection& collection = (**it_Coll);
        chanBegin  = chanEnd;
        chanEnd   += collection.size();
        Muon::MuonPRD_Collection_p1& pcollection = persCont->m_collections[collIndex];
        pcollection.m_id    = collection.identify().get_identifier32().get_compact();
        pcollection.m_hashId = (unsigned int) collection.identifyHash();
        pcollection.m_begin = chanBegin;
        pcollection.m_end   = chanEnd;
        // Add in channels
        persCont->m_PRD.resize(chanEnd);
        for (unsigned int i = 0; i < collection.size(); ++i) {
            const Muon::MdtPrepData* chan = collection[i];
            persCont->m_PRD[i + chanBegin] = toPersistent((CONV**)0, chan, log );
        }
    }
    if (log.level() <= MSG::DEBUG) log << MSG::DEBUG  << " ***  Writing MdtPrepDataContainer ***" << endreq;
}

void  Muon::MdtPrepDataContainerCnv_p1::persToTrans(const Muon::MuonPRD_Container_p1* persCont, Muon::MdtPrepDataContainer* transCont, MsgStream &log) 
{

    // The transient model has a container holding collections and the
    // collections hold channels.
    //
    // The persistent model flattens this so that the persistent
    // container has two vectors:
    //   1) all collections, and
    //   2) all channels
    //
    // The persistent collections, then only maintain indexes into the
    // container's vector of all channels. 
    //
    // So here we loop over all collection and extract their channels
    // from the vector.


    Muon::MdtPrepDataCollection* coll = 0;

    MdtPrepDataCnv_p1  chanCnv;
    typedef ITPConverterFor<Trk::PrepRawData> CONV;

    if (log.level() <= MSG::DEBUG) log << MSG::DEBUG  << " Reading " << persCont->m_collections.size() << "Collections" << endreq;
    for (unsigned int icoll = 0; icoll < persCont->m_collections.size(); ++icoll) {

        // Create trans collection - is NOT owner of MdtPrepData (SG::VIEW_ELEMENTS)
    // IDet collection don't have the Ownership policy c'tor
        const Muon::MuonPRD_Collection_p1& pcoll = persCont->m_collections[icoll];        
        Identifier collID(Identifier(pcoll.m_id));
        IdentifierHash collIDHash(IdentifierHash(pcoll.m_hashId));
        coll = new Muon::MdtPrepDataCollection(collIDHash);
        coll->setIdentifier(Identifier(pcoll.m_id));
        unsigned int nchans           = pcoll.m_end - pcoll.m_begin;
        coll->resize(nchans);
//        MuonDD::MdtReadoutElement * de = m_muonDetMgr->getDetectorElement(collIDHash);
// No hash based lookup for Mdts?
        // Fill with channels
        for (unsigned int ichan = 0; ichan < nchans; ++ ichan) {
            const TPObjRef pchan = persCont->m_PRD[ichan + pcoll.m_begin];
            Muon::MdtPrepData* chan = dynamic_cast<Muon::MdtPrepData*>(createTransFromPStore((CONV**)0, pchan, log ) );
            if (chan!=0) {
                const MuonGM::MdtReadoutElement * de = m_muonDetMgr->getMdtReadoutElement(chan->identify());
                chan->m_detEl = de;
                (*coll)[ichan] = chan;
            } else {
                log<<MSG::WARNING<<"Dynamic cast or createTransFromPStore failed - channel missing"<<endreq;
            }
        }

        // register the rdo collection in IDC with hash - faster addCollection
        StatusCode sc = transCont->addCollection(coll, collIDHash);
        if (sc.isFailure()) {
            throw std::runtime_error("Failed to add collection to ID Container");
        }
        if (log.level() <= MSG::DEBUG) {
            log << MSG::DEBUG << "AthenaPoolTPCnvIDCont::persToTrans, collection, hash_id/coll id = " << (int) collIDHash << " / " << 
                collID.get_compact() << ", added to Identifiable container." << endreq;
        }
    }

    if (log.level() <= MSG::DEBUG) log << MSG::DEBUG  << " ***  Reading MdtPrepDataContainer" << endreq;
}



//================================================================
Muon::MdtPrepDataContainer* Muon::MdtPrepDataContainerCnv_p1::createTransient(const Muon::MuonPRD_Container_p1* persObj, MsgStream& log) 
{
    if(!m_isInitialized) {
        if (this->initialize(log) != StatusCode::SUCCESS) {
            log << MSG::FATAL << "Could not initialize MdtPrepDataContainerCnv_p1 " << endreq;
            return 0;
        } 
    }
    std::auto_ptr<Muon::MdtPrepDataContainer> trans(new Muon::MdtPrepDataContainer(m_MdtId->module_hash_max()));
    persToTrans(persObj, trans.get(), log);
    return(trans.release());
}


