/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetEventTPCnv/InDetRIO_OnTrack/SCT_ClusterOnTrackCnv_p2.h"

#include "AthLinks/ElementLink.h"
#include "Identifier/Identifier.h" 
#include "InDetPrepRawData/SCT_ClusterContainer.h"
#include "InDetRIO_OnTrack/SCT_ClusterOnTrack.h"
#include "TrkEventTPCnv/helpers/EigenHelpers.h"

void SCT_ClusterOnTrackCnv_p2::persToTrans( const InDet::SCT_ClusterOnTrack_p2 *persObj, InDet::SCT_ClusterOnTrack *transObj, MsgStream &log ) {

  ElementLinkToIDCSCT_ClusterContainer rio;
  m_elCnv.persToTrans(&persObj->m_prdLink,&rio,log);  

  Trk::LocalParameters localParams;
  fillTransFromPStore( &m_localParCnv, persObj->m_localParams, &localParams, log );

  Trk::ErrorMatrix dummy;
  Amg::MatrixX localCovariance;
  fillTransFromPStore( &m_errorMxCnv, persObj->m_localErrMat, &dummy, log );
  EigenHelpers::vectorToEigenMatrix(dummy.values, localCovariance, "SCT_ClusterOnTrackCnv_p2");
    
  // When reading in 32-bit id, must cast to unsigned int
  const Identifier::value_type upper = 0XFFFFFFFF00000000LL;
  const Identifier::value_type lower = 0X00000000FFFFFFFFLL;
  const Identifier::value_type testUpper = persObj->m_id & upper;
  const Identifier::value_type testLower = persObj->m_id & lower;
  Identifier id;
  if ( testUpper == 0 && testLower > 0) {
    Identifier32::value_type id1 = persObj->m_id;
    id = id1;
  } else {
    id = persObj->m_id;
  }

  *transObj = InDet::SCT_ClusterOnTrack (rio,
                                         localParams,
                                         localCovariance,
                                         persObj->m_idDE,
                                         id,
                                         persObj->m_isbroad,
                                         persObj->m_positionAlongStrip);

  // Attempt to call supertool to fill in detElements
  m_eventCnvTool->recreateRIO_OnTrack(transObj);
  if (transObj->detectorElement()==nullptr) {
    log << MSG::WARNING<<"Unable to reset DetEl for this RIO_OnTrack, "
        << "probably because of a problem with the Identifier/IdentifierHash : ("
        << transObj->identify()<<"/"<<transObj->idDE()<<endmsg;
  }
}

void SCT_ClusterOnTrackCnv_p2::transToPers(const InDet::SCT_ClusterOnTrack* transObj, InDet::SCT_ClusterOnTrack_p2* persObj, MsgStream& log) {  
  if (transObj==nullptr or persObj==nullptr) return;

  persObj->m_id = transObj->identify().get_compact();
  persObj->m_localParams = toPersistent( &m_localParCnv, &transObj->localParameters(), log );
  Trk::ErrorMatrix pMat;
  EigenHelpers::eigenMatrixToVector(pMat.values, transObj->localCovariance(), "SCT_ClusterOnTrackCnv_p2");
  persObj->m_localErrMat = toPersistent( &m_errorMxCnv, &pMat, log );
  persObj->m_idDE = transObj->idDE();
  persObj->m_isbroad = transObj->isBroadCluster();
  persObj->m_positionAlongStrip = static_cast<float>(transObj->positionAlongStrip());

  static const SG::InitializedReadHandleKey<InDet::SCT_ClusterContainer> sctClusContName ("SCT_Clusters");
  ElementLink<InDet::SCT_ClusterContainer>::index_type hashAndIndex{0};
  bool isFound{m_eventCnvTool->getHashAndIndex<InDet::SCT_ClusterContainer, InDet::SCT_ClusterOnTrack>(transObj, sctClusContName, hashAndIndex)};
  if(m_eventCnvTool->doTrackOverlay()){
    persObj->m_prdLink.m_contName = (isFound ? "Bkg_SCT_Clusters" : "");
    if(!isFound){ //in this case the input collection is called Bkg_SCT_Clusters as well
      static const SG::InitializedReadHandleKey<InDet::SCT_ClusterContainer> sctClusContName("Bkg_SCT_Clusters");
      isFound=m_eventCnvTool->getHashAndIndex<InDet::SCT_ClusterContainer, InDet::SCT_ClusterOnTrack>(transObj, sctClusContName, hashAndIndex);
      persObj->m_prdLink.m_contName = (isFound ? "Bkg_SCT_Clusters" : "");
    }
  }
  else persObj->m_prdLink.m_contName = (isFound ? sctClusContName.key() : "");
  persObj->m_prdLink.m_elementIndex = hashAndIndex;
}
