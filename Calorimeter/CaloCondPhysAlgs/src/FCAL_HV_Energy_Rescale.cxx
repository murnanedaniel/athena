///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

// FCAL_HV_Energy_Rescale.cxx 
// Implementation file for class FCAL_HV_Energy_Rescale
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 

// CaloCondPhysAlgs includes
#include "FCAL_HV_Energy_Rescale.h"
#include "CaloIdentifier/CaloCell_ID.h"
#include "StoreGate/ReadCondHandle.h"
#include "AthenaKernel/errorcheck.h"

#include "LArElecCalib/ILArHVScaleCorr.h"
#include "AthenaPoolUtilities/AthenaAttributeList.h"
#include "CaloCondBlobObjs/CaloCondBlobFlt.h"

#include "CoolKernel/types.h"
#include "CoolKernel/Record.h"
#include "CoralBase/AttributeListSpecification.h"

FCAL_HV_Energy_Rescale::FCAL_HV_Energy_Rescale( const std::string& name, 
						ISvcLocator* pSvcLocator ) : 
  AthAlgorithm( name, pSvcLocator )
{
  declareProperty("Folder",m_folder="/LAR/CellCorrOfl/EnergyCorr");
}


FCAL_HV_Energy_Rescale::~FCAL_HV_Energy_Rescale()
= default;


StatusCode FCAL_HV_Energy_Rescale::initialize()
{
  ATH_CHECK( m_cablingKey.initialize());
  ATH_CHECK( m_scaleCorrKey.initialize());
  return StatusCode::SUCCESS;
}

StatusCode FCAL_HV_Energy_Rescale::finalize()
{
  return StatusCode::SUCCESS;
}

StatusCode FCAL_HV_Energy_Rescale::execute()
{  
  return StatusCode::SUCCESS;
}

StatusCode FCAL_HV_Energy_Rescale::stop()
{  
  const EventContext& ctx = Gaudi::Hive::currentContext();
  SG::ReadCondHandle<ILArHVScaleCorr> scaleCorr (m_scaleCorrKey, ctx);

  const CaloCell_ID*    calocell_id;	
  CHECK(detStore()->retrieve(calocell_id,"CaloCell_ID"));
 
  IdentifierHash hashMin,hashMax;
  calocell_id->calo_cell_hash_range(CaloCell_ID::LARFCAL, hashMin,hashMax);

  SG::ReadCondHandle<LArOnOffIdMapping> cablingHdl{m_cablingKey};
  const LArOnOffIdMapping* cabling{*cablingHdl};
  if(!cabling) {
     ATH_MSG_ERROR("Do not have cabling mapping from key " << m_cablingKey.key() );
     return StatusCode::FAILURE;
  }

  ATH_MSG_INFO( "Working on hash range 0 to " << hashMax  );

  coral::AttributeListSpecification* spec = new coral::AttributeListSpecification();
  spec->extend("CaloCondBlob16M","blob");// , cool::StorageType::Blob16M);
  AthenaAttributeList* attrList=new AthenaAttributeList(*spec);
  coral::Blob& blob=(*attrList)["CaloCondBlob16M"].data<coral::Blob>();
  CaloCondBlobFlt* flt= CaloCondBlobFlt::getInstance(blob);
  spec->release(); // deletes spec
  // cppcheck-suppress memleak
  spec = nullptr;

  //Blob Definition Vector
  std::vector<std::vector<float> > defVec;
  defVec.emplace_back(1,1);
  flt->init(defVec,hashMax,1);

  CHECK(detStore()->record(attrList,m_folder));

  std::vector<float> setVec(1,1);
  unsigned nFCAL=0, nSmall=0, nSet=0;
  for(unsigned h=0;h<hashMax;++h) {
    float value=1.0;
    if (calocell_id->is_fcal(h)) {
      ++nFCAL;
      Identifier id=calocell_id->cell_id(h);
      HWIdentifier hwid=cabling->createSignalChannelID(id);
      const float corr=scaleCorr->HVScaleCorr(hwid);
      
      if (fabs(1.0-corr)>0.01) {
	value=corr;
	++nSet;
      }
      else
	++nSmall;
      ATH_MSG_INFO( "FCAL module " << calocell_id->sampling(id) << " corr=" << corr << " =>" << value  );
    }
    //std::cout << h << " " << value << std::endl;
    setVec[0]=value;
    flt->setData(h,0,setVec);
  }//end loop over hash

  ATH_MSG_INFO( "Found " << nFCAL << " FCAL channels of which " << nSet << " have a correction. (" << nSmall << " below threshold)"  );
  
  return StatusCode::SUCCESS;
}

