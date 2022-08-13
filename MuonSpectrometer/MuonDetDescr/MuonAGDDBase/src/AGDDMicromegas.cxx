/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonAGDDBase/AGDDMicromegas.h"
#include "AGDDModel/AGDDParameterStore.h"
#include "AGDDKernel/AGDDDetectorStore.h"
#include "AGDDKernel/AGDDVolume.h"
#include "AGDDKernel/AGDDBuilder.h"

#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoMaterial.h"

#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "GeoModelInterfaces/StoredMaterialManager.h"

#include "MuonGeoModel/MicromegasComponent.h"
#include "MuonGeoModel/Micromegas.h"
#include "MuonGeoModel/MYSQL.h"


using MuonGM::MYSQL;

AGDDMicromegas::AGDDMicromegas(const std::string& s,
                               AGDDDetectorStore& ds,
                               AGDDVolumeStore& vs,
                               AGDDSectionStore& ss)
  : MMDetectorDescription(s,ds),AGDDVolume(s,vs,ss,true)
{
	Register();
}

void AGDDMicromegas::CreateSolid (const AGDDBuilder& /*builder*/)
{
//	std::cout<<"this is AGDDMicromegas::CreateSolid()"<<std::endl;
//	void *p=GetSolid();
//	if (!p)
//	{
//		std::cout<<" creating solid with dimensions "<<
//		m_small_x<<" "<<m_large_x<<" "<<m_y<<" "<<m_z<<std::endl;
//		GeoShape* solid=new GeoTrd(m_small_x/2.,m_large_x/2.,m_y/2.,m_y/2.,m_z/2.);
//		SetSolid(solid);
//	}

}

void AGDDMicromegas::CreateVolume (AGDDBuilder& builder)
{
//    std::cout<<"this is AGDDMicromegas::CreateVolume()"<<std::endl;
	
	MuonGM::MicromegasComponent mm_comp;
	mm_comp.name=tech;
	mm_comp.dx1=small_x();
	mm_comp.dx2=large_x();
	mm_comp.dy=y();
	mm_comp.subType=subType();
	
	MuonGM::Micromegas cham (&mm_comp);
	GeoPhysVol *vvv=(GeoPhysVol*)cham.build(builder.GetMaterialManager(), 1);

	CreateSolid (builder);

	if (!GetVolume())
	{
		SetVolume(vvv);
	}
}

