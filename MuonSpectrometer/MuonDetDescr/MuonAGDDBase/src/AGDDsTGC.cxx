/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/



#include "MuonAGDDBase/AGDDsTGC.h"
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

#include "MuonGeoModel/sTGCComponent.h"
#include "MuonGeoModel/sTGC.h"
#include "MuonGeoModel/MYSQL.h"

#include <sstream>

using MuonGM::MYSQL;


AGDDsTGC::AGDDsTGC(const std::string& s,
                   AGDDDetectorStore& ds,
                   AGDDVolumeStore& vs,
                   AGDDSectionStore& ss)
  : sTGCDetectorDescription(s,ds),AGDDVolume(s,vs,ss,true)
{
    Register();
}

void AGDDsTGC::CreateSolid (const AGDDBuilder& /*builder*/)
{

}

void AGDDsTGC::CreateVolume (AGDDBuilder& builder)
{
	
	MuonGM::sTGCComponent stgc_comp;
	stgc_comp.name=tech;
	stgc_comp.dx1=small_x();
	stgc_comp.dx2=large_x();
	stgc_comp.dy=y();
	stgc_comp.subType=subType();
	stgc_comp.yCutout=yCutout();
	stgc_comp.yCutoutCathode=yCutoutCathode();
	
	MuonGM::sTGC cham(&stgc_comp);
	GeoPhysVol *vvv=(GeoPhysVol*)cham.build(builder.GetMaterialManager(), 1);

	CreateSolid (builder);

	if (!GetVolume())
	{
		SetVolume(vvv);
	}
}

