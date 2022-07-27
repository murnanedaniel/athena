/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRT_GEOMODEL_TRT_DETDESCRDB_PARAMETERINTERFACE_H
#define TRT_GEOMODEL_TRT_DETDESCRDB_PARAMETERINTERFACE_H


#include "TRTParameterInterface.h"
#include "InDetGeoModelUtils/InDetDDAthenaComps.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"

class IRDBRecord;
class TopLevelPlacements;



class TRT_DetDescrDB_ParameterInterface : public TRTParameterInterface {

public:

   
  // Only allowed constructor
  TRT_DetDescrDB_ParameterInterface(InDetDD::AthenaComps * athenaComps);

  ~TRT_DetDescrDB_ParameterInterface();
  // delete copy c'tor
  TRT_DetDescrDB_ParameterInterface(const TRT_DetDescrDB_ParameterInterface &) = delete;
  // delete assignment
  TRT_DetDescrDB_ParameterInterface & operator=(const TRT_DetDescrDB_ParameterInterface &) = delete;
  //
  virtual const GeoTrf::Transform3D & partTransform(const std::string & partName) const;
  virtual bool partPresent(const std::string & partName) const;
  virtual InDetDD::DistortedMaterialManager * distortedMatManager();

  const std::string& name() const { 
    static const std::string n("TRT_GeoModel::TRT_DetDescrDB_ParameterInterface"); 
    return n;
  }

  IRDBRecordset_ptr  scalingTable() const {return m_scalingTable;}

  void SetValues();

  MsgStream& msg (MSG::Level lvl) { return m_athenaComps->msg(lvl); }

  
private:

  InDetDD::AthenaComps * m_athenaComps;

  InDetDD::DistortedMaterialManager * m_distortedMatManager;
  TopLevelPlacements * m_placements;
  IRDBRecordset_ptr m_scalingTable;
};

#endif
