/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

// std
#include <sys/stat.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

// other packages
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/MsgStream.h"
#include "Identifier/Identifier.h"
#include "MuonCablingData/MuonMDT_CablingMap.h"
#include "MuonIdHelpers/MdtIdHelper.h"
#include "PathResolver/PathResolver.h"
#include "StoreGate/DataHandle.h"
#include "StoreGate/StoreGateSvc.h"

// this package
#include "MdtCalibSvc/MdtCalibrationShiftMapBase.h"

#include "TTree.h"
// TRandom3 for random shifting
// should later be replaced with std
#include "TRandom3.h"

//
// private helper functions
//

MdtCalibrationShiftMapBase::MdtCalibrationShiftMapBase(const std::string &name,
                                                       ISvcLocator *pSvcLocator)
    : base_class(name, pSvcLocator),
      m_cablingSvc("MuonMDT_CablingSvc", name),
      m_mapIsInitialized(false),
      m_mapFileName(""),
      m_centralValue(0),
      m_sigma(10),
      m_forceMapRecreate(false) {
  declareProperty("MapFile", m_mapFileName);
  declareProperty("CentralValue", m_centralValue = 0);
  declareProperty("Sigma", m_sigma = 10);
  declareProperty("ForceMapRecreate", m_forceMapRecreate = false);
}

MdtCalibrationShiftMapBase::~MdtCalibrationShiftMapBase() { ; }

// return failure if not overloaded
StatusCode MdtCalibrationShiftMapBase::initializeMap() {
  return StatusCode::FAILURE;
}

StatusCode MdtCalibrationShiftMapBase::dumpMapAsFile() {
  /* initialize map if it's not there */
  if (!m_mapIsInitialized) {
    initializeMap();
  }

  /* write the map to a file */
  {
    std::ofstream file(m_mapFileName.c_str(), std::ios::binary);

    /* see if opening the file was successful */
    if (!file.is_open()) {
      ATH_MSG_FATAL(
          "Cannot open map output file for writing. Tried accessing file at "
          "\"./"
          << m_mapFileName.c_str() << "\"");
      return StatusCode::FAILURE;
    }

    /* dump map contents */
    for (auto shift : m_shiftValues) {
      file.write(reinterpret_cast<const char *>(&(shift.first)),
                 sizeof(Identifier));
      file.write(reinterpret_cast<const char *>(&(shift.second)),
                 sizeof(float));
    }
  }  // '}' flushes file
  return StatusCode::SUCCESS;
}

StatusCode MdtCalibrationShiftMapBase::loadMapFromFile() {
  /* check if map was already initialized */
  if (m_mapIsInitialized) {
    ATH_MSG_WARNING("Map already initialized, won't overwrite it.");
    return StatusCode::SUCCESS;
  }

  /* resolve path */
  std::string fileWithPath = PathResolver::find_file(m_mapFileName, "DATAPATH");
  if (fileWithPath == "") {
    ATH_MSG_ERROR("Cannot find file " << m_mapFileName << " in $DATAPATH");
    return StatusCode::FAILURE;
  }

  /* check if the file exists */
  struct stat buffer;
  if (stat(fileWithPath.c_str(), &buffer) != 0) {
    ATH_MSG_ERROR("Cannot stat the file \""
                  << fileWithPath.c_str()
                  << "\" -> map can not be initialized from this file.");
    return StatusCode::FAILURE;
  }

  /* check if the map is already stored as a binary file */
  std::ifstream fin(fileWithPath.c_str(), std::ios::binary);
  Identifier id;
  float shift;
  while (fin.read(reinterpret_cast<char *>(&id), sizeof(Identifier)) &&
         fin.read(reinterpret_cast<char *>(&shift), sizeof(float))) {
    m_shiftValues[id] = shift;
  }
  m_mapIsInitialized = true;
  ATH_MSG_INFO("Successfully initialized shift map from file \""
               << fileWithPath.c_str() << "\"");
  return StatusCode::SUCCESS;
}

float MdtCalibrationShiftMapBase::getValue(const Identifier &id) {
  if (!m_mapIsInitialized) {
    StatusCode sc = initializeMap();
    if (sc.isFailure()) {
      ATH_MSG_FATAL("Could not initialize the shift map.");
      throw GaudiException("Cannot run without shift map.", "", sc);
    }
  }

  return m_shiftValues[id];
}
