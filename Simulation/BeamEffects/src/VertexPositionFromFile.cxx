/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// class header include
#include "VertexPositionFromFile.h"

// HepMC includes
#include "AtlasHepMC/GenEvent.h"
// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
// Athena headers
#include "StoreGate/ReadHandle.h"

namespace Simulation
{

  /** Constructor */
  VertexPositionFromFile::VertexPositionFromFile( const std::string& t,
                                                  const std::string& n,
                                                  const IInterface* p )
    : base_class(t,n,p)
  {
  }

  /** Athena algtool's Hooks */
  StatusCode VertexPositionFromFile::initialize()
  {
    ATH_MSG_VERBOSE("Initializing ...");

    // read-in and cache the content of the VertexPositionsFile, throw error if:
    //  * no VertexPositionsFile is given
    //  * or something goes wrong in the file read-in
    if ( m_vertexPositionFile.empty() || readVertexPosFile().isFailure() ) {
      ATH_MSG_ERROR("Something went wrong with setting up the vertex positioning from file '"
                    << m_vertexPositionFile << "'");
      return StatusCode::FAILURE;
    }

    // if a filename is given for a event number overrides
    //   --> read in this file and cache its content locally
    if ( !m_runEventNumbersFile.empty() && readRunEventNumFile().isFailure()) {
      ATH_MSG_ERROR("Something went wrong with setting up the run/event number overrides from file'"
                    << m_runEventNumbersFile << "'");
      return StatusCode::FAILURE;
    }

    ATH_CHECK(m_eventInfoKey.initialize());

    // everything set up properly
    return StatusCode::SUCCESS;
  }

  /** Athena algtool's Hooks */
  StatusCode VertexPositionFromFile::finalize()
  {
    ATH_MSG_VERBOSE("Finalizing ...");
    return StatusCode::SUCCESS;
  }

  /** read-in and cache vertex positions from file */
  StatusCode VertexPositionFromFile::readVertexPosFile()
  {
    ATH_MSG_INFO("Will read in vertex positions from file.");
    FILE *vfile = fopen( m_vertexPositionFile.value().c_str(),"r");
    if (!vfile) {
      ATH_MSG_ERROR("Could not open vertex position file: " << m_vertexPositionFile);
      return StatusCode::FAILURE;
    }
    ATH_MSG_DEBUG("Opened vertex position file: " << m_vertexPositionFile);
    int          vrun(0);           // run number
    int          vevent(0);         // event number
    double       vx(0.); // vertex coordinates
    double       vy(0.); // vertex coordinates
    double       vz(0.); // vertex coordinates
    unsigned int numReadIn(0);      // number of vertex overrides read in
    // read in file
    while (true) {
      // fill local variables with values given in file
      int r = fscanf(vfile, "%i %i %lf %lf %lf\n", &vrun, &vevent, &vx, &vy, &vz);
      // if read-in was successful
      if (r>0) {
        ATH_MSG_VERBOSE( "Read "<<r<<" vertex position values from file: "<<vrun
                         <<"/"<<vevent<<" "<<vx<<","<<vy<<","<<vz);
        // get the corresponding (#run,#event) entry in the m_vertexPositionMap
        RunEventPair    curRunEvt(vrun, vevent);
        XYZCoordinates &curCoordinates = m_vertexPositionMap[curRunEvt];
        // check if (vrun,vevent) combination already filled
        if ( curCoordinates.size()!=0) {
          ATH_MSG_WARNING( "Already position information for run/event "<<vrun<<"/"<<vevent
                           << ", size=" << curCoordinates.size() );
        }
        else {
          curCoordinates.resize(3);
        }
        // store the (x,y,z) coordinates in the vertexPositionMap:
        curCoordinates[0] = vx;
        curCoordinates[1] = vy;
        curCoordinates[2] = vz;
        // use this trick to only allocate the amount of memory
        // actually used by curXYZ
        //curCoordinates.swap( curCoordinates);
        ++numReadIn;
      }
      // nothing read-in
      else {
        ATH_MSG_VERBOSE("Got "<<r<<" from fscanf, stopping");
        break;
      }
    } // loop over lines in file
    // close file
    fclose(vfile);
    ATH_MSG_VERBOSE("Read " << numReadIn << " vertex position entries from file.");
    return StatusCode::SUCCESS;
  }

  /** read-in and cache run/event number overrides locally for vertex positioning */
  StatusCode VertexPositionFromFile::readRunEventNumFile()
  {
    FILE *vefile = fopen( m_runEventNumbersFile.value().c_str(),"r");
    if ( !vefile) {
      ATH_MSG_ERROR("Could not open vertex positioning run/event number file: "<< m_runEventNumbersFile);
      return StatusCode::FAILURE;
    }
    ATH_MSG_VERBOSE("Opened vertex positioning run/event number file: " << m_runEventNumbersFile);
    //svcMgr.EvtIdModifierSvc.add_modifier(run_nbr=167776, evt_nbr=22, time_stamp=1299948350, lbk_nbr=130, nevts=1)
    int verun(0);     // run number
    int veevent(0);   // event number
    int vetime(0);    // time stamp
    int velbn(0);     // lumi block nr
    int ven(0);       // num events
    int numReadIn(0); // number of read-ins
    // read in file
    while (true) {
      // fill local variables with values given in file
      int r = fscanf(vefile, "svcMgr.EvtIdModifierSvc.add_modifier(run_nbr=%i, evt_nbr=%i, time_stamp=%i, lbk_nbr=%i, nevts=%i)\n",
                     &verun, &veevent, &vetime, &velbn, &ven);

      // if read-in was successful
      if (r>0) {
        ATH_MSG_DEBUG( "Read "<<r<<" vertex positioning run/event values: "
                       <<verun<<"/"<<veevent<<" "<<vetime<<"," <<velbn<<","<<ven );
        // store the run and event number locally
        m_vertexPositionRunNum.push_back( verun);
        m_vertexPositionEventNum.push_back( veevent);
        ++numReadIn;
      }
      // nothing read-in
      else {
        ATH_MSG_VERBOSE("Got "<<r<<" from fscanf, stopping");
        break;
      }
    } // loop over lines in file
    // close file
    fclose(vefile);
    ATH_MSG_VERBOSE("Read " << numReadIn <<" vertex positioning run/event entries from file.");
    return StatusCode::SUCCESS;
  }

  /** modifies (displaces) the given GenEvent */
  CLHEP::HepLorentzVector *VertexPositionFromFile::generate(const EventContext& ctx) const
  {
    unsigned int runNumber(0), eventNumber(0);
    // override the run/event number from file
    if (!m_runEventNumbersFile.empty()) {
      // This works because we iterate over the file exactly once
      static std::atomic<size_t> runEventNumbersIndex(0);
      ATH_MSG_DEBUG("Retrieving event info from event file, position " << runEventNumbersIndex);
      runNumber   = m_vertexPositionRunNum[runEventNumbersIndex];
      eventNumber = m_vertexPositionEventNum[runEventNumbersIndex];
      ++runEventNumbersIndex;
    }
    // use run/event numbers from EventInfo class in storegate
    else {
      ATH_MSG_DEBUG("Retrieving event info from SG");
      SG::ReadHandle<xAOD::EventInfo> eventInfo(m_eventInfoKey, ctx);
      if (eventInfo.isValid()) {
        // read out run/event number
        runNumber   = eventInfo->runNumber();
        eventNumber = eventInfo->eventNumber();
      }
      else {
        ATH_MSG_ERROR("Could not retrieve event info from SG");
        return nullptr;
      }
    }
    ATH_MSG_DEBUG("Got run/event: " << runNumber << "/" << eventNumber);
    // read the (x,y,z) coordinates for the current (run,event)
    const RunEventPair curRunEvtPair(runNumber, eventNumber);
    const XYZCoordinates &updatedVertexPosition = m_vertexPositionMap.at(curRunEvtPair);
    ATH_MSG_DEBUG("Got vertex offset: " << updatedVertexPosition[0] << " " <<
                  updatedVertexPosition[1] << " " << updatedVertexPosition[2]);
    // no (x,y,z) coordinates given for the current (run,event) numbers
    if (updatedVertexPosition.size()!=3) {
      ATH_MSG_ERROR("Vertex position requested, but no info found in map for run/event: " <<
                    runNumber << "/" << eventNumber);
      return nullptr;
    }
    // store the actual vertex offset
    CLHEP::HepLorentzVector *vertexOffset =
      new CLHEP::HepLorentzVector( updatedVertexPosition[0], updatedVertexPosition[1],
                                   updatedVertexPosition[2], 0. );
    // and return it
    return vertexOffset;
  }

} // namespace Simulation
