/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef BEAMEFFECTS_VERTEXPOSITIONFROMFILE_H
#define BEAMEFFECTS_VERTEXPOSITIONFROMFILE_H 1

// STL includes
#include <string>
#include <utility>
#include <vector>
#include <map>
// FrameWork includes
#include "GaudiKernel/ServiceHandle.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "StoreGate/ReadHandleKey.h"
// EventInfo
#include "xAODEventInfo/EventInfo.h"

#include "HepMC_Interfaces/ILorentzVectorGenerator.h"

namespace Simulation {

  /**
   *  @class VertexPositionFromFile
   *
   *  Returns a VertexShift based on what is given in the input file.
   *
   *  based on:
   *  https://svnweb.cern.ch/trac/atlasoff/browser/Simulation/G4Atlas/G4AtlasUtilities/trunk/src/VertexPositioner.cxx
   *
   *  @author Elmar.Ritsch -at- cern.ch
   */

  class VertexPositionFromFile : public extends<AthAlgTool, ILorentzVectorGenerator> {

    public:
      /** Constructor with parameters */
      VertexPositionFromFile( const std::string& t, const std::string& n, const IInterface* p );

      /** Destructor */
      virtual ~VertexPositionFromFile() = default;

      /** Athena algtool's Hooks */
      StatusCode  initialize() override final;
      StatusCode  finalize() override final;

      /** returns current shift */
      CLHEP::HepLorentzVector *generate(const EventContext& ctx) const override final;

    private:

      /** using typedefs to make code more readable for (#run,#event) -> (x,y,z) mapping */
      typedef std::pair< int, int >                    RunEventPair;
      typedef std::vector< double >                    XYZCoordinates;
      typedef std::map< RunEventPair, XYZCoordinates>  EventCoordinatesMap;

      /** read-in and cache vertex positions from file */
      StatusCode  readVertexPosFile();
      /** read-in and cache run/event number overrides locally for vertex positioning */
      StatusCode  readRunEventNumFile();

      /** set vertex position by file */
      Gaudi::Property<std::string>                     m_vertexPositionFile{this, "VertexPositionsFile", ""};        //!< Name of file containing vertex position overrides
      EventCoordinatesMap             m_vertexPositionMap{};         //!< vertex position for (#run,#event) pairs

      /** run and event number overrides according to file
          (to be used optionally in combination with 'set vertex by file') */
      Gaudi::Property<std::string>                     m_runEventNumbersFile{this, "RunAndEventNumbersFile", ""};       //!< Name of file containing event info overrides for pos overrides
      std::vector<int>                m_vertexPositionRunNum{};      //!< run number override vector
      std::vector<int>                m_vertexPositionEventNum{};    //!< event number override vector
      SG::ReadHandleKey<xAOD::EventInfo>    m_eventInfoKey{this, "EventInfoKey", "EventInfo"};              //!< Name of the EventInfo object in StoreGate
  };

}

#endif //> !BEAMEFECTS_VERTEXPOSITIONFROMFILE_H
