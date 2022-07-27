// Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

// System include(s):
#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <cstdlib>
#include <sstream>

// ROOT include(s):
#include <TApplication.h>
#include <TClass.h>
#include <TError.h>
#include <TSystem.h>

// Local include(s):
#include "xAODRootAccess/Init.h"

namespace xAOD {

   /// Pointer to the original error handler function if there was one
   static ErrorHandlerFunc_t sErrorHandler = 0;
   /// Function filtering the warnings coming from ROOT
   void ErrorHandler( Int_t level, Bool_t abort, const char* location,
                      const char* message );

   /// Internal status flag
   static bool sInitialised = false;

   /// Width of the message source strings
   static size_t sMessageSourceWidth = 25;

   StatusCode Init( const char* appname ) {

      return Init( appname, 0, 0 );
   }

   StatusCode Init( const char* appname, int* argc, char** argv ) {

      // Check if we need to do anything:
      if( sInitialised ) return StatusCode::SUCCESS;

      // Set up our own error handler function:
      sErrorHandler = ::SetErrorHandler( ErrorHandler );
      if( ! sErrorHandler ) {
         std::cerr << "<xAOD::Init> ERROR Couldn't set up ROOT message "
                   << "filtering" << std::endl;
         return StatusCode::FAILURE;
      }

      // Create an application. This is needed to ensure the auto-loading
      // of the xAOD dictionaries.
      if( ! gApplication ) {
         if( argc && argv ) {
            [[maybe_unused]]
            static ::TApplication sApplication( appname, argc, argv );
         } else {
            ::TApplication::CreateApplication();
         }
      }

      // Load the libraries in a carefully selected order.
      // This is a temporary work-around (26 Oct 20) until the current
      // xAOD dictionary issues are worked out.
      for (const char *name : {
            "xAOD::TruthParticle_v1",
            "xAOD::MuonRoI_v1",
            "xAOD::CaloCluster_v1",
            "xAOD::TrackParticle_v1",
            "xAOD::Electron_v1",
            "xAOD::Muon_v1",
            "xAOD::Jet_v1",
            "xAOD::TauJet_v1",
            "xAOD::PFO_v1",
            "xAOD::TrigElectron_v1",
            "xAOD::L2CombinedMuon_v1",
            "xAOD::Particle_v1"}) {
        // silently ignore missing classes, because this gets used in
        // all projects, and not all projects contain all xAOD classes
        static constexpr Bool_t LOAD = kTRUE;
        static constexpr Bool_t SILENT = kTRUE;
        TClass::GetClass( name, LOAD, SILENT );
      }

      // Let the user know what happened:
      ::Info( appname, "Environment initialised for data access" );

      // Return gracefully:
      sInitialised = true;
      return StatusCode::SUCCESS;
   }

   void SetMessageSourceWidth( size_t value ) {

      sMessageSourceWidth = value;
      return;
   }

   /// All ROOT/xAOD messages are set up to pass through this function.
   /// It is a good place to selectively turn off printing some warnings
   /// coming from ROOT that the xAOD users should just not be bothered
   /// with.
   ///
   /// @param level The message level
   /// @param abort Whether the job needs to abort
   /// @param location The source of the message
   /// @param message The message that needs to be printed
   ///
   void ErrorHandler( Int_t level, Bool_t abort, const char* location,
                      const char* message ) {

      // Check if we need to print anything for this level:
      if( level < gErrorIgnoreLevel ) {
         return;
      }

      // Source of the missing dictionary warnings:
      static const char* const DICT_WARNING_SOURCE = "TClass::Init";

      // Filter out warnings about missing dictionaries. As these are relatively
      // common. In case a problem occurs because of a missing dictionary, some
      // other piece of code will complain anyway.
      if( ( level == kWarning ) &&
          ( ! strcmp( location, DICT_WARNING_SOURCE ) ) ) {
         return;
      }

      // Construct a string version of the message's level:
      const char* msgLevel = 0;
      if( level >= kFatal ) {
         msgLevel = "FATAL  ";
      } else if( level >= kError ) {
         msgLevel = "ERROR  ";
      } else if( level >= kWarning ) {
         msgLevel = "WARNING";
      } else {
         msgLevel = "INFO   ";
      }

      // Make sure that the message's source/location is not longer than a
      // pre-set maximum value:
      std::string source( location );
      if( source.size() > sMessageSourceWidth ) {
         source = source.substr( 0, sMessageSourceWidth - 3 );
         source += "...";
      }

      // Print the message to stdout/std::cout:
      std::ostringstream output;
      output << std::setiosflags( std::ios::left )
             << std::setw( sMessageSourceWidth ) << source << " " << msgLevel
             << " " << message;
      std::cout << output.str() << std::endl;

      // If we don't need to abort, return now:
      if( ! abort ) {
         return;
      }

      // Abort with a stack trace if possible:
      std::cout << std::endl << "Aborting..." << std::endl;
      if( gSystem ) {
         gSystem->StackTrace();
         gSystem->Abort();
      } else {
         std::abort();
      }

      return;
   }

} // namespace xAOD
