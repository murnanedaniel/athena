// Dear emacs, this is -*- c++ -*-
//
// Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
//
#ifndef XAODCORECNV_ROOTHEADERLOADERSVC_H
#define XAODCORECNV_ROOTHEADERLOADERSVC_H

// Local include(s).
#include "xAODCoreCnv/IROOTHeaderLoaderSvc.h"

// Framework include(s).
#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/Property.h"

// System include(s).
#include <string>
#include <vector>

namespace xAODMaker {

   class ROOTHeaderLoaderSvc : public extends< AthService,
                                               IROOTHeaderLoaderSvc > {

   public:
      // Inherit the base class's constructor(s).
      using extends::extends;

      /// Function initialising the service
      virtual StatusCode initialize() override;

      /// @name Implementation of the @c xAODMaker::IEventFormatSvc interface
      /// @{

      /// (Force-)Load one particular header
      virtual StatusCode
      loadHeader( const std::string& headerName ) const override;

      /// @}

   private:
      /// Names of the headers to auto-load during initialisation
      Gaudi::Property< std::vector< std::string > > m_headerNames{ this,
         "HeaderNames", {},
         "Names of the headers to auto-load during initialisation" };

   }; // class ROOTHeaderLoaderSvc

} // namespace xAODMaker

#endif // XAODCORECNV_ROOTHEADERLOADERSVC_H
