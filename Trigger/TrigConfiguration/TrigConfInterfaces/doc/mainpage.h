/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
   @mainpage TrigConfInterfaces Package

   @author Joerg Stelzer        <Joerg.Stelzer@cern.ch>
   @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>

   $Revision: 580580 $
   $Date: 2014-01-29 09:44:43 +0100 (Wed, 29 Jan 2014) $

   @section TrigConfInterfacesOverview Overview

   This package holds all the interfaces through which the trigger
   configuration services should be accessed both in- and outside of
   Athena.

   @section TrigConfInterfacesClasses Classes

   The following are the interfaces that can be used outside of Athena:
     - TrigConf::IILVL1ConfigSvc: Interface to services providing LVL1
       configuration data.
     - TrigConf::IIHLTConfigSvc: Interface to services providing HLT
       configuration data.

   The following are the Athena-specific interfaces that should be used
   as the template argument for ServiceHandle objects in order to retrieve
   configuration services in Athena:
     - TrigConf::ILVL1ConfigSvc: Interface providing LVL1 configuration data.
     - TrigConf::IHLTConfigSvc: Interface providing HLT configuration data.
     - TrigConf::IL1TopoConfigSvc: Interface providing the configuration of
       the L1Topo hardware.
     - TrigConf::ITrigConfigSvc: Interface providing all aspects of the
       trigger configuration for offline usage. This is the most often used
       interface for analysis purposes.

   @htmlinclude used_packages.html

   @include requirements
*/
