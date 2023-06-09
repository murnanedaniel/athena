# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

from AthenaCommon.CfgGetter import addTool

# this returns three different tools, depending on the runtime config
addTool("G4CosmicFilter.G4CosmicFilterConfig.getCosmicFilter", "G4CosmicFilter")
addTool("G4CosmicFilter.G4CosmicFilterConfig.getCosmicFilterTool", "G4UA::G4CosmicFilterTool")
#addTool("G4CosmicFilter.G4CosmicFilterConfig.getCosmicAndFilter", "G4CosmicAndFilter")
#addTool("G4CosmicFilter.G4CosmicFilterConfig.getCosmicOrFilter", "G4CosmicOrFilter")
