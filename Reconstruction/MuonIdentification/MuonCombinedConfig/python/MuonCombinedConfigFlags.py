# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.AthConfigFlags import AthConfigFlags
from AthenaConfiguration.Enums import BeamType

def createMuonCombinedConfigFlags(): 
    mcf=AthConfigFlags()
    # This is based on the following from the old configuration:
    # https://gitlab.cern.ch/atlas/athena/blob/release/22.0.8/Reconstruction/MuonIdentification/MuonCombinedRecExample/python/MuonCombinedRecFlags.py
    mcf.addFlag("MuonCombined.doCosmicSplitTracks",False)
    mcf.addFlag("MuonCombined.doMuGirl",True)
    mcf.addFlag("MuonCombined.doCombinedFit",True)
    mcf.addFlag("MuonCombined.doStatisticalCombination",lambda prevFlags: prevFlags.Beam.Type is BeamType.Collisions)
    mcf.addFlag("MuonCombined.doMuonSegmentTagger",True)
    # 'silicon-associated'muons, or muons which rely on special ID reconstruction because they're outside the usual acceptance.
    mcf.addFlag("MuonCombined.doSiAssocForwardMuons",lambda prevFlags : prevFlags.Detector.GeometryID)      
    # Switch on/off algorithms that make Muons for the CaloMuonCollection
    mcf.addFlag("MuonCombined.doCaloTrkMuId",True) 
    # Switch on/off algorithms that make Muons for the MuGirlLowBetaMuonCollection         
    mcf.addFlag("MuonCombined.doMuGirlLowBeta",lambda prevFlags : prevFlags.MuonCombined.doMuGirl)
    mcf.addFlag("MuonCombined.writeUnAssocSegments", True)

    return mcf
    