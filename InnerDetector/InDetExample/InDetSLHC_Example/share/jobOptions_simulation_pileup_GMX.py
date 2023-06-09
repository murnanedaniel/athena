#--------------------------------------------------------------
# Template jobOptions: SLHC 
#  - Simulation of pileup (minbias) events
#--------------------------------------------------------------

execfile("ID_only.py")
execfile("geometry.py")
# SimFlags.EventFilter.set_Value({'VertexPositioners': True, 'VertexRangeChecker': True, 'EtaPhiFilters': False})
# SimFlags.KinematicsMode = 'EventGenerator'
# SimFlags.GeneratorPath = 'jobOptions_pythia_SLHC.py'

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
athenaCommonFlags.PoolHitsOutput = "pileup_slhc_Hits.pool.root"
athenaCommonFlags.EvtMax = 1000
# Don't look for a file to read events from:
#athenaCommonFlags.PoolEvgenInput.set_Off()
athenaCommonFlags.PoolEvgenInput = ["mc12_14TeV.119996.Pythia8_A2MSTW2008LO_minbias_inelastic_high.evgen.EVNT.e1133/"]

## Algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
topSeq = AlgSequence()

from G4AtlasApps.PyG4Atlas import PyG4AtlasAlg
topSeq += PyG4AtlasAlg()

MessageSvc.OutputLevel = WARNING
MessageSvc.defaultLimit = 500
MessageSvc.Format = "% F%40W%S%7W%R%T %0W%M"

from G4AtlasApps import SimKernel

