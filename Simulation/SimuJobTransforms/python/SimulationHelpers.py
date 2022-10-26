# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from SimulationConfig.SimEnums import BeamPipeSimMode, CalibrationRun, CavernBackground, LArParameterization


def getDetectorsFromRunArgs(ConfigFlags, runArgs):
    """Generate detector list based on runtime arguments."""
    if hasattr(runArgs, 'detectors'):
        detectors = runArgs.detectors
    else:
        from AthenaConfiguration.AutoConfigFlags import getDefaultDetectors
        detectors = getDefaultDetectors(ConfigFlags.GeoModel.AtlasVersion, includeForward=False)

    # Support switching on Forward Detectors
    if hasattr(runArgs, 'LucidOn'):
        detectors = detectors+['Lucid']
    if hasattr(runArgs, 'ZDCOn'):
        detectors = detectors+['ZDC']
    if hasattr(runArgs, 'AFPOn'):
        detectors = detectors+['AFP']
    if hasattr(runArgs, 'ALFAOn'):
        detectors = detectors+['ALFA']
    if hasattr(runArgs, 'FwdRegionOn'):
        detectors = detectors+['FwdRegion']
    # TODO here support switching on Cavern geometry
    # if hasattr(runArgs, 'CavernOn'):
    #     detectors = detectors+['Cavern']

    # Fatras does not support simulating the BCM, so have to switch that off
    from SimulationConfig.SimEnums import SimulationFlavour
    if ConfigFlags.Sim.ISF.Simulator in [SimulationFlavour.ATLFASTIIFMT, SimulationFlavour.ATLFASTIIF_G4MS, SimulationFlavour.ATLFAST3F_G4MS]:
        try:
            detectors.remove('BCM')
        except ValueError:
            pass

    return detectors


def enableFrozenShowersFCalOnly(ConfigFlags):
    """Turns on GFlash shower parametrization for FCAL"""
    ConfigFlags.Sim.LArParameterization = LArParameterization.FrozenShowersFCalOnly
    ConfigFlags.Sim.CalibrationRun = CalibrationRun.Off


def enableBeamPipeKill(ConfigFlags):
    ConfigFlags.Sim.BeamPipeCut = 0.
    ConfigFlags.Sim.BeamPipeSimMode = BeamPipeSimMode.FastSim


def enableTightMuonStepping(ConfigFlags):
    ConfigFlags.Sim.TightMuonStepping = True


def enableG4SignalCavern(ConfigFlags):
    """Set ConfigFlags to take care of Neutron BG"""
    ConfigFlags.Sim.CavernBackground = CavernBackground.Signal


def enableCalHits(ConfigFlags):
    """Turns on calibration hits for LAr and Tile"""
    ConfigFlags.Sim.CalibrationRun = CalibrationRun.LArTile
    # deactivate incompatible optimizations
    ConfigFlags.Sim.LArParameterization = LArParameterization.NoFrozenShowers
    ConfigFlags.Sim.NRRThreshold = False
    ConfigFlags.Sim.NRRWeight = False
    ConfigFlags.Sim.PRRThreshold = False
    ConfigFlags.Sim.PRRWeight = False


def enableParticleID(ConfigFlags):
    """Mods to have primary particle barcode signature on for calorimeter calibration hits."""
    ConfigFlags.Sim.ParticleID=True


def enableVerboseSelector(ConfigFlags):
    """ """
    ConfigFlags.Sim.OptionalUserActionList += ['G4DebuggingTools.G4DebuggingToolsConfig.VerboseSelectorToolCfg']
