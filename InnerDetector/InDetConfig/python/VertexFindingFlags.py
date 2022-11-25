# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

import AthenaCommon.SystemOfUnits as Units
from AthenaConfiguration.Enums import BeamType, LHCPeriod


def createSecVertexingFlags():
    from AthenaConfiguration.AthConfigFlags import AthConfigFlags
    flags = AthConfigFlags()
    flags.addFlag("Fitter.IterationNumber", 10)
    flags.addFlag("Fitter.MakeExtendedVertex", False)
    flags.addFlag("Fitter.FirstMeasuredPoint", True)
    flags.addFlag("Fitter.Robustness", 6)
    flags.addFlag("Fitter.InputParticleMasses", [139.57018, 139.57018])
    flags.addFlag("Fitter.VertexForConstraint", [0., 0., 0.])
    flags.addFlag("Fitter.CovVrtForConstraint", [
                  0.015*0.015, 0., 0.015*0.015, 0., 0., 10000.*10000.])
    flags.addFlag("Fitter.FirstMeasuredPointLimit", False)
    flags.addFlag("Fitter.usePhiCnst", False)
    flags.addFlag("Fitter.useThetaCnst", False)
    # Track selector tool
    flags.addFlag("TrkSel.maxSiD0", 80.0)
    flags.addFlag("TrkSel.maxTrtD0", 10000.)
    flags.addFlag("TrkSel.maxSiZ0", 10000.)
    flags.addFlag("TrkSel.maxTrtZ0", 10000.)
    flags.addFlag("TrkSel.minPt",
                  lambda pcf: pcf.InDet.Tracking.MainPass.minPT)
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): Ntrt < 15
    flags.addFlag("TrkSel.RatioCut1", 0.0)
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): 15 < Ntrt < 25
    flags.addFlag("TrkSel.RatioCut2", 0.0)
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): Ntrt > 25
    flags.addFlag("TrkSel.RatioCut3", 0.0)
    # e-prob cut for TRT conversion tracks (affects 1TRT, 2TRT, SiTRT) (see also below)
    flags.addFlag("TrkSel.RatioTRT", 0.0)
    # eta bins (10) for eta-dep cuts on TRT conversion tracks
    flags.addFlag("TrkSel.TRTTrksEtaBins", [
                  999, 999, 999, 999, 999, 999, 999, 999, 999, 999])
    # eta-dep e-prob cut for TRT conversion tracks
    flags.addFlag("TrkSel.TRTTrksBinnedRatioTRT", [
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0])
    flags.addFlag("TrkSel.IsConversion", False)
    flags.addFlag("TrkSel.significanceD0_Si", 2.)
    flags.addFlag("TrkSel.RatioV0", 0.95)
    # Track pairs selector
    flags.addFlag("TrkPairSel.MaxFirstHitRadius", 510.)
    flags.addFlag("TrkPairSel.MaxDistBetweenTracks", [4., 90., 30.])
    flags.addFlag("TrkPairSel.MaxEta", [1.2, 0.9, 0.4])
    flags.addFlag("TrkPairSel.MaxInitDistance", [100000., 640., 80.])
    flags.addFlag("TrkPairSel.MinTrackAngle", 0.6)
    # Vertex point estimator
    flags.addFlag("VtxPt.MinDeltaR", [-10000., -10000., -10000.])
    flags.addFlag("VtxPt.MaxDeltaR", [10000., 10000., 10000.])
    flags.addFlag("VtxPt.MaxPhi", [10000., 10000., 10000.])
    # Secondary Vertex post selector
    flags.addFlag("SecVtxPost.MaxChi2Vtx", [15., 50., 25.])
    flags.addFlag("SecVtxPost.MaxInvariantMass", [10000., 50., 25.])
    flags.addFlag("SecVtxPost.MinFitMomentum", [1100., 2000., 2000.])
    flags.addFlag("SecVtxPost.MinRadius", [0., 35., 350.])
    flags.addFlag("SecVtxPost.MinPt", 0.0)       # off
    flags.addFlag("SecVtxPost.MaxdR", -10000.)   # off
    flags.addFlag("SecVtxPost.MaxPhiVtxTrk", 10000.)
    # Single track secondary vertex tool
    flags.addFlag("SingleTrk.MaxBLayerHits", 0)
    flags.addFlag("SingleTrk.MinInitialHitRadius", 0.0)
    flags.addFlag("SingleTrk.MinInitialHitRadius_noBlay", 0.0)
    # e-prob cut for 1TRT and 1Si converisons
    flags.addFlag("SingleTrk.MinRatioOfHLhits", 0.0)
    #  InDetSecVtxFinderTool
    flags.addFlag("Finder.RemoveTrtTracks", True)
    flags.addFlag("Finder.IsConversion", False)
    flags.addFlag("Finder.MinDistVtxHit", -200.)
    flags.addFlag("Finder.MaxDistVtxHit", 50.)
    flags.addFlag("Finder.MinFlightAngle", -0.5)
    flags.addFlag("Finder.MinInitVtxR", 5.)
    return flags


def createEGammaPileUpSecVertexingFlags():
    flags = createSecVertexingFlags()
    # InDetSecVxFitterTool
    flags.Fitter.IterationNumber = 30    # huge number of iterations ?
    flags.Fitter.MakeExtendedVertex = True
    flags.Fitter.FirstMeasuredPoint = True
    flags.Fitter.Robustness = 6
    flags.Fitter.InputParticleMasses = [0.511, 0.511]
    flags.Fitter.VertexForConstraint = [0., 0., 0.]
    flags.Fitter.CovVrtForConstraint = [
        0.015*0.015, 0., 0.015*0.015, 0., 0., 10000.*10000.]
    flags.Fitter.FirstMeasuredPointLimit = True
    flags.Fitter.usePhiCnst = True
    flags.Fitter.useThetaCnst = True
    # Track selector tool
    flags.TrkSel.maxSiD0 = 10000.  # off
    flags.TrkSel.maxTrtD0 = 10000.  # off
    flags.TrkSel.maxSiZ0 = 10000.  # off
    flags.TrkSel.maxTrtZ0 = 10000.  # off
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): Ntrt < 15
    flags.TrkSel.RatioCut1 = 0.0
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): 15 < Ntrt < 25
    flags.TrkSel.RatioCut2 = 0.1
    # e-prob for Si conversion tracks (affects 1Si, 2Si, SiTRT): Ntrt > 25
    flags.TrkSel.RatioCut3 = 0.1
    # e-prob cut for TRT conversion tracks (affects 1TRT, 2TRT, SiTRT) (see also below)
    flags.TrkSel.RatioTRT = 0.51
    # eta bins (10) for eta-dep cuts on TRT conversion tracks
    flags.TrkSel.TRTTrksEtaBins = [
        0.7,   0.8,   0.9,  1.2,  1.3,  1.6,  1.7,  1.8,  1.9,  999]
    # eta-dep e-prob cut for TRT conversion tracks
    flags.TrkSel.TRTTrksBinnedRatioTRT = lambda pcf: (
        [0.51,  0.80,  0.90, 0.80, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51]
        if pcf.GeoModel.Run is LHCPeriod.Run3
        else [0.60,  0.80,  0.90, 0.80, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51])
    flags.TrkSel.IsConversion = True
    flags.TrkSel.significanceD0_Si = -1.    # V0 only cuts
    flags.TrkSel.RatioV0 = -1.    # V0 only cuts
    # Track pairs selector
    # hacky way to determine if TRT only of SI
    flags.TrkPairSel.MaxFirstHitRadius = 500.
    # minimum distance cut at vertex
    flags.TrkPairSel.MaxDistBetweenTracks = [10., 50., 50.]
    # delta cotan(theta) cuts, no cut in barrel for TRT only in code
    flags.TrkPairSel.MaxEta = [0.3, 0.5, 0.5]
    flags.TrkPairSel.MaxInitDistance = [10000.0, 10000.0, 10000.0]  # off
    flags.TrkPairSel.MinTrackAngle = 0.0                          # off, V0 cut
    # Vertex point estimator
    flags.VtxPt.MinDeltaR = [-5., -25., -50.]              # D-R1-R2 min cut
    flags.VtxPt.MaxDeltaR = [5., 10., 10.]                 # D-R1-R2 max cut
    flags.VtxPt.MaxPhi = [0.05, 0.2, 0.2]             # dphi cut at vertex
    # Secondary Vertex post selector
    # MaxChi2Vtx was [40.,100.,80.]
    flags.SecVtxPost.MaxChi2Vtx = [50., 50., 50.]
    # MaxInvariantMass off, was [60.,60.,30.]
    flags.SecVtxPost.MaxInvariantMass = [10000., 10000., 10000.]
    # MinFitMomentum off, was [2000.,2000.,2000.]
    flags.SecVtxPost.MinFitMomentum = [0., 0., 0.]
    flags.SecVtxPost.MinRadius = [20.0, 70.0, 250.0]
    flags.SecVtxPost.MinPt = 0.                           # off
    flags.SecVtxPost.MaxdR = -250.                        # off, -250.
    flags.SecVtxPost.MaxPhiVtxTrk = 0.2
    # Single track secondary vertex tool
    flags.SingleTrk.MaxBLayerHits = 0
    flags.SingleTrk.MinInitialHitRadius = 70.
    flags.SingleTrk.MinInitialHitRadius_noBlay = 120.
    # e-prob cut for 1TRT and 1Si converisons
    flags.SingleTrk.MinRatioOfHLhits = lambda pcf: (
        0.40 if pcf.GeoModel.Run is LHCPeriod.Run3 else 0.51)
    #  InDetSecVtxFinderTool
    flags.Finder.RemoveTrtTracks = False
    flags.Finder.IsConversion = True
    flags.Finder.MinDistVtxHit = -350.0     # off
    flags.Finder.MaxDistVtxHit = 250.0      # off
    flags.Finder.MinFlightAngle = 0.0       # off
    flags.Finder.MinInitVtxR = 0.0          # off
    return flags


def createPriVertexingFlags():
    from AthenaConfiguration.AthConfigFlags import AthConfigFlags
    flags = AthConfigFlags()

    # Turn on the primary vertex reconstruction
    flags.addFlag("doVertexFinding", lambda prevFlags: prevFlags.Beam.Type is not BeamType.Cosmics)
    flags.addFlag("maxAbsEta", 9999.0)
    flags.addFlag("maxD0", 4.0 * Units.mm)
    flags.addFlag("maxNPixelHoles", 1)
    flags.addFlag("maxSigmaD0", 5.0 * Units.mm)
    flags.addFlag("maxSigmaZ0SinTheta", 10.0 * Units.mm)
    flags.addFlag("maxZ0", 1000.0 * Units.mm)
    flags.addFlag("maxZ0SinTheta", 1000.0 * Units.mm)
    flags.addFlag("minNInnermostLayerHits", 0)
    flags.addFlag("minNPixelHits", 1)
    flags.addFlag("minNSctHits", 4)
    flags.addFlag("minNSiHits", 6)
    flags.addFlag("minNTrtHits", 0)
    flags.addFlag("minPt", 500.0 * Units.MeV)
    flags.addFlag("maxZinterval", 3)
    # MaxTracks cuts are specific to the IterativeFinding config
    flags.addFlag("doMaxTracksCut", True)
    flags.addFlag("MaxTracks", 3000)
    # string to store the setup for primary vertexing.
    #  Possible types:
    # 'ActsGaussAdaptiveMultiFinding', 'GaussAdaptiveMultiFinding',
    # 'GaussIterativeFinding', 'AdaptiveMultiFinding', 'IterativeFinding'
    flags.addFlag("setup", "ActsGaussAdaptiveMultiFinding")
    # string to store the type of sorting algorithm to separate signal and pile-up vertices.
    # Possible types:
    # 'SumPt2Sorting','SumPtSorting'
    flags.addFlag("sortingSetup", "SumPt2Sorting")
    flags.addFlag("useBeamConstraint", True)

    return flags


def createITkPriVertexingFlags():
    flags = createPriVertexingFlags()

    flags.minNTrtHits = 0
    flags.minPt = 900.0 * Units.MeV
    flags.maxD0 = 1.0 * Units.mm
    flags.maxSigmaD0 = 0.35 * Units.mm
    flags.maxSigmaZ0SinTheta = 2.5 * Units.mm
    flags.minNPixelHits = 3
    flags.minNSctHits = 0
    flags.minNSiHits = 7
    flags.maxZinterval = 0.5

    return flags
