# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.AthConfigFlags import AthConfigFlags
from AthenaConfiguration.AutoConfigFlags import GetFileMD, DetDescrInfo
from AthenaConfiguration.Enums import LHCPeriod, ProductionStep, Project

def createGeoModelConfigFlags():
    gcf=AthConfigFlags()

    gcf.addFlag('GeoModel.Layout', 'atlas') # replaces global.GeoLayout

    gcf.addFlag("GeoModel.AtlasVersion",
                lambda prevFlags : (GetFileMD(prevFlags.Input.Files).get("GeoAtlas", None)
                                    or "ATLAS-R2-2016-01-00-01"))

    gcf.addFlag("GeoModel.Align.Dynamic",
                lambda prevFlags : not prevFlags.Input.isMC and prevFlags.Common.ProductionStep not in [ProductionStep.Simulation, ProductionStep.Overlay])

    gcf.addFlag("GeoModel.Align.LegacyConditionsAccess",
                lambda prevFlags : prevFlags.Common.Project is Project.AthSimulation or prevFlags.Common.ProductionStep is ProductionStep.Simulation)
                # Mainly for G4 which still loads alignment on initialize

    gcf.addFlag("GeoModel.Run",
                lambda prevFlags : LHCPeriod(DetDescrInfo(prevFlags.GeoModel.AtlasVersion)['Common']['Run']),
                enum=LHCPeriod)
                # Based on CommonGeometryFlags.Run

    gcf.addFlag("GeoModel.Type",
                lambda prevFlags : DetDescrInfo(prevFlags.GeoModel.AtlasVersion)['Common']['GeoType'])
                # Geometry type in {ITKLoI, ITkLoI-VF, etc...}

    gcf.addFlag("GeoModel.IBLLayout",
                lambda prevFlags : DetDescrInfo(prevFlags.GeoModel.AtlasVersion)['Pixel']['IBLlayout'])
                # IBL layer layout  in {"planar", "3D", "noIBL"}

    return gcf
