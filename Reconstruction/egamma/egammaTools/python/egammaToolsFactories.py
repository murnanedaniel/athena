# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

__doc__ = """ToolFactories to instantiate
all egammaTools with default configuration"""
__author__ = "Bruno Lenzi"


from .EMPIDBuilderBase import EMPIDBuilderPhotonBase
from .EMPIDBuilderBase import EMPIDBuilderElectronBase
from ElectronPhotonSelectorTools import ElectronPhotonSelectorToolsConf
from egammaTrackTools.egammaTrackToolsFactories import EMExtrapolationTools
from egammaMVACalib.egammaMVACalibFactories import egammaMVASvc


from egammaTools import egammaToolsConf
from egammaRec.Factories import ToolFactory
from egammaRec import egammaKeys
# to set jobproperties.egammaRecFlags
from egammaRec.egammaRecFlags import jobproperties


_clusterTypes = dict(
    Ele35='ele35', Ele55='ele55', Ele37='ele37',
    Gam35='gam35_unconv', Gam55='gam55_unconv', Gam37='gam37_unconv',
    Econv35='gam35_conv', Econv55='gam55_conv', Econv37='gam37_conv'
)


# Configure corrections for superclusters.
def configureSuperClusterCorrections(swTool):
    """Add attributes ClusterCorrectionToolsXX to egammaSwTool
       object for corrections for superclusters."""
    from CaloClusterCorrection.CaloSwCorrections import make_CaloSwCorrections
    from CaloRec.CaloRecMakers import _process_tools

    for attrName, clName in _clusterTypes.items():
        n = 'ClusterCorrectionToolsSuperCluster' + attrName
        if not hasattr(swTool, n) or getattr(swTool, n):
            continue

        setattr(swTool, n, _process_tools(
            swTool,
            make_CaloSwCorrections(
                clName,
                suffix='EGSuperCluster',
                version=jobproperties.egammaRecFlags.superClusterCorrectionVersion(),
                cells_name=egammaKeys.caloCellKey())))


egammaSwSuperClusterTool = ToolFactory(
    egammaToolsConf.egammaSwTool,
    postInit=[configureSuperClusterCorrections])


EMClusterTool = ToolFactory(
    egammaToolsConf.EMClusterTool,
    OutputClusterContainerName=egammaKeys.outputClusterKey(),
    MVACalibSvc=egammaMVASvc
)


EMConversionBuilder = ToolFactory(
    egammaToolsConf.EMConversionBuilder,
    ConversionContainerName=egammaKeys.outputConversionKey(),
    ExtrapolationTool=EMExtrapolationTools)

EGammaAmbiguityTool = ToolFactory(
    ElectronPhotonSelectorToolsConf.EGammaAmbiguityTool)

EMFourMomBuilder = ToolFactory(egammaToolsConf.EMFourMomBuilder)

egammaLargeClusterMakerTool = ToolFactory(
    egammaToolsConf.egammaLargeClusterMaker,
    name="egammaLCMakerTool",
    InputClusterCollection=egammaKeys.ClusterKey(),
    CellsName=egammaKeys.caloCellKey()
)

egammaLargeFWDClusterMakerTool = ToolFactory(
    egammaToolsConf.egammaLargeClusterMaker,
    name="egammaLCFWDMakerTool",
    InputClusterCollection=egammaKeys.FwdClusterKey(),
    CellsName=egammaKeys.caloCellKey(),
    doFWDelesurraundingWindows=True
)

# Electron Selectors
ElectronPIDBuilder = ToolFactory(
    EMPIDBuilderElectronBase,
    name="ElectronPIDBuilder")

# Photon Selectors
PhotonPIDBuilder = ToolFactory(
    EMPIDBuilderPhotonBase,
    name="PhotonPIDBuilder")

# -------------------------

# Import the factories that are not defined here
from .EMShowerBuilder import EMShowerBuilder            # noqa: F401
