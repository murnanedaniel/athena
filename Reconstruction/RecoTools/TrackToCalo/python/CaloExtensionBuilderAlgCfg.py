# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentFactory import CompFactory

def CaloExtensionBuilderAlgCfg(flags, **kwargs):
    from TrkConfig.AtlasExtrapolatorConfig import AtlasExtrapolatorCfg    
    result = AtlasExtrapolatorCfg(flags)
    kwargs.setdefault("LastCaloExtentionTool", CompFactory.Trk.ParticleCaloExtensionTool(Extrapolator=result.popPrivateTools()))

    # P->T conversion extra dependencies
    if flags.Detector.GeometryITk:
        kwargs.setdefault("ExtraInputs", [
            ("InDetDD::SiDetectorElementCollection", "ConditionStore+ITkPixelDetectorElementCollection"),
            ("InDetDD::SiDetectorElementCollection", "ConditionStore+ITkStripDetectorElementCollection"),
        ])
    else:
        kwargs.setdefault("ExtraInputs", [
            ("InDetDD::SiDetectorElementCollection", "ConditionStore+PixelDetectorElementCollection"),
            ("InDetDD::SiDetectorElementCollection", "ConditionStore+SCT_DetectorElementCollection"),
        ])

    result.addEventAlgo(CompFactory.Trk.CaloExtensionBuilderAlg(**kwargs))
    return result
