# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from IOVDbSvc.IOVDbSvcConfig import addFoldersSplitOnline


def TRTAlignCondAlgCfg(flags, name="TRTAlignCondAlg", **kwargs):
    """Return a ComponentAccumulator for TRTAlignCondAlg algorithm"""
    from TRT_GeoModel.TRT_GeoModelConfig import TRT_GeoModelCfg
    acc = TRT_GeoModelCfg(flags)
    acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/Calib/DX", "/TRT/Calib/DX"))

    if flags.GeoModel.Align.Dynamic:
        acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/AlignL1/TRT", "/TRT/AlignL1/TRT", className="CondAttrListCollection"))
        acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/AlignL2", "/TRT/AlignL2", className="AlignableTransformContainer"))
        kwargs.setdefault("ReadKeyDynamicGlobal", "/TRT/AlignL1/TRT")
        kwargs.setdefault("ReadKeyDynamicRegular", "/TRT/AlignL2")
    else:
        acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/Align", "/TRT/Align", className="AlignableTransformContainer"))

    kwargs.setdefault("UseDynamicFolders", flags.GeoModel.Align.Dynamic)

    acc.addCondAlgo(CompFactory.TRTAlignCondAlg(name, **kwargs))
    return acc


def TRTStrawCondAlgCfg(flags, name="TRTStrawCondAlg", **kwargs):
    """Return a ComponentAccumulator for TRTStrawCondAlg algorithm"""
    acc = TRTAlignCondAlgCfg(flags)
    from TRT_ConditionsServices.TRT_ConditionsServicesConfig import TRT_StrawStatusSummaryToolCfg
    StrawStatusTool = acc.popToolsAndMerge(TRT_StrawStatusSummaryToolCfg(flags))
    acc.addPublicTool(StrawStatusTool)  # public as it is has many clients to save some memory
    acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/Cond/Status", "/TRT/Cond/Status", className="TRTCond::StrawStatusMultChanContainer"))
    kwargs.setdefault("TRTStrawStatusSummaryTool", StrawStatusTool)
    # Alive straws algorithm
    acc.addCondAlgo(CompFactory.TRTStrawCondAlg(name, **kwargs))
    return acc


def TRTActiveCondAlgCfg(flags, name="TRTActiveCondAlg", **kwargs):
    """Return a ComponentAccumulator for TRTActiveCondAlg algorithm"""
    acc = TRTAlignCondAlgCfg(flags)
    from TRT_ConditionsServices.TRT_ConditionsServicesConfig import TRT_StrawStatusSummaryToolCfg
    StrawStatusTool = acc.popToolsAndMerge(TRT_StrawStatusSummaryToolCfg(flags))
    acc.addPublicTool(StrawStatusTool)  # public as it is has many clients to save some memory
    acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/Cond/Status", "/TRT/Cond/Status", className="TRTCond::StrawStatusMultChanContainer"))
    kwargs.setdefault("TRTStrawStatusSummaryTool", StrawStatusTool)
    acc.addCondAlgo(CompFactory.TRTActiveCondAlg(name, **kwargs))
    return acc


def TRTPhaseCondCfg(flags, name="TRTPhaseCondAlg", **kwargs):
    """Return a ComponentAccumulator for TRTPhaseCondAlg algorithm"""
    acc = ComponentAccumulator()
    from TRT_ConditionsServices.TRT_ConditionsServicesConfig import TRT_CalDbToolCfg
    CalDbTool = acc.popToolsAndMerge(TRT_CalDbToolCfg(flags))
    acc.addPublicTool(CalDbTool)  # public as it is has many clients to save some memory
    kwargs.setdefault("TRTCalDbTool", CalDbTool)
    # Average T0 CondAlg
    acc.merge(addFoldersSplitOnline(flags, "TRT", "/TRT/Onl/Calib/T0", "/TRT/Calib/T0",
                                    className='TRTCond::StrawT0MultChanContainer'))
    acc.addCondAlgo(CompFactory.TRTPhaseCondAlg(name, **kwargs))
    return acc
