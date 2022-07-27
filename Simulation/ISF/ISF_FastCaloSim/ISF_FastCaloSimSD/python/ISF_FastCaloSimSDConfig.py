# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaCommon import CfgMgr

def getFCS_StepInfoSensitiveDetector(name="FCS_StepInfoSensitiveDetector", **kwargs):
    ## Main configuration
    kwargs.setdefault("StacVolumes",["LArMgr::LAr::EMB::STAC"])
    kwargs.setdefault("PresamplerVolumes",["LArMgr::LAr::Barrel::Presampler::Module"])
    kwargs.setdefault("NegIWVolumes",["LArMgr::LAr::EMEC::Neg::InnerWheel"])
    kwargs.setdefault("NegOWVolumes",["LArMgr::LAr::EMEC::Neg::OuterWheel"])
    kwargs.setdefault("NegBOBarretteVolumes",["LArMgr::LAr::EMEC::Neg::BackOuterBarrette::Module::Phidiv"])
    kwargs.setdefault("MiniVolumes",["LArMgr::MiniFCAL::Wafer"])
    kwargs.setdefault("PosIWVolumes",["LArMgr::LAr::EMEC::Pos::InnerWheel"])
    kwargs.setdefault("PosOWVolumes",["LArMgr::LAr::EMEC::Pos::OuterWheel"])
    kwargs.setdefault("PosBOBarretteVolumes",["LArMgr::LAr::EMEC::Pos::BackOuterBarrette::Module::Phidiv"])
    kwargs.setdefault("PresVolumes", ["LArMgr::LAr::Endcap::Presampler::LiquidArgon"])
    kwargs.setdefault("SliceVolumes",["LArMgr::LAr::HEC::Module::Depth::Slice"])
    kwargs.setdefault("FCAL1Volumes",["LArMgr::LAr::FCAL::Module1::Gap"])
    kwargs.setdefault("FCAL2Volumes",["LArMgr::LAr::FCAL::Module2::Gap"])
    kwargs.setdefault("FCAL3Volumes",["LArMgr::LAr::FCAL::Module3::Gap"])
    kwargs.setdefault("TileVolumes",["Tile::Scintillator"])
    kwargs.setdefault("OutputCollectionNames", ["EventSteps"])
    # TODO Add extra configuration here!!
    return CfgMgr.FCS_Param__FCS_StepInfoSDTool(name, **kwargs)
