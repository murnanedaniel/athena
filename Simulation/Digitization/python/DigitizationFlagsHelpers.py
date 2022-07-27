#!/usr/bin/env python
"""Digitization legacy configuration helpers

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""


def setupDigitizationLegacyDetectorFlags(detectors):
    """Setup digitization legacy detector flags"""
    from AthenaCommon.DetFlags import DetFlags

    # Truth is always on
    DetFlags.Truth_setOn()

    # Other subdetectors
    if not detectors or 'BCM' in detectors or 'ID' in detectors:
        DetFlags.BCM_setOn()
    if not detectors or 'DBM' in detectors or 'ID' in detectors:
        DetFlags.DBM_setOn()
    if not detectors or 'Pixel' in detectors or 'ID' in detectors:
        DetFlags.pixel_setOn()
    if not detectors or 'SCT' in detectors or 'ID' in detectors:
        DetFlags.SCT_setOn()
    if not detectors or 'TRT' in detectors or 'ID' in detectors:
        DetFlags.TRT_setOn()
    if not detectors or 'LAr' in detectors or 'Calo' in detectors or 'L1Calo' in detectors:
        DetFlags.LAr_setOn()
    if not detectors or 'Tile' in detectors or 'Calo' in detectors or 'L1Calo' in detectors:
        DetFlags.Tile_setOn()
    if not detectors or 'L1Calo' in detectors:
        DetFlags.LVL1_setOn()
    if not detectors or 'CSC' in detectors or 'Muon' in detectors:
        DetFlags.CSC_setOn()
    if not detectors or 'MDT' in detectors or 'Muon' in detectors:
        DetFlags.MDT_setOn()
    if not detectors or 'RPC' in detectors or 'Muon' in detectors:
        DetFlags.RPC_setOn()
    if not detectors or 'TGC' in detectors or 'Muon' in detectors:
        DetFlags.TGC_setOn()
    if not detectors or 'sTGC' in detectors or 'Muon' in detectors:
        DetFlags.sTGC_setOn()
    if not detectors or 'MM' in detectors or 'Muon' in detectors:
        DetFlags.MM_setOn()

    return DetFlags
