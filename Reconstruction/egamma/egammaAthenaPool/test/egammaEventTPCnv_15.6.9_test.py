#!/usr/bin/env python
"""
Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
"""

from AthenaPoolUtilities.TPCnvTestConfig import TPCnvTest

if __name__ == "__main__":
    infile = 'aod/AOD-15.6.9/AOD-15.6.9-full.pool.root'
    keys = [
        #ElectronContainer_p2
        'AtlfastIsoElectronCollection',
        'AtlfastElectronCollection',
        'ElectronAODCollection',

        #PhotonContainer_p2
        'AtlfastIsoPhotonCollection',
        'AtlfastPhotonCollection',
        'PhotonAODCollection',

        #egammaContainer_p2
        'HLT_egamma',
        'HLT_egamma_Photons',
        'HLT_egamma_Electrons',

        #egDetailContainer_p2
        # Dumped as part of the above.
    ]

    TPCnvTest(infile, keys, useGeoModelSvc=True, doPixel=True, doSCT=True, doTRT=True, doMuon=True, doTracks=True)

