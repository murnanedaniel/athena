#!/usr/bin/env python
"""
Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
"""

from AthenaPoolUtilities.TPCnvTestConfig import TPCnvTest

if __name__ == "__main__":

    infile = 'aod/AOD-21.0.79/AOD-21.0.79-full.pool.root'
    keys = [
        #xAOD::CaloCluster
        'CaloCalTopoClusters',
        'ForwardElectronClusters',
        'HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex',
        'InDetTrackParticlesAssociatedClusters',
        'MuonClusterCollection',
        'TauPi0Clusters',
        'egammaClusters',
             ]

    TPCnvTest(infile, keys)
