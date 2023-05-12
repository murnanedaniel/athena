#!/bin/sh
#
# art-description: RecoTrf
# art-type: grid
# art-include: master/Athena
# art-include: 22.0-mc20/Athena
# art-output: log.*

Reco_tf.py \
--AMI=q441 \
--athenaopts='--threads=1' \
--conditionsTag 'all:OFLCOND-MC16-SDR-RUN2-09' \
--inputHITSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/q441/22.0/HITS.12560240._000299.pool.root.1 \
--inputRDO_BKGFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/q441/22.0/RDO.17190395._000013.pool.root.1 \
--maxEvents=100 \
--outputAODFile=q441.AOD.pool.root --outputESDFile=q441.ESD.pool.root

#This test currently has the muon isolation reconstruction switched off. It should be switched back on at a later date. 
echo "art-result: $?"
