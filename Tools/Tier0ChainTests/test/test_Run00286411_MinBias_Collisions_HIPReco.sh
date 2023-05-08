#!/bin/sh
#
# art-description: HIP RecoTrf on 2015 pp 5 TeV data, MinBias stream (Run286411 LB660)
# art-type: grid
# art-include: 21.0/Athena
# art-include: 21.0-mc16a/Athena
# art-include: 21.0-mc16d/Athena
# art-include: 21.0-TrigMC/Athena
# art-include: master/Athena
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-output: log.*

Reco_tf.py --outputHISTFile=myMergedMonitoring_0.root --outputESDFile=myESD_0.pool.root --outputAODFile=myAOD_0.AOD.pool.root --outputTAGFile=myTAG_0.root --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data15_5TeV.00286411.physics_MinBias.daq.RAW._lb0660._SFO-4._0001.data --AMI=q315 --maxEvents=25 --imf False
echo "art-result: $? Reco"

ArtRef=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/$1/TCT_21.0-mc16a_references/$2
art.py compare ref --entries 10 . $ArtRef
echo "art-result: $? Diff"

