#!/bin/sh
#
# art-description: RecoTrf
# art-type: grid
# art-include: 21.0/Athena
# art-include: 21.0-mc16a/Athena
# art-include: 21.0-mc16d/Athena
# art-include: 21.0-TrigMC/Athena
# art-include: master/Athena
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-output: log.*

Reco_tf.py  --outputHISTFile=myMergedMonitoring_express_3.root --AMI=f628 --maxEvents=30 --outputESDFile=myESD_express_3.pool.root --outputAODFile=myAOD_express_3.AOD.pool.root --outputTAGFile=myTAG_express_3.root --ignoreErrors=False --steering=doRAWtoALL --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data15_13TeV.00276689.express_express.merge.RAW._lb0226._SFO-ALL._0001.1,/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data15_13TeV.00276689.express_express.merge.RAW._lb0227._SFO-ALL._0001.1 --imf False
echo "art-result: $? Reco"

ArtRef=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/$1/TCT_21.0-mc16a_references/$2
art.py compare ref --entries 10 . $ArtRef
echo "art-result: $? Diff"
