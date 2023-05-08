#!/bin/sh
#
# art-description: HIP RecoTrf on 2013 pPb 5 TeV data, HardProbes stream (Run218048 LB1102)
# art-type: grid
# art-include: 21.0/Athena
# art-include: 21.0-mc16a/Athena
# art-include: 21.0-mc16d/Athena
# art-include: 21.0-TrigMC/Athena
# art-include: master/Athena
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-output: log.*

Reco_tf.py --outputHISTFile=myMergedMonitoring_0.root --outputESDFile=myESD_0.pool.root --outputAODFile=myAOD_0.AOD.pool.root --outputTAGFile=myTAG_0.root --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data13_hip.00218048.physics_HardProbes.merge.RAW._lb1102._SFO-8._0001.1 --AMI=q315 --maxEvents=100 --conditionsTag "all:COMCOND-BLKPA-RUN1-09" --geometryVersion "all:ATLAS-R1-2012-02-00-00" --preExec "all:rec.doHIP.set_Value_and_Lock(True);DQMonFlags.doLVL1CaloMon=False;DQMonFlags.doCTPMon=False" --steering=doRAWtoALL --imf False
echo "art-result: $? Reco"

ArtRef=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/$1/TCT_21.0-mc16a_references/$2
art.py compare ref --entries 10 . $ArtRef
echo "art-result: $? Diff"
