#!/bin/sh
#
# art-description: Reco_tf runs on 7TeV collision data with all streams 2011
# art-athena-mt: 4
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena

Reco_tf.py --athenaopts="--threads=8" --inputBSFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/RecJobTransformTests/high_mu-data11_7TeV.00179725.physics_JetTauEtmiss.merge.RAW._lb0021.data --maxEvents 300 --autoConfiguration everything --conditionsTag="COMCOND-BLKPA-RUN1-07" --postExec 'all:from IOVDbSvc.CondDB import conddb;conddb.addOverride("/TRT/Calib/PID_NN", "TRTCalibPID_NN_v1");conddb.addOverride("/TRT/Onl/Calib/PID_NN", "TRTCalibPID_NN_v1");' --outputESDFile myESD.pool.root --outputAODFile myAOD.pool.root --outputHISTFile myHist.root

RES=$?
echo "art-result: $RES Reco"
