#!/bin/sh
#
# art-description: Reco_tf runs with AMI configTag q122
# art-athena-mt: 4
# art-type: grid
# art-include: 22.0/Athena
# art-include: master/Athena

Reco_tf.py --athenaopts="--threads=8" --AMI=q122 --DataRunNumber 00191920 --outputESDFile=myESD.pool.root --outputAODFile=myAOD.pool.root --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/RecJobTransformTests/data11_7TeV.00191920.physics_JetTauEtmiss.merge.RAW._lb0257._SFO-9._0001.1.10evts --postExec 'all:from IOVDbSvc.CondDB import conddb;conddb.addOverride("/TRT/Calib/PID_NN", "TRTCalibPID_NN_v1");conddb.addOverride("/TRT/Onl/Calib/PID_NN", "TRTCalibPID_NN_v1");'

RES=$?
echo "art-result: $RES Reco"
