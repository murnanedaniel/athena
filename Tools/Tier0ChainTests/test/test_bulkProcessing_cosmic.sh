#!/bin/sh
#
# art-description: Express processing at Tier0
# art-type: grid
# art-include: master/Athena
# art-include: 23.0/Athena
# art-athena-mt: 8

Reco_tf.py  \
--AMI f1287  \
--inputBSFile="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_cos.00421739.express_express.merge.RAW._lb0534._SFO-12._0001.1" \
--preExec='all:from RecExConfig.RecFlags import rec; rec.doZdc.set_Value_and_Lock(False); from AthenaConfiguration.AllConfigFlags import ConfigFlags; ConfigFlags.Trigger.triggerConfig="DB";  ConfigFlags.DQ.Steering.HLT.doBjet=True; ConfigFlags.DQ.Steering.HLT.doInDet=True; ConfigFlags.DQ.Steering.HLT.doBphys=True; ConfigFlags.DQ.Steering.HLT.doCalo=True; ConfigFlags.DQ.Steering.HLT.doEgamma=False; ConfigFlags.DQ.Steering.HLT.doMET=True; ConfigFlags.DQ.Steering.HLT.doJet=True; ConfigFlags.DQ.Steering.HLT.doMinBias=True; ConfigFlags.DQ.Steering.HLT.doMuon=True; ConfigFlags.DQ.Steering.HLT.doTau=True; ConfigFlags.Tile.doTimingHistogramsForGain=0;' \
--outputAODFile="AOD.pool.root" \
--outputESDFile="ESD.pool.root" \
--outputHISTFile="HIST.root" \
--imf False

rc1=$?
echo "art-result: $rc1 Reco"

rc2=-9999
if [ ${rc1} -eq 0 ]
then
  ArtPackage=$1
  ArtJobName=$2
  art.py compare grid --entries 30 ${ArtPackage} ${ArtJobName} --mode=semi-detailed --order-trees --ignore-exit-code diff-pool
  rc2=$?
fi
echo  "art-result: ${rc2} (against previous nightly)"

rc3=-9999
if [ ${rc1} -eq 0 ]
then
  art.py compare ref . /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3-22.0_references_for_comparison/test_bulkProcessing_cosmic_2022-05-23T2101 \
  --entries 100 --mode=semi-detailed --order-trees --ignore-exit-code diff-pool
  rc3=$?
fi
echo  "art-result: ${rc3} (against reference)"
