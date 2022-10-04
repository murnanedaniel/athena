#!/bin/bash
#
# art-description: Reco_tf.py q449 RAWtoALL in MT and ComponentAccumulator mode
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-athena-mt: 8

mkdir ca
cd ca
Reco_tf.py --CA \
  --AMI q449  \
  --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_13p6TeV.00431493.physics_Main.daq.RAW._lb0525._SFO-16._0001.data \
  --conditionsTag=CONDBR2-BLKPA-2022-07 \
  --geometryVersion=ATLAS-R3S-2021-03-00-00 \
  --multithreaded="True" \
  --steering "doRAWtoALL" \
  --outputAODFile myAOD_ca.pool.root \
  --outputESDFile myESD_ca.pool.root \
  --outputHISTFile myHIST_ca.root \
  --preExec "all:from AthenaConfiguration.AllConfigFlags import ConfigFlags; ConfigFlags.Jet.WriteToAOD=True; ConfigFlags.MET.WritetoAOD=True; from AthenaConfiguration.AllConfigFlags import ConfigFlags; ConfigFlags.Trigger.triggerConfig='DB'" \
  --imf="False" \
  --maxEvents 100

rc1=$?
echo "art-result: ${rc1} Reco_tf_q449_r2a_ca_mt" 

cd ..
mkdir def
cd def
Reco_tf.py \
  --AMI q449  \
  --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_13p6TeV.00431493.physics_Main.daq.RAW._lb0525._SFO-16._0001.data \
  --conditionsTag=CONDBR2-BLKPA-2022-07 \
  --geometryVersion=ATLAS-R3S-2021-03-00-00 \
  --multithreaded="True" \
  --steering "doRAWtoALL" \
  --outputAODFile myAOD_def.pool.root \
  --outputESDFile myESD_def.pool.root \
  --outputHISTFile myHIST_def.root \
  --preExec "all:from JetRec.JetRecFlags import jetFlags; jetFlags.writeJetsToAOD.set_Value_and_Lock(True); from METReconstruction.METRecoFlags import metFlags; metFlags.WriteMETAssocToOutput.set_Value_and_Lock(True); metFlags.WriteMETToOutput.set_Value_and_Lock(True); from AthenaConfiguration.AllConfigFlags import ConfigFlags; ConfigFlags.Trigger.triggerConfig='DB'" \
  --imf="False" \
  --maxEvents 100

rc2=$?
echo "art-result: ${rc2} Reco_tf_q449_r2a_mt" 

cd ..

# Check for FPEs in the logiles
# test_trf_check_fpe.sh      # currently disabled since FPEAuditor does not work in CA mode
#fpeStat=$?

#echo "art-result: ${fpeStat} FPEs in logfiles"

echo "============ checkxAOD myAOD_ca.pool.root"
checkxAOD ca/myAOD_ca.pool.root
rc3=$?
echo "art-result: ${rc3} checkxAOD myAOD_ca.pool.root"

echo "============ checkxAOD myAOD_def.pool.root"
checkxAOD def/myAOD_def.pool.root
rc4=$?
echo "art-result: ${rc4} checkxAOD myAOD_def.pool.root"

echo "============ xAODDigest.py --extravars myAOD_ca.pool.root"
xAODDigest.py --extravars ca/myAOD_ca.pool.root myAOD_ca.txt
rc5=$?
echo "art-result: ${rc5} xAODDigest.py --extravars myAOD_ca.pool.root"
echo "============ myAOD_ca.txt"
cat myAOD_ca.txt
echo "============ myAOD_ca.txt"

echo "============ xAODDigest.py --extravars myAOD_def.pool.root"
xAODDigest.py --extravars def/myAOD_def.pool.root myAOD_def.txt
rc6=$?
echo "art-result: ${rc6} xAODDigest.py --extravars myAOD_def.pool.root"
echo "============ myAOD_def.txt"
cat myAOD_def.txt
echo "============ myAOD_def.txt"

echo "============ comparexAODDigest.py myAOD_def.txt myAOD_ca.txt"
comparexAODDigest.py myAOD_def.txt myAOD_ca.txt
rc7=$?
echo "art-result: ${rc7} comparexAODDigest.py myAOD_def.txt myAOD_ca.txt"

echo "============ hist_diff.sh ca/myHIST_ca.root def/myHIST.root -i -x (TIME_execute|LAr/Coverage|MismatchEventNumbers|L1Calo/Overview/Errors)"
hist_diff.sh ca/myHIST_ca.root def/myHIST_def.root -i -x "(TIME_execute|LAr/Coverage|MismatchEventNumbers|L1Calo/Overview/Errors)"
rc8=$?
echo "art-result: ${rc8} hist_diff.sh ca/myHIST_ca.root def/myHIST.root"

echo "============ diff-root def/myAOD_def.pool.root ca/myAOD_ca.pool.root"
acmd.py diff-root def/myAOD_def.pool.root ca/myAOD_ca.pool.root --nan-equal --error-mode resilient --ignore-leaves RecoTimingObj_p1_HITStoRDO_timings RecoTimingObj_p1_RAWtoESD_mems RecoTimingObj_p1_RAWtoESD_timings RAWtoESD_mems RAWtoESD_timings ESDtoAOD_mems ESDtoAOD_timings HITStoRDO_timings RAWtoALL_mems RAWtoALL_timings RecoTimingObj_p1_RAWtoALL_mems RecoTimingObj_p1_RAWtoALL_timings RecoTimingObj_p1_EVNTtoHITS_timings EVNTtoHITS_timings RecoTimingObj_p1_Bkg_HITStoRDO_timings index_ref  --order-trees --entries 50 --mode semi-detailed
rc9=$?
echo "art-result: ${rc9} diff-root def/myAOD_def.pool.root ca/myAOD_ca.pool.root"
