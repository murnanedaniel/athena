#!/bin/bash
#
# art-description: Reco_tf.py data17 RAWtoALL in MT mode
# art-type: grid
# art-include: master/Athena
# art-include: 22.0-mc20/Athena
# art-include: 22.0/Athena
# art-athena-mt: 8

timeout 64800 Reco_tf.py \
  --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_13TeV.00330470.physics_Main.daq.RAW._lb0310._SFO-1._0001.data \
  --outputAODFile=myAOD.pool.root \
  --outputHISTFile=myHIST.root \
  --outputDESDM_MCPFile=myDESDM_MCP.pool.root \
  --outputDRAW_ZMUMUFile=myDRAW_ZMUMU.data \
  --outputDAOD_IDTIDEFile=myIDTIDE.pool.root \
  --multithreaded='True' \
  --preExec 'all:from AthenaMonitoring.DQMonFlags import DQMonFlags; DQMonFlags.doHLTMon=False' \
  --postExec 'FPEAuditor.NStacktracesOnFPE=10' \
  --autoConfiguration='everything' \
  --conditionsTag 'all:CONDBR2-BLKPA-RUN2-09' --geometryVersion='default:ATLAS-R2-2016-01-00-01' \
  --runNumber='357750' --steering='doRAWtoALL' --maxEvents='-1'

rc1=$?
echo "art-result: ${rc1} Reco_tf_data17_mt"

# Check for FPEs in the logiles
test_trf_check_fpe.sh
fpeStat=$?

echo "art-result: ${fpeStat} FPEs in logfiles"
