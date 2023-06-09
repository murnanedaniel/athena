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

Reco_tf.py --AMI=f739 --imf False --maxEvents=500 --ignoreErrors=False --conditionsTag=CONDBR2-BLKPA-2017-05 --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data16_13TeV.00304008.physics_Main.daq.RAW._lb0838._SFO-6._0001.data --outputDESDM_IDALIGNFile=IDALIGN_Main_0.pool.root --outputDESDM_SLTTMUFile=SLTTMU_Main_0.pool.root --outputDRAW_ZMUMUFile=ZMUMU_Main_0.data --outputDAOD_IDTIDEFile=IDTIDE_Main_0.AOD.pool.root --outputDESDM_MCPFile=MCP_Main_0.pool.root --outputDRAW_TOPSLMUFile=TOPSLMU_Main_0.data --outputDRAW_EGZFile=EGZ_Main_0.data --outputDESDM_PHOJETFile=PHOJET_Main_0.pool.root --outputDESDM_TILEMUFile=TILEMU_Main_0.pool.root --outputDRAW_RPVLLFile=RPVLL_Main_0.data --outputDAOD_IDTRKVALIDFile=IDTRKVALID_Main_0.AOD.pool.root --outputHISTFile=myMergedMonitoring_Main_0.root --outputDESDM_EXOTHIPFile=EXOTHIP_Main_0.pool.root --outputDESDM_SGLELFile=SGLEL_Main_0.pool.root --outputESDFile=myESD_Main_0.pool.root --outputDESDM_EGAMMAFile=EGAMMA_Main_0.pool.root --outputAODFile=myAOD_Main_0.AOD.pool.root --outputTAGFile=myTAG_Main_0.root --outputDRAW_TAUMUHFile=TAUMUH_Main_0.data --outputDESDM_CALJETFile=CALJET_Main_0.pool.root
rc1=$?
echo "art-result: $rc1 Reco"

rc2=-9999
if [ ${rc1} -eq 0 ]
then
    ArtRef=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/$1/TCT_21.0-mc16a_references/$2
    art.py compare ref --entries 2 . $ArtRef
    rc2=$?
fi
echo  "art-result: ${rc2} Diff"


