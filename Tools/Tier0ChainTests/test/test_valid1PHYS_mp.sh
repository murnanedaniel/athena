#!/bin/sh

# art-include: master/Athena
# art-include: 22.0-mc20/Athena
# art-description: DAOD building PHYS valid1 MP
# art-type: grid
# art-output: *.pool.root
# art-output: checkFile.txt
# art-output: checkxAOD.txt
# art-output: checkIndexRefs*.txt
# art-athena-mt: 8

set -e

Reco_tf.py --athenaopts='--nprocs=2' --athenaMPMergeTargetSize 'DAOD_*:0' --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/valid1.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.AOD.e5112_s3227_r13145/AOD.27253804._000061.pool.root.1 --outputDAODFile art.pool.root --reductionConf PHYS --maxEvents -1 --preExec 'from AthenaCommon.DetFlags import DetFlags; DetFlags.detdescr.all_setOff(); DetFlags.BField_setOn(); DetFlags.digitize.all_setOff(); DetFlags.detdescr.Calo_setOn(); DetFlags.simulate.all_setOff(); DetFlags.pileup.all_setOff(); DetFlags.overlay.all_setOff();'

echo "art-result: $? reco"

DAODMerge_tf.py --inputDAOD_PHYSFile DAOD_PHYS.art.pool.root --outputDAOD_PHYS_MRGFile art_merged.pool.root

echo "art-result: $? merge"

checkFile.py DAOD_PHYS.art.pool.root > checkFile.txt

echo "art-result: $?  checkfile"

checkxAOD.py DAOD_PHYS.art.pool.root > checkxAOD.txt

echo "art-result: $?  checkxAOD"

checkIndexRefs.py DAOD_PHYS.art.pool.root > checkIndexRefs_PHYS.txt 2>&1

echo "art-result: $?  checkIndexRefs PHYS"
