#!/bin/sh

# art-include: master/Athena
# art-description: DAOD building PHYS and PHYSLITE data18 MP w/ SharedWriter
# art-type: grid
# art-output: *.pool.root
# art-output: checkFile*.txt
# art-output: checkxAOD*.txt
# art-output: checkIndexRefs*.txt
# art-athena-mt: 8

ATHENA_CORE_NUMBER=8 Reco_tf.py \
  --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data18_13TeV.00357772.physics_Main.recon.AOD.r13286/AOD.27654050._000557.pool.root.1 \
  --outputDAODFile art.pool.root \
  --reductionConf PHYS PHYSLITE \
  --maxEvents -1 \
  --sharedWriter True \
  --multiprocess True \
  --preExec 'from AthenaCommon.DetFlags import DetFlags; DetFlags.detdescr.all_setOff(); DetFlags.BField_setOn(); DetFlags.digitize.all_setOff(); DetFlags.detdescr.Calo_setOn(); DetFlags.simulate.all_setOff(); DetFlags.pileup.all_setOff(); DetFlags.overlay.all_setOff();' \
  --postExec 'from DerivationFrameworkJetEtMiss.JetCommon import swapAlgsInSequence; swapAlgsInSequence(topSequence,"jetalg_ConstitModCorrectPFOCSSKCHS_GPFlowCSSK", "UFOInfoAlgCSSK" )'

echo "art-result: $? reco"

checkFile.py DAOD_PHYS.art.pool.root > checkFile_PHYS.txt

echo "art-result: $?  checkfile PHYS"

checkFile.py DAOD_PHYSLITE.art.pool.root > checkFile_PHYSLITE.txt

echo "art-result: $?  checkfile PHYSLITE"

checkxAOD.py DAOD_PHYS.art.pool.root > checkxAOD_PHYS.txt

echo "art-result: $?  checkxAOD PHYS"

checkxAOD.py DAOD_PHYSLITE.art.pool.root > checkxAOD_PHYSLITE.txt

echo "art-result: $?  checkxAOD PHYSLITE"

checkIndexRefs.py DAOD_PHYS.art.pool.root > checkIndexRefs_PHYS.txt 2>&1

echo "art-result: $?  checkIndexRefs PHYS"

checkIndexRefs.py DAOD_PHYSLITE.art.pool.root > checkIndexRefs_PHYSLITE.txt 2>&1

echo "art-result: $?  checkIndexRefs PHYSLITE"
