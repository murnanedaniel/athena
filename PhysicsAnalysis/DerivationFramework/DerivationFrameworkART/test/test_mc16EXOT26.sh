#!/bin/sh

# art-include: 21.2/AthDerivation
# art-description: DAOD building EXOT26 mc16
# art-type: grid
# art-output: *.pool.root
# art-output: checkFile.txt
# art-output: checkxAOD.txt

set -e

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/AOD.14795494._005958.pool.root.1 --outputDAODFile art.pool.root --reductionConf EXOT26 --maxEvents -1 --preExec 'rec.doApplyAODFix.set_Value_and_Lock(True);from BTagging.BTaggingFlags import BTaggingFlags;BTaggingFlags.CalibrationTag = "BTagCalibRUN12-08-40" ' 

echo "art-result: $? reco"

DAODMerge_tf.py --inputDAOD_EXOT26File DAOD_EXOT26.art.pool.root --outputDAOD_EXOT26_MRGFile art_merged.pool.root

echo "art-result: $? merge"

checkFile.py DAOD_EXOT26.art.pool.root > checkFile.txt

echo "art-result: $?  checkfile"

checkxAOD.py DAOD_EXOT26.art.pool.root > checkxAOD.txt

echo "art-result: $?  checkxAOD"
