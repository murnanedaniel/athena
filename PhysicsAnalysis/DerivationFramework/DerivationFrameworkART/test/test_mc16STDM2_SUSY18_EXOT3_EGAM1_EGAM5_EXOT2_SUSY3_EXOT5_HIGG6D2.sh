#!/bin/sh

# art-include
# art-description: DAOD building STDM2 SUSY18 EXOT3 EGAM1 EGAM5 EXOT2 SUSY3 EXOT5 HIGG6D2 mc16
# art-type: grid
# art-output: *.pool.root

set -e

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/AOD.14795494._005958.pool.root.1 --outputDAODFile art.pool.root --reductionConf STDM2 SUSY18 EXOT3 EGAM1 EGAM5 EXOT2 SUSY3 EXOT5 HIGG6D2 --maxEvents 500  --preExec 'rec.doApplyAODFix.set_Value_and_Lock(True);from BTagging.BTaggingFlags import BTaggingFlags;BTaggingFlags.CalibrationTag = "BTagCalibRUN12-08-47" '  --passThrough True 
