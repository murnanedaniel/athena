#!/bin/sh

# art-description: DAOD building HIGG4D2 mc16
# art-type: grid
# art-output: *.pool.root

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/AOD.11866988._000378.pool.root.1 --outputDAODFile art.pool.root --reductionConf HIGG4D2 --maxEvents 5000

DAODMerge_tf.py --maxEvents 5 --inputDAOD_HIGG4D2File DAOD_HIGG4D2.art.pool.root --outputDAOD_HIGG4D2_MRGFile art_merged.pool.root
