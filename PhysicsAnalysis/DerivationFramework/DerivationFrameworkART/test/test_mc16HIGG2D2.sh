#!/bin/sh

# art-description: DAOD building HIGG2D2 mc16
# art-type: grid
# art-output: *.pool.root

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/AOD.11866988._000378.pool.root.1 --outputDAODFile art.pool.root --reductionConf HIGG2D2 --maxEvents 5000

DAODMerge_tf.py --maxEvents 5 --inputDAOD_HIGG2D2File DAOD_HIGG2D2.art.pool.root --outputDAOD_HIGG2D2_MRGFile art_merged.pool.root
