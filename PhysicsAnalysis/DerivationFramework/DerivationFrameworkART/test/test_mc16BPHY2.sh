#!/bin/sh

# art-description: DAOD building BPHY2 mc16
# art-type: grid
# art-output: *.pool.root

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/AOD.11866988._000378.pool.root.1 --outputDAODFile art.pool.root --reductionConf BPHY2 --maxEvents 5000

DAODMerge_tf.py --maxEvents 5 --inputDAOD_BPHY2File DAOD_BPHY2.art.pool.root --outputDAOD_BPHY2_MRGFile art_merged.pool.root
