#!/bin/sh

# art-include: master--HepMC3/Athena
# art-include: master/Athena
# art-description: Pythia event generation -- min_bias
# art-type: grid

Gen_tf.py --ecmEnergy=13600. --maxEvents=10000 --randomSeed=123456 --outputEVNTFile=EVNT.root --jobConfig=421113 --AMIConfig=e8462

echo "art-result: $? Gen_tf"

EVNTMerge_tf.py --inputEVNTFile=EVNT.root --maxEvent=10000 --skipEvents=0 --outputEVNT_MRGFile=EVNT.MERGE_pool.root --AMIConfig=e8455

echo "art-result: $? ENVTMerge_tf"



Sim_tf.py --AMIConfig=s3879 --maxEvents=10000 --inputEVNTFile=EVNT.MERGE_pool.root --outputHITSFile=HITS_421113.pool.root #--jobNumber=113

echo "art-result: $? Sim_tf"

HITSMerge_tf.py  --AMIConfig=s3875 --inputHITSFile=HITS_421113.pool.root --outputHITS_MRGFile=HITS_Mrg_421113_13p6.pool.root --maxEvents=10000 --skipEvents=0

echo "art-results: $? HITSMerge_tf"



Reco_tf.py --AMIConfig=r13829 --inputHITSFile=HITS_Mrg_421113_13p6.pool.root --maxEvents=10000 --jobNumber=113 --outputAODFile=AOD_421111_13p6.pool.root --inputRDO_BKGFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/mc21_13p6TeV/RDO.29616548._005230.pool.root.1

echo "art-result: $? Reco_tf"

AODMerge_tf.py  --AMIConfig=r13831 --inputAODFile=AOD_421113.pool.root --outputAOD_MRGFile=AOD_merge_421113.pool.root --skipEvents=0

echo "art-result: $? AODMerge_tf"
