#!/bin/sh
#
# art-description: RecoTrf
# art-type: grid
# art-include: master/Athena
# art-athena-mt: 8                                                                                                                                    

Reco_tf.py \
--AMI=f628 \
--athenaopts="--threads=8" \
--maxEvents=300 \
--ignoreErrors=False \
--inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/expressstream_input_data18/data18_13TeV.00357750.express_express.merge.RAW._lb0093._SFO-ALL._0001.1 \
--conditionsTag='CONDBR2-BLKPA-RUN2-02' \
--geometryVersion='ATLAS-R2-2016-01-00-01' \
--outputESDFile=myESD_express_0.pool.root --outputAODFile=myAOD_express_0.AOD.pool.root --outputHISTFile=myMergedMonitoring_express_0.root --imf False

echo "art-result: $?"

ArtPackage=$1
ArtJobName=$2
art.py compare grid --entries 30 ${ArtPackage} ${ArtJobName}
echo "art-result: $?"
