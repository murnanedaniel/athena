#!/bin/sh
#
# art-description: RecoTrf
# art-type: grid
# art-include: 21.0/Athena
# art-include: 21.0-TrigMC/Athena
# art-include: master/Athena
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-athena-mt: 8

Reco_tf.py \
--AMI=q220 \
--athenaopts='--nprocs=2' \
--outputAODFile=myAOD.pool.root \
--maxEvents=100 \
--outputESDFile=myESD.pool.root --outputHISTFile=myHIST.root --imf False

echo "art-result: $?"

ArtPackage=$1
ArtJobName=$2
art.py compare grid --entries 30 ${ArtPackage} ${ArtJobName}
echo "art-result: $?"
