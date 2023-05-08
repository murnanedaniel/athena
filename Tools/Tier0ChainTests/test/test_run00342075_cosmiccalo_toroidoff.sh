#!/bin/sh
#
# art-description: RecoTrf
# art-type: grid
# art-include: 21.0/Athena
# art-include: 21.0-mc16a/Athena
# art-include: 21.0-mc16d/Athena
# art-include: 21.0-TrigMC/Athena
# art-include: master/Athena
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-output: log.*

Reco_tf.py --conditionsTag all:CONDBR2-BLKPA-2017-14 --ignoreErrors 'False' --autoConfiguration='everything' --maxEvents '500' --AMITag 'f908' --preExec  'r2a:from InDetRecExample.InDetJobProperties import InDetFlags; InDetFlags.useDynamicAlignFolders.set_Value_and_Lock(True);TriggerFlags.AODEDMSet="AODFULL";rec.doTrigger=False;' --geometryVersion all:ATLAS-R2-2016-01-00-01 --steering='doRAWtoALL' --imf False --outputHISTFile=myMergedMonitoring_CosmicCalo_0.root --outputESDFile=myESD_CosmicCalo_0.pool.root --outputAODFile=myAOD_CosmicCalo_0.AOD.pool.root --outputTAGFile=myTAG_CosmicCalo_0.root --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_cos.00342075.physics_CosmicCalo.merge.RAW._lb1545._SFO-ALL._0001.1,/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_cos.00342075.physics_CosmicCalo.merge.RAW._lb1546._SFO-ALL._0001.1,/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_cos.00342075.physics_CosmicCalo.merge.RAW._lb1547._SFO-ALL._0001.1
echo "art-result: $? Reco"

ArtRef=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/$1/TCT_21.0-mc16d_references/$2
art.py compare ref --entries 10 . $ArtRef
echo "art-result: $? Diff"
