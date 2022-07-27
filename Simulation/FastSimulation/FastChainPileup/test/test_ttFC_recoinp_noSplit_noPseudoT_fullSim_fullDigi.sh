#!/usr/bin/env bash

# art-description: test job ttFC_fullSim_fullDigi + ttFC_reco_noSplit_noPseudoT_fullSim_fullDigi (was for grid, master/Athena)
# art-output: config.txt
# art-output: RAWtoESD_config.txt
# art-output: *.root
# art-output: dcube-id
# art-html: dcube-id


inputRefDir="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/FastChainPileup/DCube-refs/${AtlasBuildBranch}/test_ttFC_recoinp_noSplit_noPseudoT_fullSim_fullDigi"
inputXmlDir="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/FastChainPileup/DCube-configs/${AtlasBuildBranch}"
dcubeXmlID="${inputXmlDir}/physval-noSplit_noPseudoT_fullSim_fullDigi.xml"
dcubeRefID="${inputRefDir}/physval-noSplit_noPseudoT_fullSim_fullDigi.root"
aodFile="AOD_noSplit_noPseudoT_fullSim_fullDigi.pool.root"
ntupFile="physval-noSplit_noPseudoT_fullSim_fullDigi.root"

# RDO input from 21.3/Apr 27 nightly test_ttFC_reco_noSplit_noPseudoT_fullSim_fullDigi.sh
FastChain_tf.py --maxEvents 500 \
    --skipEvents 0 \
    --geometryVersion ATLAS-R2-2015-03-01-00 \
    --conditionsTag OFLCOND-MC16-SDR-RUN2-09 \
    --inputRDOFile '/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/FastChainPileup/RecoInputRefs/RDO_pileup_fullsim_fulldigi.pool.root' \
    --outputAODFile ${aodFile} \
    --preExec "all:rec.doTrigger.set_Value_and_Lock(False)" \
    --postExec 'RAWtoESD:from AthenaCommon.ConfigurationShelve import saveToAscii;saveToAscii("RAWtoESD_config.txt")' \
    --imf False
rc2=$?

rc3=999
rc4=999
rc5=999
if [ ${rc2} -eq 0 ]
then
    # NTUP prod.
    Reco_tf.py --inputAODFile ${aodFile} --maxEvents '-1' \
           --outputNTUP_PHYSVALFile ${ntupFile} \
           --ignoreErrors True \
           --validationFlags 'doInDet' \
           --valid 'True'
    rc3=$?

    # Regression test
    ArtPackage=$1
    ArtJobName=$2
    art.py compare grid --entries 10 ${ArtPackage} ${ArtJobName} --mode=summary
    rc4=$?

    if [ ${rc3} -eq 0 ]
    then
        # Histogram comparison with DCube
        $ATLAS_LOCAL_ROOT/dcube/current/DCubeClient/python/dcube.py \
        -p -x dcube-id \
        -c ${dcubeXmlID} -r ${dcubeRefID} ${ntupFile}
        rc5=$?
    fi
fi

echo  "art-result: ${rc2} RDOtoAOD"
echo  "art-result: ${rc3} AODtoNTUP"
echo  "art-result: ${rc4} regression"
echo  "art-result: ${rc5} dcubeID"
