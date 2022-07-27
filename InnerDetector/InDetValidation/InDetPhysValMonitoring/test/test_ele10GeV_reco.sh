#!/bin/bash
# art-description: art job for InDetPhysValMonitoring, Single ele 10GeV
# art-type: grid
# art-input: user.keli:user.keli.mc16_13TeV.422029.ParticleGun_single_ele_Pt10GeV_Rel22073
# art-input-nfiles: 10
# art-cores: 4
# art-memory: 4096
# art-include: master/Athena
# art-include: 22.0/Athena
# art-output: physval*.root
# art-output: *.xml 
# art-output: art_core_0
# art-output: dcube*

#RDO is made at rel 22.0.73
#reference plots are made at rel 22.0.73

set -x

echo "ArtProcess: $ArtProcess"
lastref_dir=last_results
script="`basename \"$0\"`"
success_run=0

case $ArtProcess in
  "start")
    echo "Starting"
    echo "List of files = " ${ArtInFile}
    ;;
  "end")
    echo "Ending"
    if [ ${success_run} -eq 0 ]  ;then
      echo "download latest result"
      art.py download --user=artprod --dst="$lastref_dir" "$ArtPackage" "$ArtJobName"
      ls -la "$lastref_dir"
      echo "Merging physval.root"
      hadd  physval.root art_core_*/physval.ntuple.root
      echo "postprocess"
      postProcessIDPVMHistos physval.root

      dcubeXml="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/InDetPhysValMonitoring/dcube/config/IDPVMPlots_R22_GSF.xml"
      dcubeRef="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/InDetPhysValMonitoring/ReferenceHistograms/physval_ele10GeV_reco_r22.root"
      echo "compare with R22.0.73"
      $ATLAS_LOCAL_ROOT/dcube/current/DCubeClient/python/dcube.py \
	   -p -x dcube \
	   -c ${dcubeXml} \
	   -r ${dcubeRef} \
	   physval.root
      echo "art-result: $? dcube"

      echo "compare with last build"
      $ATLAS_LOCAL_ROOT/dcube/current/DCubeClient/python/dcube.py \
	   -p -x dcube_last \
	   -c ${dcubeXml} \
	   -r last_results/physval.root \
	   physval.root
      echo "art-result: $? dcube_last"
    else
      echo "reco failed"
    fi
    ;;
  *)
    echo "Test $ArtProcess"
    mkdir "art_core_${ArtProcess}"
    cd "art_core_${ArtProcess}"
    IFS=',' read -r -a file <<< "${ArtInFile}"
    file=${file[${ArtProcess}]}
    x="../$file"
    echo "Unsetting ATHENA_NUM_PROC=${ATHENA_NUM_PROC} and ATHENA_PROC_NUMBER=${ATHENA_PROC_NUMBER}"
    unset  ATHENA_NUM_PROC
    unset  ATHENA_PROC_NUMBER

    Reco_tf.py \
      --inputRDOFile $x \
      --outputNTUP_PHYSVALFile physval.ntuple.root \
      --outputAODFile   physval.AOD.root \
      --conditionsTag   'OFLCOND-MC16-SDR-RUN2-08' \
      --steering        doRAWtoALL \
      --checkEventCount False \
      --ignoreErrors    True \
      --maxEvents       -1 \
      --skipEvents      0 \
      --valid           True \
      --validationFlags doInDet \
      --autoConfiguration everything \
      --postExec 'condSeq.TileSamplingFractionCondAlg.G4Version = -1' \
      --preExec 'from InDetRecExample.InDetJobProperties import InDetFlags; \
      InDetFlags.doSlimming.set_Value_and_Lock(False); rec.doTrigger.set_Value_and_Lock(False); \
      from InDetPhysValMonitoring.InDetPhysValJobProperties import InDetPhysValFlags; \
      InDetPhysValFlags.doValidateTightPrimaryTracks.set_Value_and_Lock(True); \
      InDetPhysValFlags.doValidateTracksInJets.set_Value_and_Lock(False); \
      InDetPhysValFlags.doValidateGSFTracks.set_Value_and_Lock(True); \
      InDetPhysValFlags.doExpertOutput.set_Value_and_Lock(True); \
      InDetPhysValFlags.doHitLevelPlots.set_Value_and_Lock(True); \
      rec.doDumpProperties=True; rec.doCalo=True; rec.doEgamma=True; \
      rec.doForwardDet=False; rec.doInDet=True; rec.doJetMissingETTag=True; \
      rec.doLArg=True; rec.doLucid=True; rec.doMuon=True; rec.doMuonCombined=True; \
      rec.doSemiDetailedPerfMon=True; rec.doTau=True; rec.doTile=True; \
      from ParticleBuilderOptions.AODFlags import AODFlags; \
      AODFlags.ThinGeantTruth.set_Value_and_Lock(False);  \
      AODFlags.ThinNegativeEnergyCaloClusters.set_Value_and_Lock(False); \
      AODFlags.ThinNegativeEnergyNeutralPFOs.set_Value_and_Lock(False);\
      AODFlags.ThinInDetForwardTrackParticles.set_Value_and_Lock(False) '
    rec_tf_exit_code=$?
    echo "art-result: $rec_tf_exit_code reco_${file}"
    if [ $rec_tf_exit_code -ne 0 ]  ;then
      success_run=$rec_tf_exit_code
    fi
    ls -lR
    ;;
esac
