#!/bin/sh

# art-description: Run FastChain with Simulation (ATLFAST3F_G4MS) and MC+MC Overlay in one job without reco for MC21a (RUN3), ttbar
# art-type: grid
# art-include: 22.0/Athena
# art-output: *.root
# art-output: config.txt
# art-output: RAWtoESD_config.txt
# art-output: ESDtoAOD_config.txt
# art-architecture: '#x86_64-intel'

events=25
EVNT_File="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.evgen.EVNT.e4993.EVNT.08166201._000012.pool.root.1"
RDO_BKG_File="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayTests/PresampledPileUp/22.0/Run3/v3/mc21a_presampling.RDO.pool.root"
RDO_File="MC_plus_MC.RDO.pool.root"
AOD_File="MC_plus_MC.AOD.pool.root"
NTUP_File="MC_plus_MC.NTUP.pool.root"


FastChain_tf.py \
  --simulator ATLFAST3F_G4MS \
  --steering doFCwOverlay \
  --physicsList FTFP_BERT_ATL \
  --useISF True \
  --randomSeed 123 \
  --inputEVNTFile ${EVNT_File} \
  --inputRDO_BKGFile ${RDO_BKG_File} \
  --outputRDOFile ${RDO_File} \
  --maxEvents ${events} \
  --skipEvents 0 \
  --digiSeedOffset1 511 \
  --digiSeedOffset2 727 \
  --conditionsTag 'OFLCOND-MC21-SDR-RUN3-05' \
  --geometryVersion 'ATLAS-R3S-2021-02-00-00' \
  --postInclude 'default:PyJobTransforms/UseFrontier.py' \
  --preInclude 'all:Campaigns/MC21Simulation.py,SimulationJobOptions/preInclude.FrozenShowersFCalOnly.py,Campaigns/MC21a.py' \
  --postExec 'from AthenaCommon.ConfigurationShelve import saveToAscii;saveToAscii("config.txt")' \
  --DataRunNumber '330000' \
  --imf False

rc=$?
echo "art-result: ${rc} EVNTtoRDOwOverlay"


rc2=999
rc3=999
rc4=999
if [ ${rc} -eq 0 ]
then
    # Reconstruction
    Reco_tf.py --inputRDOFile ${RDO_File} --maxEvents '-1' \
               --autoConfiguration=everything \
               --outputAODFile ${AOD_File} \
               --steering 'doRDO_TRIG' \
               --athenaopts "all:--threads=1" \
               --postExec 'RAWtoESD:from AthenaCommon.ConfigurationShelve import saveToAscii;saveToAscii("RAWtoESD_config.txt")' 'ESDtoAOD:from AthenaCommon.ConfigurationShelve import saveToAscii;saveToAscii("ESDtoAOD_config.txt")' \
               --imf False

     rc2=$?
     if [ ${rc2} -eq 0 ]
     then
         # NTUP prod.
         Reco_tf.py --inputAODFile ${AOD_File} --maxEvents '-1' \
                    --outputNTUP_PHYSVALFile ${NTUP_File} \
                    --ignoreErrors True \
                    --validationFlags 'doInDet' \
                    --valid 'True'
         rc3=$?

         # regression test
         art.py compare grid --entries 10 ${ArtPackage} ${ArtJobName} --mode=summary
         rc4=$?
     fi
fi


echo  "art-result: ${rc2} RDOtoAOD"
echo  "art-result: ${rc3} AODtoNTUP"
echo  "art-result: ${rc4} regression"
