#!/bin/sh
#
# art-description: ATLFASTIIF_G4MS test with pile-up profile
# art-type: grid
# art-include: master/Athena
# art-architecture: '#x86_64-intel'

maxevent=25
inputfile="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/ISF_Validation/mc12_valid.110401.PowhegPythia_P2012_ttbar_nonallhad.evgen.EVNT.e3099.01517252._000001.pool.root.1"
HighPtMinbiasHitsFiles="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.merge.HITS.e4981_s3087_s3089/*"
LowPtMinbiasHitsFiles="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.merge.HITS.e4981_s3087_s3089/*"
rdoFile="RDO.pool.root"
aodFile="AOD.pool.root"
ntupFile="physval_g4ms.root"

FastChain_tf.py \
    --CA \
    --simulator ATLFASTIIF_G4MS \
    --useISF True \
    --digiSteeringConf "StandardSignalOnlyTruth" \
    --randomSeed 123 \
    --enableLooperKiller True \
    --physicsList 'FTFP_BERT_ATL' \
    --jobNumber 1 \
    --digiSeedOffset1 '1' \
    --digiSeedOffset2 '2' \
    --inputEVNTFile ${inputfile} \
    --outputRDOFile ${rdoFile} \
    --maxEvents ${maxevent} \
    --skipEvents 0 \
    --geometryVersion default:ATLAS-R2-2016-01-00-01 \
    --conditionsTag default:OFLCOND-MC16-SDR-RUN2-09 \
    --preSimExec 'from TrkDetDescrSvc.TrkDetDescrJobProperties import TrkDetFlags;TrkDetFlags.TRT_BuildStrawLayers=True;' \
    --preInclude 'Campaigns.MC16a' \
    --postInclude='PyJobTransforms.UseFrontier' \
    --postExec 'from IOVDbSvc.IOVDbSvcConfig import addOverride;cfg.merge(addOverride(ConfigFlags,"/TILE/OFL02/CALIB/SFR","TileOfl02CalibSfr-SIM-05"))' \
    --inputHighPtMinbiasHitsFile ${HighPtMinbiasHitsFiles} \
    --inputLowPtMinbiasHitsFile ${LowPtMinbiasHitsFiles} \
    --pileupFinalBunch '6' \
    --numberOfHighPtMinBias '0.116075313' \
    --numberOfLowPtMinBias '44.3839246425' \
    --numberOfCavernBkg 0 \
    --imf False
rc1=$?
echo  "art-result: ${rc1} EVNTtoRDO"

# RDO -> AOD and AOD -> NTUP stages
rc1_1=999 
rc1_2=999 
if [ ${rc1} -eq 0 ]
then
    echo "Running Reco_tf.py:  RDO to AOD"
    Reco_tf.py --inputRDOFile ${rdoFile} --maxEvents '-1' \
               --skipEvents '0' --conditionsTag 'default:OFLCOND-MC16-SDR-RUN2-09' \
               --geometryVersion 'default:ATLAS-R2-2016-01-00-01' \
               --outputAODFile ${aodFile} \
               --steering 'doRDO_TRIG' \
               --athenaopts "all:--threads=1" \
               --imf False
     rc1_1=$?
     if [ ${rc1_1} -eq 0 ]
     then
	 echo "Running Reco_tf.py:  AOD to NTUP"
         # NTUP prod.
         Reco_tf.py --inputAODFile ${aodFile} --maxEvents '-1' \
                    --outputNTUP_PHYSVALFile ${ntupFile} \
                    --ignoreErrors True \
                    --validationFlags 'doInDet' \
                    --valid 'True'
         rc1_2=$?
     fi
fi
echo "art-result: ${rc1_1} RDOtoAOD"
echo "art-result: ${rc1_2} AODtoNTUP"

rc2=999
if [ ${rc1_1} -eq 0 ]
then
    art.py compare grid --entries 10 ${ArtPackage} ${ArtJobName} --mode=summary
    rc2=$?
fi
echo "art-result: ${rc2} regression"
