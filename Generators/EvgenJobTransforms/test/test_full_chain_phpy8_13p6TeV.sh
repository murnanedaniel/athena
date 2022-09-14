#!/bin/sh

# art-include: master/Athena
# art-description: Pythia event generation -- ttbar
# art-architecture: '#x86_64-intel'
# art-type: grid

Gen_tf.py --ecmEnergy=13600. --maxEvents=10 --randomSeed=123456 --outputEVNTFile=test_ttbar_EVNT.pool.root --jobConfig=950070

echo "art-result: $? Gen_tf"

EVNTMerge_tf.py --inputEVNTFile=test_ttbar_EVNT.pool.root --maxEvent=10 --skipEvents=0 --outputEVNT_MRGFile=EVNT.MERGE_pool.root

echo "art-result: $? ENVTMerge_tf"

Sim_tf.py --AMIConfig s3822 --inputEVNTFile test_ttbar_EVNT.pool.root --maxEvents 10 --runNumber 950070 --firstEvent 759001 --randomSeed 760 --outputHITSFile HITS.28723600._001356.pool.root.1 --physicsList FTFP_BERT_ATL_VALIDATION --truthStrategy MC15aPlus --conditionsTag 'OFLCOND-MC21-SDR-RUN3-05' --geometryVersion 'ATLAS-R3S-2021-02-00-00_VALIDATION' --preExec 'simFlags.CalibrationRun.set_Off();' --postExec 'conddb.addOverride("/Indet/Beampos","IndetBeampos-RunDep-MC21-BestKnowledge-002");' --preInclude 'SimulationJobOptions/preInclude.BeamPipeKill.py' 'SimulationJobOptions/preInclude.G4Optimizations.py' 'SimulationJobOptions/preInclude.ExtraParticles.py' 'SimulationJobOptions/preInclude.G4ExtraProcesses.py' --postInclude 'RecJobTransforms/UseFrontier.py' 'Campaigns/postInclude.MC21BirksConstantTune.py'


echo "art-result: $? Sim_tf"

export ATHENA_CORE_NUMBER=8
Reco_tf.py --inputHITSFile=HITS.28723600._001356.pool.root.1 --athenaMPEventsBeforeFork=1 --athenaopts "HITtoRDO:--nprocs=$ATHENA_CORE_NUMBER" "HITtoRDO:--threads=0" --deleteIntermediateOutputfiles=True --maxEvents=10 --multithreaded=True --postExec "all:conddb.addOverride(\"/Indet/Beampos\",\"IndetBeampos-RunDep-MC21-BestKnowledge-002\")" "all:from IOVDbSvc.CondDB import conddb; conddb.addOverride(\"/TRT/Calib/PID_NN\", \"TRTCalibPID_NN_v2\")" "all:conddb.addOverride(\"/PIXEL/ChargeCalibration\",\"PixelChargeCalibration-SIM-MC16-000-12\")" --postInclude "default:PyJobTransforms/UseFrontier.py" --preExec "HITtoRDO:userRunLumiOverride={\"run\":330000, \"startmu\":25.0,\"endmu\":52.0,\"stepmu\":1,\"startlb\":1,\"timestamp\":1560000000};from Digitization.DigitizationFlags import digitizationFlags;digitizationFlags.doPixelPlanarRadiationDamage=True" "all:from AthenaConfiguration.AllConfigFlags import ConfigFlags; ConfigFlags.Trigger.AODEDMSet='AODFULL'" --preInclude "all:Campaigns/MC21a.py" "HITtoRDO:Digitization/ForceUseOfPileUpTools.py,SimulationJobOptions/preInlcude.PileUpBunchTrainsMC16c_2017_Config1.py,RunDependentSimData/configLumi_muRange.py" --skipEvents=0 --autoConfiguration=everything --valid=True --conditionsTag "default:OFLCOND-MC21-SDR-RUN3-04" --geometryVersion="default:ATLAS-R3S-2021-02-00-00" --runNumber=950070 --digiSeedOffset1=10 --digiSeedOffset2=10 --digiSteeringConf='StandardSignalOnlyTruth'  --steering "doRDO_TRIG" "doTRIGtoALL" --inputHighPtMinbiasHitsFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CampaignInputs/mc21/mc21_13p6TeV.800831.Py8EG_minbias_inelastic_highjetphotonlepton.merge.HITS.e8453_e8455_s3879_s3880 --inputLowPtMinbiasHitsFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CampaignInputs/mc21/mc21_13p6TeV.900311.Epos_minbias_inelastic_lowjetphoton.merge.HITS.e8453_s3879_s3880 --numberOfHighPtMinBias=0.103 --numberOfLowPtMinBias=52.397 --pileupFinalBunch=6 --pileupInitialBunch=-32 --outputAODFile=AOD.28723600._001356.pool.root.1 --jobNumber=10 


echo "art-result: $? Reco_tf"
