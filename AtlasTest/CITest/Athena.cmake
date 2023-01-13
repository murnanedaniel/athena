# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#
# CI test definitions for the Athena project
# --> README.md before you modify this file
#

#################################################################################
# General
#################################################################################
atlas_add_citest( DuplicateClass
   SCRIPT python -c 'import ROOT'
   PROPERTIES FAIL_REGULAR_EXPRESSION "class .* is already in" )

#################################################################################
# Digitization/Simulation
#################################################################################

atlas_add_citest( FastChain
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/FastChain.sh )

atlas_add_citest( SimulationRun2AF3
   SCRIPT RunWorkflowTests_Run2.py --CI -s -w AF3 )

atlas_add_citest( SimulationRun4FullSim
   SCRIPT RunWorkflowTests_Run4.py --CI -s -w FullSim -e '--maxEvents 5' --no-output-checks
   LOG_IGNORE_PATTERN "WARNING FPE INVALID" )  # ignore FPEs from Geant4

atlas_add_citest( PileUpPresamplingRun2
   SCRIPT RunWorkflowTests_Run2.py --CI -p -w PileUpPresampling -e '--maxEvents 5' --no-output-checks )

atlas_add_citest( PileUpPresamplingRun3
   SCRIPT RunWorkflowTests_Run3.py --CI -p -w PileUpPresampling -e '--maxEvents 5' --no-output-checks )

atlas_add_citest( OverlayRun2MC
   SCRIPT RunWorkflowTests_Run2.py --CI -o -w MCOverlay )

atlas_add_citest( OverlayRun2Data
   SCRIPT RunWorkflowTests_Run2.py --CI -o -w DataOverlay )

atlas_add_citest( OverlayRun3MC
   SCRIPT RunWorkflowTests_Run3.py --CI -o -w MCOverlay )

atlas_add_citest( OverlayRun3MC_CAConfig
   SCRIPT RunWorkflowTests_Run3.py --CI -o -w MCOverlay -e '--CA True' )

#################################################################################
# Standard reconstruction workflows
#################################################################################

atlas_add_citest( RecoRun2Data
   SCRIPT RunWorkflowTests_Run2.py --CI -r -w DataReco -e '--maxEvents 25' )

atlas_add_citest( RecoRun2Data_CAConfig
   SCRIPT RunWorkflowTests_Run2.py --CI -r -w DataReco -e '--CA --maxEvents 25' --no-output-checks )

atlas_add_citest( RecoRun2Data_LegacyVsCA
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/RecoLegacyVsCA.sh RecoRun2Data q442
   DEPENDS_SUCCESS RecoRun2Data RecoRun2Data_CAConfig )

atlas_add_citest( RecoRun2Data_DAODPHYS
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/DAODPhys.sh PHYS ../RecoRun2Data/run_q442/myAOD.pool.root
   DEPENDS_SUCCESS RecoRun2Data )

atlas_add_citest( RecoRun2Data_DAODPHYSLite
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/DAODPhys.sh PHYSLITE ../RecoRun2Data/run_q442/myAOD.pool.root
   DEPENDS_SUCCESS RecoRun2Data )

atlas_add_citest( RecoRun2MC
   SCRIPT RunWorkflowTests_Run2.py --CI -r -w MCReco --threads 0 -e '--maxEvents 25' )

atlas_add_citest( RecoRun2MC_CAConfig
   SCRIPT RunWorkflowTests_Run2.py --CI -r -w MCReco -e '--CA --maxEvents 25' --no-output-checks )

atlas_add_citest( RecoRun2MC_LegacyVsCA
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/RecoLegacyVsCA.sh RecoRun2MC q443
   DEPENDS_SUCCESS RecoRun2MC RecoRun2MC_CAConfig )

atlas_add_citest( RecoRun2MC_PileUp
   SCRIPT RunWorkflowTests_Run2.py --CI -p -w MCPileUpReco -e '--maxEvents 5 --inputRDO_BKGFile=../../PileUpPresamplingRun2/run_d1730/myRDO.pool.root' --no-output-checks  # go two levels up as the test runs in a subfolder
   DEPENDS_SUCCESS PileUpPresamplingRun2 )

atlas_add_citest( RecoRun3Data
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a q449 --threads 8 -e '--maxEvents 100' --run-only
   PROPERTIES PROCESSORS 8 )

atlas_add_citest( RecoRun3Data_Checks
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a q449 --threads 8 -e '--maxEvents 100' --checks-only --output-path ../RecoRun3Data
   DEPENDS_SUCCESS RecoRun3Data )

atlas_add_citest( RecoRun3Data_CAConfig
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -e '--CA --maxEvents 100' --threads 8 --no-output-checks )

atlas_add_citest( RecoRun3Data_LegacyVsCA
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/RecoLegacyVsCA.sh RecoRun3Data q449
   DEPENDS_SUCCESS RecoRun3Data RecoRun3Data_CAConfig )

atlas_add_citest( RecoRun3Data_Bulk
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a f1262 --threads 8 -e '--skipEvents 100 --maxEvents 500 --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_13p6TeV.00431493.physics_Main.daq.RAW._lb0525._SFO-16._0001.data --postExec="all:FPEAuditor.NStacktracesOnFPE=20"' --run-only --no-output-checks
   PROPERTIES PROCESSORS 8 )

atlas_add_citest( RecoRun3Data_Bulk_Checks
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a f1262 --threads 8 -e '--skipEvents 100 --maxEvents 500 --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_13p6TeV.00431493.physics_Main.daq.RAW._lb0525._SFO-16._0001.data --postExec="all:FPEAuditor.NStacktracesOnFPE=20"' --checks-only --output-path ../RecoRun3Data_Bulk --no-output-checks
   DEPENDS_SUCCESS RecoRun3Data_Bulk )

atlas_add_citest( RecoRun3Data_Express
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a x686 -e '--maxEvents 25 --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/TCT_Run3/data22_13p6TeV.00428353.express_express.merge.RAW._lb0800._SFO-ALL._0001.1' --no-output-checks )

atlas_add_citest( RecoRun3Data_Cosmics
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a q450 -e '--maxEvents 25' --no-output-checks )

atlas_add_citest( RecoRun3Data_Calib
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w DataReco -a q451 -e '--maxEvents 25' --no-output-checks )

atlas_add_citest( RecoRun3MC
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/RecoRun3MC.sh )

atlas_add_citest( RecoRun3MC_CAConfig
   SCRIPT RunWorkflowTests_Run3.py --CI -r -w MCReco -e '--CA --maxEvents 25' --no-output-checks )

atlas_add_citest( RecoRun3MC_LegacyVsCA
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/RecoLegacyVsCA.sh RecoRun3MC q445
   DEPENDS_SUCCESS RecoRun3MC RecoRun3MC_CAConfig )

atlas_add_citest( RecoRun3MC_PileUp
   SCRIPT RunWorkflowTests_Run3.py --CI -p -w MCPileUpReco -e '--maxEvents 5 --inputRDO_BKGFile=../../PileUpPresamplingRun3/run_d1760/myRDO.pool.root' --no-output-checks  # go two levels up as the test runs in a subfolder
   DEPENDS_SUCCESS PileUpPresamplingRun3 )

atlas_add_citest( RecoRun4MC
   SCRIPT RunWorkflowTests_Run4.py --CI -r -w MCReco -e '--maxEvents 5 --inputHITSFile=../../SimulationRun4FullSim/run_s3761/myHITS.pool.root' --no-output-checks  # go two levels up as the test runs in a subfolder
   LOG_IGNORE_PATTERN "WARNING FPE" 
   DEPENDS_SUCCESS SimulationRun4FullSim )


atlas_add_citest( RecoRun4MC_DAODPHYS
  SCRIPT RunWorkflowTests_Run4.py --CI -d -e '--maxEvents 5 --inputAODFile=../../RecoRun4MC/run_q447/myAOD.pool.root' --no-output-checks  # go two levels up as the test runs in a subfolder
  LOG_IGNORE_PATTERN "WARNING FPE INVALID"
  DEPENDS_SUCCESS RecoRun4MC )

#################################################################################
# Data Quality
#################################################################################

atlas_add_citest( DataQuality_r21ESD
   SCRIPT Run3DQTestingDriver.py 'Input.Files=["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/q431/21.0/myESD.pool.root"]' Output.HISTFileName='DataQuality_r21ESD_HIST.root' DQ.Steering.doHLTMon=False --threads=1 )

atlas_add_citest( DataQuality_Run3MC
   SCRIPT Run3DQTestingDriver.py 'Input.Files=["../RecoRun3MC/run_q445/myAOD.pool.root"]' DQ.Environment=AOD --threads=1
   DEPENDS_SUCCESS RecoRun3MC )

atlas_add_citest( DataQuality_r21ESD_Postprocessing
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/DataQuality_r21ESD_Postprocessing.sh
   DEPENDS_SUCCESS DataQuality_r21ESD )

atlas_add_citest( DataQuality_Run3Data_AODtoHIST
   SCRIPT Reco_tf.py --AMI=q449 --inputAODFile="../RecoRun3Data/run_q449/myAOD.pool.root" --outputHISTFile=DataQuality_Run3Data_AODtoHIST.root  --athenaopts='--threads=1'
   DEPENDS_SUCCESS RecoRun3Data )


#################################################################################
# Special reconstruction
#################################################################################

atlas_add_citest( EgammaCAConfig
   SCRIPT Reco_tf.py --CA --steering doRAWtoALL --inputRDOFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CampaignInputs/mc20/RDO/mc20_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3681_r13145/100events.RDO.pool.root --preInclude egammaConfig.egammaOnlyFromRawFlags.egammaOnlyFromRaw --outputAODFile=AOD.pool.root --maxEvents=1 )

atlas_add_citest( Egamma
   SCRIPT ut_egammaARTJob_test.sh )

#################################################################################
# ACTS
#################################################################################

atlas_add_citest( ACTS_Propagation_ITk
   SCRIPT ActsITkTest.py )

atlas_add_citest( ACTS_Propagation_ID
   SCRIPT ActsExtrapolationAlgTest.py )

atlas_add_citest( ACTS_Workflow
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsWorkflow.sh 
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

atlas_add_citest( ACTS_ValidateClusters
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsValidateClusters.sh
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

atlas_add_citest( ACTS_ValidateSpacePoints
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsValidateSpacePoints.sh
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

atlas_add_citest( ACTS_ValidateSeeds
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsValidateSeeds.sh
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

atlas_add_citest( ACTS_ValidateOrthogonalSeeds 
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsValidateOrthogonalSeeds.sh
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

atlas_add_citest( ACTS_ActsPersistifyEDM 
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test/ActsPersistifyEDM.sh 
    LOG_IGNORE_PATTERN "WARNING FPE INVALID" )

#################################################################################
# Trigger
#################################################################################

atlas_add_citest( TriggerMC
   SCRIPT test_trigAna_RDOtoRDOTrig_v1Dev_build.py )

atlas_add_citest( TriggerData
   SCRIPT test_trig_data_v1Dev_build.py )

atlas_add_citest( TriggerDataCAConfig
   SCRIPT test_trig_data_newJO_build.py )

atlas_add_citest( Trigger_athenaHLT_v1Dev
   SCRIPT test_trigP1_v1Dev_decodeBS_build.py )

atlas_add_citest( Trigger_athenaHLT_v1PhysP1
   SCRIPT test_trigP1_v1PhysP1_build.py )

atlas_add_citest( Trigger_athenaHLT_v1Cosmic
   SCRIPT test_trigP1_v1Cosmic_build.py )
