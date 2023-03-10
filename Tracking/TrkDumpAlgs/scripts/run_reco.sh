Reco_tf.py \
    --digiSteeringConf 'StandardInTimeOnlyTruth' \
    --conditionsTag 'all:OFLCOND-MC15c-SDR-14-05' \
    --imf 'all:False' \
    --pileupFinalBunch '6' --numberOfHighPtMinBias '0.725172' --numberOfLowPtMinBias '209.2692' \
    --geometryVersion 'all:ATLAS-P2-ITK-23-00-03' \
    --preInclude 'all:InDetSLHC_Example/preInclude.NoTRT_NoBCM_NoDBM.py,InDetSLHC_Example/preInclude.SLHC_Setup.py,InDetSLHC_Example/preInclude.SLHC_Setup_Strip_GMX.py,InDetSLHC_Example/preInclude.SLHC_Calorimeter_mu200.py' 'HITtoRDO:InDetSLHC_Example/preInclude.SLHC.py,Digitization/ForceUseOfPileUpTools.py,SimulationJobOptions/preInclude.PileUpBunchTrains2012Config1_DigitConfig.py,RunDependentSimData/configLumi_muRange.py' 'default:InDetSLHC_Example/preInclude.SLHC.NoTRT_NoBCM_NoDBM.Reco.py,InDetSLHC_Example/SLHC_Setup_Reco_TrackingGeometry_GMX.py' \
    --DataRunNumber '242020' \
    --postInclude 'all:InDetSLHC_Example/postInclude.SLHC_Setup_ITK.py,LArROD/LArSuperCellEnable.py' 'HITtoRDO:InDetSLHC_Example/postInclude.SLHC_Digitization.py' 'RAWtoESD:InDetSLHC_Example/postInclude.AnalogueClustering.py' 'RAWtoESD:TrkDumpAlgs/postInclude.DumpObjects.py' \
    --preExec 'all:from AthenaCommon.GlobalFlags import globalflags; globalflags.DataSource.set_Value_and_Lock("geant4"); from InDetSLHC_Example.SLHC_JobProperties import SLHC_Flags; SLHC_Flags.doGMX.set_Value_and_Lock(True); from InDetRecExample.InDetJobProperties import InDetFlags; InDetFlags.keepAdditionalHitsOnTrackParticle.set_Value_and_Lock(True)' 'HITtoRDO:from Digitization.DigitizationFlags import digitizationFlags; digitizationFlags.doInDetNoise.set_Value_and_Lock(False); digitizationFlags.overrideMetadata+=["SimLayout","PhysicsList"]; userRunLumiOverride={"run":242020, "startmu":190.0, "endmu":210.0, "stepmu":1.0, "startlb":1, "timestamp":1412020000};digitizationFlags.experimentalDigi+=["SimpleMerge"]' 'RAWtoESD:from JetRec.JetRecFlags import jetFlags, JetContentDetail; jetFlags.detailLevel.set_Value_and_Lock(JetContentDetail.Full); from AthenaCommon.DetFlags import DetFlags;DetFlags.geometry.HGTD_setOn();DetFlags.HGTD_setOn()' 'ESDtoAOD:from JetRec.JetRecFlags import jetFlags, JetContentDetail; jetFlags.detailLevel.set_Value_and_Lock(JetContentDetail.Full)' \
    --HGTDOn 'True' \
    --postExec 'all:ServiceMgr.PixelLorentzAngleSvc.ITkL03D = True' 'HITtoRDO:CfgMgr.MessageSvc().setError+=["HepMcParticleLink"];conddb.addMarkup("MDT/Ofl/CABLING/MEZZANINE_SCHEMA","<forceRunNumber>232550</forceRunNumber>");conddb.addMarkup("MDT/Ofl/CABLING/MAP_SCHEMA","<forceRunNumber>232550</forceRunNumber>")' 'RAWtoESD:ToolSvc.InDetSCT_ClusteringTool.useRowInformation=True' 'ESDtoAOD:StreamAOD.ItemList+=["CaloCellContainer#AllCalo"]' 'RAWtoESD:dumpObjects.csvFile=True' 'RAWtoESD:dumpObjects.rootFile=True; rootfile();'\
    --inputRDOFile=/sps/l2it/stark/NewSamples/RDO/mc15_14TeV.600012.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.RDO.e8185_s3770_s3773_r13618_tid31548381_00/RDO.31548381._000001.pool.root.1 \
    --outputESDFile=ESD.31548381._000001.root \


