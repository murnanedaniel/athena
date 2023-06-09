################################################################################
# Package: TrigEventAthenaPool
################################################################################

# Declare the package name:
atlas_subdir( TrigEventAthenaPool )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Database/AthenaPOOL/AthenaPoolUtilities
                          GaudiKernel
                          PRIVATE
                          AtlasTest/TestTools
                          Database/AthenaPOOL/AthenaPoolCnvSvc
                          Database/AthenaPOOL/RootConversions
                          Trigger/TrigEvent/TrigBphysicsEvent
                          Trigger/TrigEvent/TrigCaloEvent
                          Trigger/TrigEvent/TrigCaloEventTPCnv
                          Trigger/TrigEvent/TrigCombinedEvent
                          Trigger/TrigEvent/TrigCombinedEventTPCnv
                          Trigger/TrigEvent/TrigDecisionEvent
                          Trigger/TrigEvent/TrigDecisionEventTPCnv
                          Trigger/TrigEvent/TrigInDetEvent
                          Trigger/TrigEvent/TrigInDetEventTPCnv
                          Trigger/TrigEvent/TrigMissingEtEvent
                          Trigger/TrigEvent/TrigMissingEtEventTPCnv
                          Trigger/TrigEvent/TrigMonitoringEvent
                          Trigger/TrigEvent/TrigMonitoringEventTPCnv
                          Trigger/TrigEvent/TrigMuonEvent
                          Trigger/TrigEvent/TrigMuonEventTPCnv
                          Trigger/TrigEvent/TrigParticle
                          Trigger/TrigEvent/TrigParticleTPCnv
                          Trigger/TrigEvent/TrigTopoEvent
                          Trigger/TrigEvent/TrigTopoEventTPCnv
                          Trigger/TrigFTK/TrigFTK_RawData
                          Trigger/TrigFTK/TrigFTK_RawDataTPCnv )

# Component(s) in the package:
atlas_add_poolcnv_library( TrigEventAthenaPoolPoolCnv
                           src/*.cxx
                           FILES TrigParticle/TrigL2BphysContainer.h TrigParticle/TrigElectronContainer.h TrigParticle/TrigTauContainer.h TrigParticle/TrigTau.h TrigParticle/TrigL2BjetContainer.h TrigParticle/TrigEFBjetContainer.h TrigParticle/TrigPhotonContainer.h TrigParticle/TrigEFBphysContainer.h TrigDecisionEvent/TrigDecision.h TrigInDetEvent/TrigInDetTrackCollection.h TrigInDetEvent/TrigVertexCollection.h TrigInDetEvent/TrigTauTracksInfo.h TrigInDetEvent/TrigTauTracksInfoCollection.h TrigInDetEvent/TrigHisto1D.h TrigInDetEvent/TrigHisto2D.h TrigInDetEvent/TrigSpacePointCounts.h TrigInDetEvent/TrigSpacePointCountsCollection.h TrigInDetEvent/TrigTrtHitCounts.h TrigInDetEvent/TrigTrtHitCountsCollection.h TrigInDetEvent/TrigTrackCounts.h TrigInDetEvent/TrigTrackCountsCollection.h TrigInDetEvent/TrigVertexCounts.h TrigInDetEvent/TrigVertexCountsCollection.h TrigCombinedEvent/TrigComposite.h TrigCombinedEvent/TrigCompositeContainer.h TrigMuonEvent/MuonFeature.h TrigMuonEvent/MuonFeatureContainer.h TrigMuonEvent/MuonFeatureDetails.h TrigMuonEvent/MuonFeatureDetailsContainer.h TrigMuonEvent/CombinedMuonFeature.h TrigMuonEvent/CombinedMuonFeatureContainer.h TrigMuonEvent/IsoMuonFeature.h TrigMuonEvent/IsoMuonFeatureContainer.h TrigMuonEvent/TrigMuonClusterFeature.h TrigMuonEvent/TrigMuonClusterFeatureContainer.h TrigMuonEvent/TileMuFeature.h TrigMuonEvent/TileMuFeatureContainer.h TrigMuonEvent/TileTrackMuFeature.h TrigMuonEvent/TileTrackMuFeatureContainer.h TrigMuonEvent/TrigMuonEFContainer.h TrigMuonEvent/TrigMuonEFInfoContainer.h TrigMuonEvent/TrigMuonEFIsolationContainer.h TrigMuonEvent/TrigCombDiMuonContainer.h TrigBphysicsEvent/TrigJpsi.h TrigBphysicsEvent/TrigJpsiContainer.h TrigCaloEvent/TrigEMCluster.h TrigCaloEvent/TrigEMClusterContainer.h TrigCaloEvent/TrigTauCluster.h TrigCaloEvent/TrigTauClusterContainer.h TrigCaloEvent/TrigTauClusterDetails.h TrigCaloEvent/TrigTauClusterDetailsContainer.h TrigCaloEvent/TrigT2MbtsBits.h TrigCaloEvent/TrigT2MbtsBitsContainer.h TrigCaloEvent/TrigT2ZdcSignals.h TrigCaloEvent/TrigT2ZdcSignalsContainer.h TrigCaloEvent/TrigT2Jet.h TrigCaloEvent/TrigT2JetContainer.h TrigMissingEtEvent/TrigMissingET.h TrigMissingEtEvent/TrigMissingETContainer.h TrigCaloEvent/RingerRings.h TrigCaloEvent/RingerRingsContainer.h TrigCaloEvent/TrigRNNOutput.h TrigCaloEvent/TrigRNNOutputContainer.h TrigCaloEvent/TrigCaloCluster.h TrigCaloEvent/TrigCaloClusterContainer.h TrigTopoEvent/ElectronMuonTopoInfo.h TrigTopoEvent/ElectronMuonTopoInfoContainer.h TrigMonitoringEvent/TrigMonEvent.h TrigMonitoringEvent/TrigMonEventCollection.h TrigMonitoringEvent/TrigMonConfig.h TrigMonitoringEvent/TrigMonConfigCollection.h TrigFTK_RawData/FTK_RawTrackContainer.h
                           TYPES_WITH_NAMESPACE TrigDec::TrigDecision
                           LINK_LIBRARIES AthenaPoolUtilities GaudiKernel AthenaPoolCnvSvcLib RootConversions TrigBphysicsEvent TrigCaloEvent TrigCaloEventTPCnv TrigCombinedEvent TrigCombinedEventTPCnv TrigDecisionEvent TrigDecisionEventTPCnv TrigInDetEvent TrigInDetEventTPCnv TrigMissingEtEvent TrigMissingEtEventTPCnv TrigMonitoringEvent TrigMonitoringEventTPCnv TrigMuonEvent TrigMuonEventTPCnv TrigParticle TrigParticleTPCnv TrigTopoEvent TrigTopoEventTPCnv TrigFTK_RawData TrigFTK_RawDataTPCnv )

# Install files from the package:
atlas_install_headers( TrigEventAthenaPool )
atlas_install_joboptions( share/*.py )

# Set up (a) test(s) for the converter(s):
if( IS_DIRECTORY ${CMAKE_SOURCE_DIR}/Database/AthenaPOOL/AthenaPoolUtilities )
   set( AthenaPoolUtilitiesTest_DIR
      ${CMAKE_SOURCE_DIR}/Database/AthenaPOOL/AthenaPoolUtilities/cmake )
endif()
find_package( AthenaPoolUtilitiesTest )

if( ATHENAPOOLUTILITIESTEST_FOUND )
  run_tpcnv_test( TrigEventAthenaPool_15.6.1    ESD-15.6.1 )
  run_tpcnv_test( TrigEventAthenaPool_15.6.13.2 AOD-15.6.13.2-full )
  run_tpcnv_test( TrigEventAthenaPool_16.0.3.3  AOD-16.0.3.3-full )
  run_tpcnv_test( TrigEventAthenaPool_17.0.6    AOD-17.0.6-full )
  run_tpcnv_test( TrigEventAthenaPool_18.0.0    AOD-18.0.0-full )
  run_tpcnv_test( TrigEventAthenaPool_20.1.7.2  ESD-20.1.7.2 )
else()
   message( WARNING "Couldn't find AthenaPoolUtilitiesTest. No test(s) set up." )
endif()   
