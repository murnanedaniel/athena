# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#Content included in addition to the Smart Slimming Content

ExtraContentElectrons=[
    "Electrons.Loose",
    "Electrons.Medium",
    "Electrons.Tight"
    ]

# only if DoCellReweighting is ON
ExtraContentReweightedElectrons = ["NewSwElectrons.trackParticleLinks.pt.eta.phi.m.caloClusterLinks.author.OQ.ethad1.ethad.f1.f3.f3core.e233.e237.e277.weta1.weta2.e2tsts1.fracs1.wtots1.emins1.emaxs1.etcone20.ptcone30.deltaEta1.deltaPhi1.deltaPhi2.deltaPhiRescaled2.deltaPhiFromLastMeasurement.Loose.Medium.Tight.DFCommonElectronsLHVeryLoose.DFCommonElectronsLHLoose.DFCommonElectronsLHLooseBL.DFCommonElectronsLHMedium.DFCommonElectronsLHTight.ptcone20.ptcone30.ptcone40.ptvarcone20.ptvarcone30.ptvarcone40.topoetcone20.topoetcone30.topoetcone40.charge.Reta.Rphi.Eratio.Rhad.Rhad1.DeltaE.topoetcone20ptCorrection.topoetcone30ptCorrection.topoetcone40ptCorrection.etcone20ptCorrection.etcone30ptCorrection.etcone40ptCorrection.ambiguityLink.truthParticleLink.truthOrigin.truthType.truthPdgId.firstEgMotherTruthType.firstEgMotherTruthOrigin.firstEgMotherTruthParticleLink.firstEgMotherPdgId.lastEgMotherTruthType.lastEgMotherTruthOrigin.lastEgMotherTruthParticleLink.lastEgMotherPdgId.ambiguityType.DFCommonAddAmbiguity"]
# might need to add extra variables for Min/Max variations... but not for the moment

ExtraElectronsTruth=[
    "Electrons.truthOrigin",
    "Electrons.truthType",
    "Electrons.truthParticleLink"]

ExtraContentMuons=[
    "Muons.ptcone20",
    "Muons.ptcone30",
    "Muons.ptcone40",
    "Muons.etcone20",
    "Muons.etcone30",
    "Muons.etcone40"]

ExtraMuonsTruth=[
    ]

ExtraContentPhotons=[
]

ExtraContentPrimaryVertices=["PrimaryVertices.x.y.sumPt2"]

ExtraPhotonsTruth=[
]

ExtraContentGSFConversionVertices=[
        "GSFConversionVertices.x",
        "GSFConversionVertices.y",
        "GSFConversionVertices.z",
        "GSFConversionVertices.px",
        "GSFConversionVertices.py",
        "GSFConversionVertices.pz",
        "GSFConversionVertices.pt1",
        "GSFConversionVertices.pt2",
        "GSFConversionVertices.etaAtCalo",
        "GSFConversionVertices.phiAtCalo",
        "GSFConversionVertices.trackParticleLinks"
        ]

ExtraContentHLTPhotons=[
        "HLT_xAOD__PhotonContainer_egamma_Photons.e.pt.m.author.Rhad.Rhad1.e277.Reta.Rphi.weta2.f1.fracs1.wtots1.weta1.DeltaE.Eratio.caloClusterLinks",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex.calE.calEta.calPhi.calM.e_sampl.eta_sampl.etaCalo.phiCalo.ETACALOFRAME.PHICALOFRAME"
]

from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import GainDecorator, getGainDecorations
GainDecoratorTool = GainDecorator()
ExtraContentPhotons.extend( getGainDecorations(GainDecoratorTool) )
ExtraContentElectrons.extend( getGainDecorations(GainDecoratorTool) )

ExtraContentAll=ExtraContentElectrons+ExtraContentMuons+ExtraContentPhotons+ExtraContentGSFConversionVertices+ExtraContentHLTPhotons+ExtraContentPrimaryVertices
ExtraContentAllTruth=ExtraElectronsTruth+ExtraMuonsTruth+ExtraPhotonsTruth

ExtraContainersTruth=["TruthEvents", 
                      "TruthParticles",
                      "TruthVertices",
                      "egammaTruthParticles",
                      "MuonTruthParticles"
                      ]

ExtraContainersPhotons=["Photons",
                        "GSFTrackParticles",
                        "egammaClusters",
                        "ForwardElectrons",
                        "ForwardElectronClusters"]

# for trigger studies
ExtraContainersTrigger=[
        "HLT_xAOD__ElectronContainer_egamma_Electrons",
        "HLT_xAOD__ElectronContainer_egamma_ElectronsAux.",
        "HLT_xAOD__PhotonContainer_egamma_Photons",
        "HLT_xAOD__PhotonContainer_egamma_PhotonsAux.",
        "HLT_xAOD__PhotonContainer_egamma_Iso_Photons",
        "HLT_xAOD__PhotonContainer_egamma_Iso_PhotonsAux.",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFex",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFexAux.",
        "HLT_xAOD__TrigPhotonContainer_L2PhotonFex",
        "HLT_xAOD__TrigPhotonContainer_L2PhotonFexAux.",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFexAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrig",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrigAux."
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_EFID",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_EFIDAux.",
        "LVL1EmTauRoIs",
        "LVL1EmTauRoIsAux.",
        "HLT_TrigPassBitsCollection_passbits",
        "HLT_TrigPassBitsCollection_passbitsAux.",
        "HLT_TrigPassFlagsCollection_passflags",
        "HLT_TrigPassFlagsCollection_passflagsAux.",
        "HLT_TrigRoiDescriptorCollection_initialRoI",
        "HLT_TrigRoiDescriptorCollection_initialRoIAux."
        ]

ExtraContainersTriggerDataOnly=[
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgamma",
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgammaAux.",
        "HLT_xAOD__CaloClusterContainer_TrigCaloClusterMaker",
        "HLT_xAOD__CaloClusterContainer_TrigCaloClusterMakerAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_FTF",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_FTFAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_L2ID",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_L2IDAux.",
        ]
