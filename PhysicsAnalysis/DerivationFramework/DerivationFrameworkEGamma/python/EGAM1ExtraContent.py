# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration


#Content included in addition to the Smart Slimming Content

ExtraContentElectrons=[]

ExtraContentMuons=[
    "Muons.ptcone20",
    "Muons.ptcone30",
    "Muons.ptcone40",
    "Muons.etcone20",
    "Muons.etcone30",
    "Muons.etcone40"
]   

ExtraMuonsTruth=[
    "MuonTruthParticles.e",
    "MuonTruthParticles.px",
    "MuonTruthParticles.py",
    "MuonTruthParticles.pz",
    "MuonTruthParticles.status",
    "MuonTruthParticles.pdgId",
    "MuonTruthParticles.truthOrigin",
    "MuonTruthParticles.truthType"
    ]

ExtraContentPhotons=[
    ]

ExtraContentPrimaryVertices=["PrimaryVertices.sumPt2"]

ExtraPhotonsTruth=[
    "Photons.truthOrigin",
    "Photons.truthType",
    "Photons.truthParticleLink"
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

#ExtraContentHLTElectrons=[
#        "HLT_xAOD__ElectronContainer_egamma_Electrons.e.pt.Rhad.Rhad1.e277.Reta.Rphi.weta2.f1.fracs1.wtots1.weta1.DeltaE.Eratio.caloClusterLinks"
#]

ExtraContentTrackJets=["AntiKt4PV0TrackJets.pt.eta.phi.e.m.btaggingLink.constituentLinks"]

# not in R22 yet for AntiKt4PV0TrackJets
ExtraContentBtagging=[]
#ExtraContentBtagging=["BTagging_AntiKt4Track.SV1_pb.SV1_pc.SV1_pu.IP2D_pb.IP2D_pc.IP2D_pu.IP3D_pb.IP3D_pc.IP3D_pu.JetFitter_pb.JetFitter_pc.JetFitter_pu.JetFitterCombNN_pb.JetFitterCombNN_pc.JetFitterCombNN_pu.MV2c00_discriminant.MV2c10_discriminant.MV2c20_discriminant"]


from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import GainDecorator, getGainDecorations
GainDecoratorTool = GainDecorator()
ExtraContentPhotons.extend( getGainDecorations(GainDecoratorTool) )
ExtraContentElectrons.extend( getGainDecorations(GainDecoratorTool) )

ExtraContentAll=ExtraContentElectrons+ExtraContentMuons+ExtraContentPhotons+ExtraContentGSFConversionVertices+ExtraContentPrimaryVertices+ExtraContentTrackJets+ExtraContentBtagging

ExtraContentAllTruth=ExtraMuonsTruth+ExtraPhotonsTruth

ExtraContainersTruth=["TruthEvents", 
                      "TruthParticles",
                      "TruthVertices",
                      "egammaTruthParticles"
                      ]

ExtraContainersElectrons=["Electrons",
                          "GSFTrackParticles",
                          "egammaClusters",
                          "CaloCalTopoClusters"]

ExtraContainersReweightedElectrons=[
                          "NewSwElectrons",    # only if DoCellReweighting is ON
                          "MaxVarSwElectrons", # if variations are ON
                          "MinVarSwElectrons"  # if variations are ON
                          ]

# for trigger studies
ExtraContainersTrigger=[
        "HLT_xAOD__ElectronContainer_egamma_Electrons",
        "HLT_xAOD__ElectronContainer_egamma_ElectronsAux.",
        "HLT_xAOD__PhotonContainer_egamma_Photons",
        "HLT_xAOD__PhotonContainer_egamma_PhotonsAux.",
        "HLT_xAOD__TrigRingerRingsContainer_TrigT2CaloEgamma",
        "HLT_xAOD__TrigRingerRingsContainer_TrigT2CaloEgammaAux.",
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgamma",
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgammaAux.",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFexAux.",
        "HLT_xAOD__TrigRNNOutputContainer_TrigRingerNeuralFex",
        "HLT_xAOD__TrigRNNOutputContainer_TrigRingerNeuralFexAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrig",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrigAux.",
        "HLT_xAOD__TrigPassBitsContainer_passbits",
        "HLT_xAOD__TrigPassBitsContainer_passbitsAux.",
        "LVL1EmTauRoIs",
        "LVL1EmTauRoIsAux.",
        "HLT_TrigRoiDescriptorCollection_initialRoI",
        "HLT_TrigRoiDescriptorCollection_initialRoIAux.",
        "HLT_xAOD__RoiDescriptorStore_initialRoI",
        "HLT_xAOD__RoiDescriptorStore_initialRoIAux.",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFex",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFexAux."
        ]

ExtraContainersTriggerDataOnly=[]

ExtraVariablesEventShape=[
    "TopoClusterIsoCentralEventShape.Density",
    "TopoClusterIsoForwardEventShape.Density",
]

#should slim electron/fwdelectrons/cluster collections and keep only relevant subset of variables..
