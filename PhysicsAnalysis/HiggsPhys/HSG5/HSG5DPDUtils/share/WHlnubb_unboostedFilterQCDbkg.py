'''
Created on 28 Feb 2012

@author: wvazquez
'''

from D3PDMakerConfig.D3PDMakerFlags import D3PDMakerFlags
from AthenaCommon.AlgSequence import AlgSequence,AthSequencer
from ROOT import egammaPID

sequencer = AthSequencer("HSG5WHQ_Sequencer",
                         StopOverride=False)



from HSG5DPDUtils.HSG5DPDUtilsConf import HSG5__LeptonFilter
leptonFilter=HSG5__LeptonFilter("HSG5WHQ_SingleLeptonFilter",
                                ElectronFilterNameAndType="D2PDElectronSelector/HSG5WHQ_ElectronSelector",
                                MuonFilterNameAndTypeVec=["D2PDMuonSelector/HSG5WHQ_MuidMuonSelector",
                                                          "D2PDMuonSelector/HSG5WHQ_StacoMuonSelector",
                                                          "D2PDMuonSelector/HSG5WHQ_ThirdChainMuonSelector"])

# create electron and muon selectors
from D2PDMaker.D2PDMakerConf import D2PDElectronSelector
leptonFilter += D2PDElectronSelector( "HSG5WHQ_ElectronSelector",
                                      inputCollection      = 'ElectronAODCollection',
                                      outputLinkCollection = 'HSG5WHQ_LooseElectronLinkCollection',
                                      minNumberPassed      = 1,
                                      ptMin                = 20.0*Units.GeV,
                                      clusterEtaMin        = -2.5,
                                      clusterEtaMax        = 2.5,
                                      electronIsEM         = egammaPID.ElectronMedium )

from HSG5DPDUtils.HSG5Selectors import MuonSelector
muSelector = MuonSelector( minNumberPassed = 1,
                           ptMin = 18.0,
                           absEtaMax = 2.5,
                           acceptIsCombined = True,
                           acceptIsSegmentTagged = True,
                           doRelPtCone20  = True,
                           relPtCone20Max = 0.5 )

leptonFilter += muSelector.getMuonSelector('HSG5WHQ_MuidMuonSelector','MuidMuonCollection',
                                           'HSG5WHQ_LooseMuidMuonLinkCollection')
leptonFilter += muSelector.getMuonSelector('HSG5WHQ_StacoMuonSelector','StacoMuonCollection',
                                           'HSG5WHQ_LooseStacoMuonLinkCollection')
leptonFilter += muSelector.getMuonSelector('HSG5WHQ_ThirdChainMuonSelector','Muons',
                                           'HSG5WHQ_LooseThirdChainMuonLinkCollection')

sequencer += leptonFilter

HSG5D3PD_Stream.RequireAlgs.append("HSG5WHQ_SingleLeptonFilter")

# jet selector
from AssociationComps.AssociationCompsConf import DeltaRAssociationTool
ToolSvc += DeltaRAssociationTool( "HSG5WHQ_emDeltaRAssociationTool",
                                  OutputLevel = INFO,
                                  inputAssociateToCollection = 'HSG5WHQ_LooseElectronLinkCollection',
                                  deltaRMax = 0.3,
                                  writeUserData = False)

from D2PDMaker.D2PDMakerConf import D2PDJetSelector
sequencer += D2PDJetSelector( "HSG5WHQ_JetFilter",
                              inputCollection      = 'AntiKt4TopoEMJets',
                              outputLinkCollection = 'HSG5WHQ_JetLinkCollection',
                              associationToolList = [ ToolSvc.HSG5WHQ_emDeltaRAssociationTool ],
                              outputAssociationContainerList = [ "HSG5WHQ_jetsMatchedToElectrons" ],
                              numberOfAssociationsMaxCutList = [ 0 ],
                              minNumberPassed      = 1,
                              ptMin                = 20.0*Units.GeV,
                              etaMax               = 2.5)

HSG5D3PD_Stream.RequireAlgs.append("HSG5WHQ_JetFilter")

if False:
    # (for private production ony) insert in beginning of PreD3PDSequencer
    mainSequencer = AlgSequence(D3PDMakerFlags.PreD3PDAlgSeqName(),
                                StopOverride = False)

    if not hasattr( topSequence, D3PDMakerFlags.PreD3PDAlgSeqName() ):
        topSequence += mainSequencer

    mainSequencer.insert(0,sequencer)

else:
    topSequence+=sequencer
