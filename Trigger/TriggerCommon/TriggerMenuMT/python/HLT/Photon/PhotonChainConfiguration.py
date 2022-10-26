#
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

from AthenaCommon.Logging import logging
logging.getLogger().info("Importing %s",__name__)
log = logging.getLogger(__name__)

from ..Config.ChainConfigurationBase import ChainConfigurationBase
from AthenaConfiguration.ComponentFactory import isRun3Cfg

if isRun3Cfg():
    pass
else:
    from ..CommonSequences.CaloSequences import fastCaloMenuSequence
    from ..Photon.FastPhotonMenuSequences import fastPhotonMenuSequence
    from ..Photon.PrecisionPhotonMenuSequences import precisionPhotonMenuSequence
    from ..Photon.PrecisionPhotonCaloIsoMenuSequences import precisionPhotonCaloIsoMenuSequence
    from ..Photon.PrecisionCaloMenuSequences import precisionCaloMenuSequence
    from ..Photon.HipTRTMenuSequences import hipTRTMenuSequence
    from TrigEgammaHypo.TrigEgammaHypoConf import TrigEgammaTopoHypoTool



from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool, defineHistogram
#----------------------------------------------------------------
# fragments generating configuration will be functions in New JO,
# so let's make them functions already now
#----------------------------------------------------------------
def fastPhotonCaloSequenceCfg( flags, doRinger = False ):
    return fastCaloMenuSequence(flags, 'Photon', doRinger=doRinger)
    
def fastPhotonSequenceCfg( flags ):    
    return fastPhotonMenuSequence( flags )

def _diPhotonComboHypoToolFromDict(chainDict, lowermass=80000,uppermass=-999,dphi=1.5,applymass=False,applydphi=False): 
    name = chainDict['chainName']
    monTool = GenericMonitoringTool("MonTool_"+name)
    monTool.Histograms = [
        defineHistogram('DphiOfAccepted', type='TH1F', path='EXPERT', title="PrecisionCalo Hypo entries per Phi;Phi", xbins=128, xmin=-3.2, xmax=3.2),
        defineHistogram('MassOfAccepted', type='TH1F', path='EXPERT', title="Mass in accepted combinations [MeV]", xbins=75, xmin=0, xmax=150000)
    ]
    tool= TrigEgammaTopoHypoTool(name)
    tool.AcceptAll = False
    tool.ApplyMassCut = applymass
    tool.LowerMassEgammaClusterCut = lowermass
    tool.UpperMassEgammaClusterCut = uppermass
    tool.ApplyDPhiCut = applydphi
    tool.ThresholdDPhiCut = dphi
    monTool.HistPath = 'EgammaMassHypo/'+tool.getName()
    tool.MonTool = monTool
    return tool

def diphotonDPhiHypoToolFromDict(chainDict):
    return _diPhotonComboHypoToolFromDict(chainDict,lowermass=80000,uppermass=-999,dphi=1.5,applymass=False,applydphi=True)

def diphotonDPhiMassHypoToolFromDict(chainDict):
    return _diPhotonComboHypoToolFromDict(chainDict,lowermass=80000,uppermass=-999,dphi=1.5,applymass=True,applydphi=True)


#----------------------------------------------------------------
# Class to configure chain
#----------------------------------------------------------------
class PhotonChainConfiguration(ChainConfigurationBase):

    def __init__(self, chainDict):
        ChainConfigurationBase.__init__(self,chainDict)
        self.chainDict = chainDict

    # ----------------------
    # Prepare the sequence
    # ----------------------
    def prepareSequence(self):
        # This function prepares the list of step names from which assembleChainImpl would make the chain assembly from.
        
        # --------------------
        # define here the names of the steps and obtain the chainStep configuration
        # --------------------

        stepNames = [] # This will contain the name of the steps we will want to configure
        # Put first fast Calo. Two possible variants now: 
        # Step 1
        stepNames += ['getFastCalo']

        # OK now, unless its a HipTRT chain we need to do fastPhoton here:
        if self.chainPart['extra'] == 'hiptrt':
            stepNames += ['getHipTRT']
            # for hiptrt chains, there is noprecision Calo nor precision Photon so returning sequence here:
            return stepNames

        # Now we do fast Photon
        # Step 2
        stepNames += ['getFastPhoton']


        # After we need to place precisionCalo. There is no chain (except hiptrt) that does not require precision calo. Otherwise please insert logic here:
        # Step 3
        stepNames += ['getPrecisionCaloPhoton']

        # And we will do precisionPhoton UNLESS its an etcut chain
        if 'etcut' in self.chainPart['IDinfo']:
            # if its an etcut chain we return the sequence up to here
            return stepNames


        #Now its the turn of precision Photon. Here we apply cutbased identification
        #Step 4
        stepNames += ['getPrecisionPhoton']

        # Finally we need to run isolaton *IF* its an isolated chain
        # Step 5
        if 'noiso' in self.chainPart['isoInfo'] or 'icaloloose' in self.chainPart['isoInfo'] or 'icalomedium' in self.chainPart['isoInfo'] or 'icalotight' in self.chainPart['isoInfo']:
            stepNames += ['getPhotonCaloIso']

        return stepNames



    # ----------------------
    # Assemble the chain depending on information from chainName
    # ----------------------
    def assembleChainImpl(self):
        log.debug("Assembling chain for %s", self.chainName)
        # This will contain the name of the steps we will want to configure
        steps = self.prepareSequence()

        chainSteps = []

        log.debug("stepNames: %s", steps)
        for step in steps:
            log.debug('Adding photon trigger step %s', step)
            chainstep = getattr(self, step)()
            chainSteps+=[chainstep]

        myChain = self.buildChain(chainSteps)
        return myChain


    # --------------------
    # Configuration of steps
    # --------------------
    def getFastCalo(self):
        stepName = "PhotonFastCalo"
        doRinger = 'ringer' in self.chainPart['L2IDAlg']

        return self.getStep(1,stepName,[ fastCaloMenuSequence], name = 'Photon', doRinger = doRinger)

    def getFastPhoton(self):
        stepName = "FastPhoton"
        return self.getStep(2,stepName,[ fastPhotonMenuSequence])

    def getPrecisionCaloPhoton(self):
        do_ion = 'ion' in self.chainPart['extra']
        if do_ion:
            stepName = "PhotonPrecisionHICalo"
        else:
            stepName = "PhotonPrecisionCalo"

        return self.getStep(3,stepName,[ precisionCaloMenuSequence], name = 'Photon', ion=do_ion)
    
    def getHipTRT(self):
        stepName = "hipTRT"
        return self.getStep(2,stepName,[ hipTRTMenuSequence])

    def getPrecisionPhoton(self):

        stepName = "precision_photon"
        do_ion = 'ion' in self.chainPart['extra'] == 'ion'

        if do_ion:
            stepName += '_ion'
        

        return self.getStep(4,stepName,sequenceCfgArray=[precisionPhotonMenuSequence], name = 'Photon',  ion=do_ion)

    def getPhotonCaloIso(self):

        stepName = "precision_photon_CaloIso"
        comboTools = []
        do_ion = 'ion' in self.chainPart['extra']

        if do_ion:
            stepName += '_ion'
        
        if "dPhi15" in self.chainDict["topo"]:
            stepName+='_dPhi15'
            if "m80" in self.chainDict["topo"]:
                stepName+= '_m80'
                comboTools.append(diphotonDPhiMassHypoToolFromDict)
            else:
                comboTools.append(diphotonDPhiHypoToolFromDict)

       

        return self.getStep(5,stepName,sequenceCfgArray=[precisionPhotonCaloIsoMenuSequence], name = 'Photon', comboTools=comboTools, ion=do_ion)
    
