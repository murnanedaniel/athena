# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration

from AthenaCommon.SystemOfUnits import GeV
from AthenaConfiguration.AllConfigFlags import ConfigFlags
from ROOT import egammaPID

#
# photon hypo alg
#
def createTrigEgammaPrecisionPhotonHypoAlg(name, sequenceOut):

  from TrigEgammaHypo.TrigEgammaHypoConf import TrigEgammaPrecisionPhotonHypoAlg
  thePrecisionPhotonHypo = TrigEgammaPrecisionPhotonHypoAlg(name)
  thePrecisionPhotonHypo.IsEMNames = ['tight','medium','loose']
  thePrecisionPhotonHypo.PhotonIsEMSelectorTools = createTrigEgammaPrecisionPhotonSelectors()
  thePrecisionPhotonHypo.Photons = sequenceOut
  return thePrecisionPhotonHypo


def same( val , tool):
  return [val]*( len( tool.EtaBins ) - 1 )
  

#
# For photons
#
class TrigEgammaPrecisionPhotonHypoToolConfig:

  __operation_points  = [  'tight'    , 
                           'medium'   , 
                           'loose'    , 
                           ]

  def __init__(self, name, monGroups, cpart, tool=None):

    from AthenaCommon.Logging import logging
    self.__log = logging.getLogger('TrigEgammaPrecisionPhotonHypoTool')
    self.__name       = name
    self.__threshold  = float(cpart['threshold']) 
    self.__sel        = cpart['addInfo'][0] if cpart['addInfo'] else cpart['IDinfo']
    self.__monGroups = monGroups
    
    if not tool:
      from TrigEgammaHypo.TrigEgammaHypoConf import TrigEgammaPrecisionPhotonHypoTool    
      tool = TrigEgammaPrecisionPhotonHypoTool( name )
     
    tool.EtaBins        = [0.0, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47]
    tool.ETthr          = same( self.__threshold*GeV , tool)
    tool.dETACLUSTERthr = 0.1
    tool.dPHICLUSTERthr = 0.1
    tool.PidName        = ""

    self.__tool = tool
    self.__log.debug( 'Chain     :%s', self.__name )
    self.__log.debug( 'Threshold :%s', self.__threshold )
    self.__log.debug( 'Pidname   :%s', self.__sel )


  def chain(self):
    return self.__name
  
  def pidname( self ):
    return self.__sel

  def etthr(self):
    return self.__threshold

  def tool(self):
    return self.__tool
  

  def etcut(self):
    self.__log.debug( 'Configure etcut' )
    self.tool().ETthr          = same( ( self.etthr() -  3 )*GeV, self.tool() )
    # No other cuts applied
    self.tool().dETACLUSTERthr = 9999.
    self.tool().dPHICLUSTERthr = 9999.


  def nominal(self):
    if not self.pidname() in self.__operation_points:
      self.__log.fatal("Bad selection name: %s" % self.pidname())
    self.tool().PidName = self.pidname()


  #
  # Compile the chain
  #
  def compile(self):

    if 'etcut' == self.pidname():
      self.etcut()

    else:
      self.nominal()

    if hasattr(self.tool(), "MonTool"):
      
      doValidationMonitoring = ConfigFlags.Trigger.doValidationMonitoring # True to monitor all chains for validation purposes
      monGroups = self.__monGroups

      if (any('egammaMon:online' in group for group in monGroups) or doValidationMonitoring):
        self.addMonitoring()


  #
  # Monitoring code
  #
  def addMonitoring(self):
    
    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
    monTool = GenericMonitoringTool("MonTool_"+self.__name,
                                    HistPath = 'PrecisionPhotonHypo/'+self.__name)
    monTool.defineHistogram('dEta', type='TH1F', path='EXPERT', title="PrecisionPhoton Hypo #Delta#eta_{EF L1}; #Delta#eta_{EF L1}", xbins=80, xmin=-0.01, xmax=0.01)
    monTool.defineHistogram('dPhi', type='TH1F', path='EXPERT', title="PrecisionPhoton Hypo #Delta#phi_{EF L1}; #Delta#phi_{EF L1}", xbins=80, xmin=-0.01, xmax=0.01)
    monTool.defineHistogram('Et_em', type='TH1F', path='EXPERT', title="PrecisionPhoton Hypo cluster E_{T}^{EM};E_{T}^{EM} [MeV]", xbins=50, xmin=-2000, xmax=100000)
    monTool.defineHistogram('Eta', type='TH1F', path='EXPERT', title="PrecisionPhoton Hypo entries per Eta;Eta", xbins=100, xmin=-2.5, xmax=2.5)
    monTool.defineHistogram('Phi', type='TH1F', path='EXPERT', title="PrecisionPhoton Hypo entries per Phi;Phi", xbins=128, xmin=-3.2, xmax=3.2)
    monTool.defineHistogram('EtaBin', type='TH1I', path='EXPERT', title="PrecisionPhoton Hypo entries per Eta bin;Eta bin no.", xbins=11, xmin=-0.5, xmax=10.5)

    cuts=['Input','#Delta #eta EF-L1', '#Delta #phi EF-L1','eta','E_{T}^{EM}']

    monTool.defineHistogram('CutCounter', type='TH1I', path='EXPERT', title="PrecisionPhoton Hypo Passed Cuts;Cut",
                            xbins=13, xmin=-1.5, xmax=12.5,  opt="kCumulative", xlabels=cuts)

    if ConfigFlags.Trigger.doValidationMonitoring:
      monTool.defineHistogram('etcone20',type='TH1F',path='EXPERT',title= "PrecisionPhoton Hypo etcone20; etcone20;", xbins=50, xmin=0, xmax=5.0)
      monTool.defineHistogram('topoetcone20',type='TH1F',path='EXPERT',title= "PrecisionPhoton Hypo; topoetcone20;", xbins=50, xmin=-10, xmax=10)
      monTool.defineHistogram('reletcone20',type='TH1F',path='EXPERT',title= "PrecisionPhoton Hypo etcone20/et; etcone20/et;", xbins=50, xmin=-0.5, xmax=0.5)
      monTool.defineHistogram('reltopoetcone20',type='TH1F',path='EXPERT',title= "PrecisionPhoton Hypo; topoetcone20/pt;", xbins=50, xmin=-0.5, xmax=0.5)

    self.tool().MonTool = monTool



def _IncTool( name,monGroups, cpart, tool=None ):
    config = TrigEgammaPrecisionPhotonHypoToolConfig(name, monGroups, cpart, tool=tool)
    config.compile()
    return config.tool()

 

def TrigEgammaPrecisionPhotonHypoToolFromDict( d , tool=None):
    """ Use menu decoded chain dictionary to configure the tool """
    cparts = [i for i in d['chainParts'] if ((i['signature']=='Electron') or (i['signature']=='Photon'))] 
    name = d['chainName']
    monGroups = d['monGroups'] 
    return _IncTool( name, monGroups, cparts[0], tool=tool )
                   
    


#
# Photon IsEM selectors
#
def createTrigEgammaPrecisionPhotonSelectors(ConfigFilePath=None):

    from ElectronPhotonSelectorTools.ConfiguredAsgPhotonIsEMSelectors import ConfiguredAsgPhotonIsEMSelector

    if not ConfigFilePath:
      ConfigFilePath = ConfigFlags.Trigger.egamma.photonPidVersion

    import collections.abc
    # Configure the IsEM selectors
    SelectorNames = collections.OrderedDict( {
            'tight'  : 'TightPhotonSelector',
            'medium' : 'MediumPhotonSelector',
            'loose'  : 'LoosePhotonSelector',
            } )
    SelectorPID = {
            'loose'  : egammaPID.PhotonIDLoose,
            'medium' : egammaPID.PhotonIDMedium,
            'tight'  : egammaPID.PhotonIDTight,
            }
    PhotonToolConfigFile = {
            'loose'  : 'PhotonIsEMLooseSelectorCutDefs.conf', 
            'medium' : 'PhotonIsEMMediumSelectorCutDefs.conf', 
            'tight'  : 'PhotonIsEMTightSelectorCutDefs.conf',
            } 
    PhotonIsEMBits = {
            'loose'  : egammaPID.PhotonLoose,
            'medium' : egammaPID.PhotonMedium,
            'tight'  : egammaPID.PhotonTight,
            }
    selectors = []
    for sel, name in SelectorNames.items():
        SelectorTool = ConfiguredAsgPhotonIsEMSelector(name, SelectorPID[sel])
        ConfigFile = ConfigFilePath + '/' + PhotonToolConfigFile[sel] 
        SelectorTool.ConfigFile = ConfigFile
        SelectorTool.ForceConvertedPhotonPID = True
        SelectorTool.isEMMask = PhotonIsEMBits[sel] 
        selectors.append(SelectorTool)

    return selectors
