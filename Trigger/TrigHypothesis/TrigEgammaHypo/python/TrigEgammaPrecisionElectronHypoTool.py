#
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#
from AthenaCommon.SystemOfUnits import GeV
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.AllConfigFlags import ConfigFlags

from ROOT import egammaPID

def same( val , tool):
  return [val]*( len( tool.EtaBins ) - 1 )

#
# Create the hypo alg with all selectors
#
def createTrigEgammaPrecisionElectronHypoAlg(name, sequenceOut):
    acc = ComponentAccumulator()
    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool, defineHistogram
    MonTool = GenericMonitoringTool("MonTool_"+name)
  
    acc_ElectronCBSelectorTools = TrigEgammaPrecisionElectronCBSelectorCfg()
    acc_ElectronLHSelectorTools = TrigEgammaPrecisionElectronLHSelectorCfg()
    acc_ElectronDNNSelectorTools = TrigEgammaPrecisionElectronDNNSelectorCfg()

    acc.merge(acc_ElectronCBSelectorTools)
    acc.merge(acc_ElectronLHSelectorTools)
    acc.merge(acc_ElectronDNNSelectorTools)
  
    thePrecisionElectronHypo = CompFactory.TrigEgammaPrecisionElectronHypoAlg(name)
    thePrecisionElectronHypo.Electrons = str(sequenceOut)
    thePrecisionElectronHypo.RunInView = True
    thePrecisionElectronHypo.ElectronCBSelectorTools = acc_ElectronCBSelectorTools.getPublicTools()
    thePrecisionElectronHypo.ElectronLHSelectorTools = acc_ElectronLHSelectorTools.getPublicTools()
    thePrecisionElectronHypo.ElectronDNNSelectorTools = acc_ElectronDNNSelectorTools.getPublicTools()
    thePrecisionElectronHypo.CBNames = ["medium", "loose", "mergedtight"] # just like the pidnames
    thePrecisionElectronHypo.LHNames = ["lhtight", "lhmedium", "lhloose", "lhvloose", 
                                        "lhtight_nopix", "lhmedium_nopix","lhloose_nopix","lhvloose_nopix",
                                        "lhtight_nogsf", "lhmedium_nogsf","lhloose_nogsf","lhvloose_nogsf",
                                        "lhtight_nogsf_nopix", "lhmedium_nogsf_nopix","lhloose_nogsf_nopix","lhvloose_nogsf_nopix"] # just like the pidnames
    thePrecisionElectronHypo.DNNNames = ["dnntight", "dnnmedium", "dnnloose"] # just like the pidnames
    MonTool.Histograms = [ 
                defineHistogram('TIME_exec', type='TH1F', path='EXPERT', title="Precision Electron Hypo Algtime; time [ us ] ; Nruns", xbins=80, xmin=0.0, xmax=8000.0),
                defineHistogram('TIME_LH_exec', type='TH1F', path='EXPERT', title="Precision Electron Hypo LH Algtime; time [ us ] ; Nruns", xbins=20, xmin=0.0, xmax=2000),
                defineHistogram('TIME_DNN_exec', type='TH1F', path='EXPERT', title="Precision Electron Hypo DNN Algtime; time [ us ] ; Nruns", xbins=20, xmin=0.0, xmax=2000),
    ]
    MonTool.HistPath = 'PrecisionElectronHypo/'+name
    thePrecisionElectronHypo.MonTool=MonTool
    #acc.addEventAlgo(thePrecisionElectronHypo)
    return thePrecisionElectronHypo, acc

def TrigEgammaPrecisionElectronHypoAlgCfg(flags, name, inputElectronCollection ):
    from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
    acc = ComponentAccumulator()
    hypo_tuple = createTrigEgammaPrecisionElectronHypoAlg( name, inputElectronCollection )
    hypo_alg = hypo_tuple[0]
    hypo_acc = hypo_tuple[1]
    acc.addEventAlgo( hypo_alg )
    acc.merge(hypo_acc)
    acc.addService( CompFactory.AthONNX.ONNXRuntimeSvc())
    return acc

class TrigEgammaPrecisionElectronHypoToolConfig:


  __operation_points  = [  'tight'    ,
                           'medium'   ,
                           'loose'    ,
                           'vloose'   ,
                           'lhtight'  ,
                           'lhmedium' ,
                           'lhloose'  ,
                           'lhvloose' ,
                           'lhtight_nopix'  ,
                           'lhmedium_nopix' ,
                           'lhloose_nopix'  ,
                           'lhvloose_nopix' ,
                           'lhtight_nogsf'  ,
                           'lhmedium_nogsf' ,
                           'lhloose_nogsf'  ,
                           'lhvloose_nogsf' ,
                           'lhtight_nogsf_nopix'  ,
                           'lhmedium_nogsf_nopix' ,
                           'lhloose_nogsf_nopix'  ,
                           'lhvloose_nogsf_nopix' ,
                           'dnntight' ,
                           'dnnmedium',
                           'dnnloose' ,
                           'mergedtight',
                           'nopid',
                           ]

  __operation_points_lhInfo = [
        'nopix'
        ]

  __operation_points_gsfInfo = [
        'nogsf'
        ]

  # isolation cuts:w
  __isolationCut = {
        None: None,
        'ivarloose': 0.1,
        'ivarmedium': 0.065,
        'ivartight': 0.05
        }

  # LRT d0 cuts
  __lrtD0Cut = {
      '': -1.,
      None: None,
      'lrtloose':2.0,
      'lrtmedium':3.0,
      'lrttight':5.
      }


  def __init__(self, name, monGroups, cpart, tool=None):

    from AthenaCommon.Logging import logging
    self.__log = logging.getLogger('TrigEgammaPrecisionElectronHypoTool')
    self.__name = name
    self.__threshold = float(cpart['threshold'])
    self.__sel = cpart['addInfo'][0] if cpart['addInfo'] else cpart['IDinfo']
    self.__iso = cpart['isoInfo']
    self.__d0  = cpart['lrtInfo']
    self.__gsfInfo = cpart['gsfInfo']
    self.__lhInfo = cpart['lhInfo']
    self.__monGroups = monGroups
    
    if not tool:
      from AthenaConfiguration.ComponentFactory import CompFactory
      tool = CompFactory.TrigEgammaPrecisionElectronHypoTool( name )

    tool.EtaBins        = [0.0, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47]
    tool.ETthr          = same( self.__threshold*GeV, tool )
    tool.dETACLUSTERthr = 0.1
    tool.dPHICLUSTERthr = 0.1
    tool.RelPtConeCut   = -999
    tool.PidName        = ""
    tool.d0Cut          = -1
    tool.AcceptAll      = False
    tool.DoNoPid	= False
    self.__tool         = tool    

    self.__log.debug( 'Electron_Chain     :%s', self.__name )
    self.__log.debug( 'Electron_Threshold :%s', self.__threshold )
    self.__log.debug( 'Electron_Pidname   :%s', self.pidname() )
    self.__log.debug( 'Electron_iso       :%s', self.__iso )
    self.__log.debug( 'Electron_d0        :%s', self.__d0 )

  def chain(self):
    return self.__name

  #
  # Get the pidname
  #
  def pidname( self ):
    # if LLH, we should append the LH extra information if exist
    pidname = self.__sel

    extra = ""
    if 'lh' in self.__sel and self.__gsfInfo and self.__gsfInfo in self.__operation_points_gsfInfo:
      extra += '_' + self.__gsfInfo
    if 'lh' in self.__sel and self.__lhInfo and self.__lhInfo in self.__operation_points_lhInfo:
      extra += '_' + self.__lhInfo

    return pidname+extra

  def etthr(self):
    return self.__threshold

  def isoInfo(self):
    return self.__iso

  def d0Info(self):
    return self.__d0

  def gsfInfo(self):
    return self.__gsfInfo

  def tool(self):
    return self.__tool
  
  def nocut(self):

    self.__log.debug( 'Configure nocut' )
    self.tool().ETthr          = same( self.etthr()*GeV, self.tool())
    self.tool().dETACLUSTERthr = 9999.
    self.tool().dPHICLUSTERthr = 9999.

  def noPid(self):

    self.tool().DoNoPid = True
    self.__log.debug( 'Configure noPid' )
    self.tool().ETthr          = same( self.etthr()*GeV, self.tool())
    # No other cuts applied
    self.tool().dETACLUSTERthr = 9999.
    self.tool().dPHICLUSTERthr = 9999.

  #
  # LRT extra cut
  #
  def addLRTCut(self):
    if not self.d0Info() in self.__lrtD0Cut:
      self.__log.fatal(f"Bad LRT selection name: {self.d0Info()}")
    self.__tool.d0Cut = self.__lrtD0Cut[self.d0Info()]

  def acceptAll(self):
     self.tool().AcceptAll = True
  #
  # Isolation extra cut
  #
  def addIsoCut(self):
    if not self.isoInfo() in self.__isolationCut:
      self.__log.fatal(f"Bad Iso selection name: {self.isoInfo()}")
    self.tool().RelPtConeCut = self.__isolationCut[self.isoInfo()]

 

  def nominal(self):
    if not self.pidname() in self.__operation_points:
      self.__log.fatal("Bad selection name: %s" % self.pidname())
    self.tool().PidName = self.pidname()

  
  #
  # Compile the chain
  #
  def compile(self):

    if 'nocut' == self.pidname():
      self.nocut()
    elif 'nopid' == self.pidname():
      self.noPid()
    else: # nominal chain using pid selection
      self.nominal()


    # secundary cut configurations
    if self.isoInfo() and self.isoInfo()!="":
      self.addIsoCut()
    if self.d0Info() and self.d0Info()!="":
      self.addLRTCut()
    

    if hasattr(self.tool(), "MonTool"):
      
      doValidationMonitoring = ConfigFlags.Trigger.doValidationMonitoring # True to monitor all chains for validation purposes
      monGroups = self.__monGroups

      if (any('egammaMon:online' in group for group in monGroups) or doValidationMonitoring):
        self.addMonitoring()


  #
  # Create the monitoring code
  #
  def addMonitoring(self):

    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool, defineHistogram
    monTool = GenericMonitoringTool("MonTool_"+self.chain())
    monTool.Histograms = [ defineHistogram('dEta', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo #Delta#eta_{EF L1}; #Delta#eta_{EF L1}", xbins=80, xmin=-0.01, xmax=0.01),
                           defineHistogram('dPhi', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo #Delta#phi_{EF L1}; #Delta#phi_{EF L1}", xbins=80, xmin=-0.01, xmax=0.01),
                           defineHistogram('Et_em', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo cluster E_{T}^{EM};E_{T}^{EM} [MeV]", xbins=50, xmin=-2000, xmax=100000),
                           defineHistogram('Eta', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo entries per Eta;Eta", xbins=100, xmin=-2.5, xmax=2.5),
                           defineHistogram('Phi', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo entries per Phi;Phi", xbins=128, xmin=-3.2, xmax=3.2),
                           defineHistogram('EtaBin', type='TH1I', path='EXPERT', title="PrecisionElectron Hypo entries per Eta bin;Eta bin no.", xbins=11, xmin=-0.5, xmax=10.5),
                           defineHistogram('LikelihoodRatio', type='TH1F', path='EXPERT', title="PrecisionElectron Hypo LH", xbins=100, xmin=-5, xmax=5),
                           defineHistogram('mu', type='TH1F', path='EXPERT', title="Average interaction per crossing", xbins=100, xmin=0, xmax=100)]

    cuts=['Input','#Delta #eta EF-L1', '#Delta #phi EF-L1','eta','E_{T}^{EM}']

    monTool.Histograms += [ defineHistogram('CutCounter', type='TH1I', path='EXPERT', title="PrecisionElectron Hypo Passed Cuts;Cut",
                                            xbins=13, xmin=-1.5, xmax=12.5,  opt="kCumulative", xlabels=cuts) ]


    if ConfigFlags.Trigger.doValidationMonitoring:
      monTool.defineHistogram('ptcone20',type='TH1F',path='EXPERT',title= "PrecisionElectron Hypo ptcone20; ptcone20;", xbins=50, xmin=0, xmax=5.0),
      monTool.defineHistogram('relptcone20',type='TH1F',path='EXPERT',title= "PrecisionElectron Hypo; ptcone20/pt;", xbins=50, xmin=0, xmax=1),
      monTool.defineHistogram('ptvarcone20',type='TH1F',path='EXPERT',title= "PrecisionElectron Hypo ptvarcone20; ptvarcone20;", xbins=50, xmin=0, xmax=5.0),
      monTool.defineHistogram('relptvarcone20',type='TH1F',path='EXPERT',title= "PrecisionElectron Hypo; ptvarcone20/pt;", xbins=50, xmin=0, xmax=0.5)
      monTool.defineHistogram('trk_d0', type="TH1F", path='EXPERT', title="PrecisionElectron Hypo Track d0; d0 [mm]", xbins=100, xmin=-1, xmax=1)

    monTool.HistPath = 'PrecisionElectronHypo/'+self.chain()
    self.tool().MonTool = monTool


def _IncTool( name, monGroups, cpart, tool=None):
    config = TrigEgammaPrecisionElectronHypoToolConfig(name, monGroups, cpart, tool=tool)
    config.compile()
    return config.tool()



def TrigEgammaPrecisionElectronHypoToolFromDict( d , tool=None):
    """ Use menu decoded chain dictionary to configure the tool """
    cparts = [i for i in d['chainParts'] if ((i['signature']=='Electron') or (i['signature']=='Electron'))]
    name = d['chainName']
    monGroups = d['monGroups']
    return _IncTool( name, monGroups, cparts[0] , tool=tool )




#
# Electron DNN Selectors
#
def TrigEgammaPrecisionElectronDNNSelectorCfg(name='TrigEgammaPrecisionElectronDNNSelector', ConfigFilePath=None):
    acc = ComponentAccumulator()
    # We should include the DNN here
    if not ConfigFilePath:
      ConfigFilePath = ConfigFlags.Trigger.egamma.dnnVersion
  
    import collections.abc
    SelectorNames = collections.OrderedDict({
          'dnntight'  :'AsgElectronDNNTightSelector',
          'dnnmedium' :'AsgElectronDNNMediumSelector',
          'dnnloose'  :'AsgElectronDNNLooseSelector',
          })

    ElectronToolConfigFile = collections.OrderedDict({
          'dnntight'  :'ElectronDNNMulticlassTight.conf',
          'dnnmedium' :'ElectronDNNMulticlassMedium.conf',
          'dnnloose'  :'ElectronDNNMulticlassLoose.conf',
          })

    for dnnname, name in SelectorNames.items():
      SelectorTool = CompFactory.AsgElectronSelectorTool(name)
      SelectorTool.ConfigFile = ConfigFilePath + '/' + ElectronToolConfigFile[dnnname]
      SelectorTool.skipDeltaPoverP = True
      acc.addPublicTool(SelectorTool)

    return acc

#
# Electron LH Selectors
#
def TrigEgammaPrecisionElectronLHSelectorCfg( name='TrigEgammaPrecisionElectronLHSelector', ConfigFilePath=None, ConfigFileNoPixPath=None, ConfigFileNoGSFPath=None, ConfigFileNoGSFNoPixPath=None):

    # Configure the LH selectors
    acc = ComponentAccumulator()

    # Must be careful that order matches LHNames at the start of the file!
    import collections.abc
    SelectorConfigFiles = collections.OrderedDict({
        'lhtight'  : 'ElectronLikelihoodTightTriggerConfig',
        'lhmedium' : 'ElectronLikelihoodMediumTriggerConfig',
        'lhloose'  : 'ElectronLikelihoodLooseTriggerConfig',
        'lhvloose' : 'ElectronLikelihoodVeryLooseTriggerConfig'
    })

    VariationConfigInfos = collections.OrderedDict({
        '_default'     : {},
        '_nopix'       : {},
        '_nogsf'       : {},
        '_nogsf_nopix' : {}
    })

    VariationConfigInfos['_default']['postfix']      = ''
    VariationConfigInfos['_nopix']['postfix']        = '_NoPix'
    VariationConfigInfos['_nogsf']['postfix']        = ''
    VariationConfigInfos['_nogsf_nopix']['postfix']  = '_NoPix'
    
    VariationConfigInfos['_default']['path']      = ConfigFilePath           if ConfigFilePath           else ConfigFlags.Trigger.egamma.electronPidVersion
    VariationConfigInfos['_nopix']['path']        = ConfigFileNoPixPath      if ConfigFileNoPixPath      else ConfigFlags.Trigger.egamma.electronNoPixPidVersion
    VariationConfigInfos['_nogsf']['path']        = ConfigFileNoGSFPath      if ConfigFileNoGSFPath      else ConfigFlags.Trigger.egamma.electronNoGSFPidVersion
    VariationConfigInfos['_nogsf_nopix']['path']  = ConfigFileNoGSFNoPixPath if ConfigFileNoGSFNoPixPath else ConfigFlags.Trigger.egamma.electronNoGSFNoPixPidVersion
    
    from AthenaCommon.Logging import logging
    log = logging.getLogger('TrigEgammaPrecisionElectronHypoTool')
    log.debug( 'TrigEgammaPrecisionElectronLHSelectorCfg, order of LH tools:' )

    for pidvar, config in VariationConfigInfos.items():
        for pidname, configfilebase in SelectorConfigFiles.items():
            fullpidname = pidname if pidvar == '_default' else pidname + pidvar
            toolname = 'AsgElectronSelector_' + fullpidname
            configfile = config['path'] + '/' + configfilebase + config['postfix'] + '.conf'
            skipdeltapcheck = True if 'nogsf' in pidvar else False

            log.debug( ' --> %s, config file: %s', fullpidname, configfile )

            SelectorTool = CompFactory.AsgElectronLikelihoodTool(toolname)
            SelectorTool.ConfigFile = configfile
            SelectorTool.usePVContainer = False 
            SelectorTool.skipDeltaPoverP = skipdeltapcheck
            acc.addPublicTool(SelectorTool)

    return acc


#
# Electron CB Selectors
#

def TrigEgammaPrecisionElectronCBSelectorCfg(name='TrigEgammaPrecisionElectronCBSelector', ConfigFilePath=None):
    acc = ComponentAccumulator()
    from ElectronPhotonSelectorTools.TrigEGammaPIDdefs import BitDefElectron

    ElectronLooseHI = (0
            | 1 << BitDefElectron.ClusterEtaRange_Electron
            | 1 << BitDefElectron.ClusterHadronicLeakage_Electron
            | 1 << BitDefElectron.ClusterMiddleEnergy_Electron
            | 1 << BitDefElectron.ClusterMiddleEratio37_Electron
            | 1 << BitDefElectron.ClusterMiddleWidth_Electron
            | 1 << BitDefElectron.ClusterStripsWtot_Electron
    )

    ElectronMediumHI = (ElectronLooseHI
            | 1 << BitDefElectron.ClusterMiddleEratio33_Electron
            | 1 << BitDefElectron.ClusterBackEnergyFraction_Electron
            | 1 << BitDefElectron.ClusterStripsEratio_Electron
            | 1 << BitDefElectron.ClusterStripsDeltaEmax2_Electron
            | 1 << BitDefElectron.ClusterStripsDeltaE_Electron
            | 1 << BitDefElectron.ClusterStripsFracm_Electron
            | 1 << BitDefElectron.ClusterStripsWeta1c_Electron
    )

    if not ConfigFilePath:
        ConfigFilePath = ConfigFlags.Trigger.egamma.electronHIPidVersion

    from collections import OrderedDict
    SelectorNames = OrderedDict({
          'medium': 'AsgElectronIsEMSelectorHIMedium',
          'loose': 'AsgElectronIsEMSelectorHILoose',
          'mergedtight'  : 'AsgElectronIsEMSelectorMergedTight',
    })

    ElectronToolConfigFile = {
          'medium': 'ElectronIsEMMediumSelectorCutDefs.conf',
          'loose': 'ElectronIsEMLooseSelectorCutDefs.conf',
          'mergedtight'  : 'ElectronIsEMMergedTightSelectorCutDefs.conf',
    }

    ElectronMaskBits = {
          'medium': ElectronMediumHI,
          'loose': ElectronLooseHI,
          'mergedtight'  : egammaPID.ElectronTightHLT,
    }

    for sel, name in SelectorNames.items():
        SelectorTool = CompFactory.AsgElectronIsEMSelector(name)
        SelectorTool.ConfigFile = ConfigFilePath + '/' + ElectronToolConfigFile[sel]
        SelectorTool.isEMMask = ElectronMaskBits[sel]
        acc.addPublicTool(SelectorTool)
    
    return acc
