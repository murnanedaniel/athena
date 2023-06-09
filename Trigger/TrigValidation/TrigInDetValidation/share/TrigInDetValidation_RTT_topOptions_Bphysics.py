#####################################################################################################
#
# top level jobOptions to run Muon chains in the RTT or standalone
# sets some global variables that adjust the execution of TrigInDetValidation_RTT_Common.py
#
# Jiri.Masik@manchester.ac.uk
#
#####################################################################################################

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

#set athenaCommonFlags.FilesInput to be able to use this job options standalone without RTT
#secondSet of files can be activated by the if statement below 

if athenaCommonFlags.FilesInput()==[]:
  athenaCommonFlags.FilesInput=[
    ]

  secondSet=True
  if secondSet: 
    athenaCommonFlags.FilesInput=[
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000005.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000006.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000018.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000019.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000030.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000033.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000047.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000048.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000049.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000064.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000069.pool.root.1",       
      "root://eosatlas//eos/atlas/atlascerngroupdisk/proj-sit/vchavda/TrigInDetValidation_muon/mc10_14TeV.106047.PythiaZmumu_no_filter.digit.RDO.e662_s1107_d459_tid285222_00/RDO.285222._000073.pool.root.1",       
    ]
  

rMC = False
if 'runMergedChain' in dir() and runMergedChain==True:
  rMC = True
rID=False
if 'doIDNewTracking' in dir() and doIDNewTracking==True:
  rID = True

(idtrigChainlist, tidaAnalysischains) = muonChains(rMC,rID)


def resetSigs():
  TriggerFlags.Slices_all_setOff()
  TriggerFlags.MuonSlice.setAll();
  TriggerFlags.MuonSlice.signatures = idTrigChainlist


PdgId=13


include("TrigInDetValidation/TrigInDetValidation_RTT_Common.py")


HLTMonManager = topSequence.HLTMonManager

from TrigIDJpsiMonitoring.TrigIDJpsiMonitoringConf import TrigJpsiMonTool
TPMonTool = TrigJpsiMonTool( name="TPMonTool")
ToolSvc += TPMonTool
HLTMonManager.AthenaMonTools += [ "TrigJpsiMonTool/TPMonTool" ]


TPMonTool.matchR = 0.1
TPMonTool.pTCutOffline = 4000.0
TPMonTool.d0Cut = 1000.0
TPMonTool.z0Cut = 2000.0
TPMonTool.d0CutOffline = 1000.0
TPMonTool.z0CutOffline = 2000.0
#TPMonTool.pTCutOffline = 4000.0
#TPMonTool.d0Cut = 1000.0
#TPMonTool.d0CutOffline = 1000.0
#TPMonTool.matchR = 1.0
#TPMonTool.z0CutOffline = 2000.0

#TPMonTool.Chain0 = "EF_2mu13_Zmumu_IDTrkNoCut"
TPMonTool.Chain0 = "EF_mu24i_tight";
#TPMonTool.Chain1 = "EF_2mu4T_Jpsimumu_IDTrkNoCut"
TPMonTool.Chain1 = ""
#TPMonTool.Chain2 = ""
#TPMonTool.Chain3 = ""

#TPMonTool.IDSCAN = "";
#TPMonTool.SiTrack = "";
TPMonTool.EF = "InDetTrigParticleCreation_Muon_EFID";
#TPMonTool.StratA = "TrigL2SiTrackFinder_MuonA";
#TPMonTool.StratB = "TrigL2SiTrackFinder_MuonB";
#TPMonTool.StratC = "TrigL2SiTrackFinder_MuonC";
TPMonTool.StratA = "";
TPMonTool.StratB = "";
TPMonTool.StratC = "";


#TPMonTool.TRTSF = "";
#TPMonTool.EFTRT = "";





