# dictionary defining a watch list for MR notifications
# keys   ... must be valid regular expression
# values ... must be sets containing Gitlab usernames as strings
WATCH_LIST = {}

WATCH_LIST['^CI$']       = set(['cgumpert'])
WATCH_LIST['Simulation'] = set(['ritsch','jchapman','vpascuzz'])
WATCH_LIST['Digitization'] = set(['jchapman'])
WATCH_LIST['Overlay'] = set(['jchapman','ahaas','tkharlam'])
WATCH_LIST['TrkiPatFitter'] = set(['pop'])
WATCH_LIST['MooPerformance'] = set(['pop'])
WATCH_LIST['JetCalibTools'] = set(['jbossios'])
WATCH_LIST['AFP'] = set(['ggach'])
WATCH_LIST['BTagging'] = set(['cpollard', 'guirriec'])
WATCH_LIST['^Database/A'] = set(['mnowak'])
WATCH_LIST['^Database/TPTools'] = set(['mnowak'])
WATCH_LIST['^Database/PersistentDataModel'] = set(['mnowak'])
WATCH_LIST['^Control/'] = set(['ssnyder'])
WATCH_LIST['^Database/'] = set(['ssnyder'])
WATCH_LIST['MuonSpectrometer'] = set(['jomeyer','wleight'])
WATCH_LIST['MuonIdentification'] = set(['jomeyer','wleight'])
WATCH_LIST['AGDD'] = set(['jomeyer'])
WATCH_LIST['(eflow)|(PFlow)|(PFO)'] = set(['mhodgkin'])
WATCH_LIST['InDetBeamSpot'] = set(['csuster'])
WATCH_LIST['MuonEfficiencyCorrections'] = set(['nkoehler','jojungge'])
WATCH_LIST['MuonTPTools'] = set(['nkoehler','jojungge'])
WATCH_LIST['MuonPerformanceAlgs'] = set(['nkoehler','jojungge','gabarone'])
WATCH_LIST['MuonPerformanceHistUtils'] = set(['nkoehler','jojungge'])
WATCH_LIST['IsolationSelection'] = set(['maklein','jojungge','jpoveda','dzhang'])
WATCH_LIST['Trigger/TrigCost'] = set(['cmcnicol'])
WATCH_LIST['TrigInDetAnalysis.*'] = set(['sutt'])
WATCH_LIST['TrigInDetValidation.*'] = set(['sutt','hartj'])
WATCH_LIST['TrigIDtrkMonitoring'] = set(['sutt'])
WATCH_LIST['.*RoiDescriptor'] = set(['sutt'])
WATCH_LIST['.*RegionSelector'] = set(['sutt'])
WATCH_LIST['RegSelLUT'] = set(['sutt'])
WATCH_LIST['(PixelMonitoring)|(PixelConditionsServices)|(PixelRawData)'] = set(['kzoch','ibragimo'])
WATCH_LIST['Trigger/TrigFTK'] = set(['karolos', 'benhoob', 'mswiatlo', 'jahreda'])
WATCH_LIST['PhysicsAnalysis/SUSYPhys'] = set(['zmarshal','szambito'])
WATCH_LIST['MuonMomentumCorrections'] =set(['gabarone'])
