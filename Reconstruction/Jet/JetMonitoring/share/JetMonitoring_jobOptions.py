# Inspired by PhysVal_jobOptions.py

jetMonMan = AthenaMonManager( "JetMonManager",
                           FileKey = DQMonFlags.monManFileKey(),
                           Environment = DQMonFlags.monManEnvironment(),
                           ManualDataTypeSetup = DQMonFlags.monManManualDataTypeSetup(),
                           DataType = DQMonFlags.monManDataType() )
topSequence += jetMonMan
 
 

from JetMonitoring.JetMonitoringHistos import athenaMonitoringTools


jetMonMan.AthenaMonTools += athenaMonitoringTools() 
