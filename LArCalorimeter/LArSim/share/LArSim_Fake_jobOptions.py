###############################################################
#
# Job options file
#
# Creates LArHits corresponding to fake clusters  in stores them
# in a separate container in StoreGate
#
#==============================================================
#--------------------------------------------------------------
# Private Application Configuration options
#--------------------------------------------------------------
theApp.Dlls += ["LArSim"] 
theApp.TopAlg += ["LArHitMaker/hitmaker1"] 
#--------------------------------------------------------------
# Algorithms Private Options
#--------------------------------------------------------------
# to communicate the Producer class name and instance name to hitmaker1
hitmaker1 = Algorithm( "hitmaker1" )
hitmaker1.HitRetrieverNameAndType = "LArFakeClusterProducer/producer1" 
# to communicate the location in the TDS where the LArHitContainer will be stored 
hitmaker1.HitContainerLocation = "LArFakeEMBHit" 
# to switch on the Debug mode , in which the LArHitContainter content is fully printed out
hitmaker1.DebugFlag = 1 
# total energy in fake cluster
producer1 = Algorithm( "producer1" )
producer1.Etot = 100. 
# spacing between 2 clusters in eta ( in middle cell units )
producer1.ClusterSpacing = 4 
# middle cell of first cluster
producer1.EtaStart = 2  
# Sampling fractions 
producer1.SampFracEMBPS = 0.2 
producer1.SampFracEMB1  = 0.2 
producer1.SampFracEMB2  = 0.2 
producer1.SampFracEMB3  = 0.2 
#==============================================================
#
# End of job options file
#
###############################################################
