# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

# $Id: xAODTrigSpacePointCountsCreator.py 578517 2014-01-15 13:48:04Z krasznaa $

# Import the package's configurables:
from xAODTrigMinBiasCnv.xAODTrigMinBiasCnvConf import *

## Helper function for creating xAOD TrigSpacePointCounts objects
def xAODTrigSpacePointCountsCreator( sequence = None, stream = None,
                             key = "HLT_spacepoints" ):
    """
    Helper function that schedules all the components that create xAOD
    TrigSpacePointCounts objects out of existing AOD objects.

    Arguments:
      sequence: Optional sequence to which the components should be added.
                If not specified, the main algorithm sequence is used.
      stream: Optional stream that the objects should be written to.
              If not specified, the objects are added to StreamAOD.
      key: Optional StoreGate key for the AOD object that is supposed to be
           translated.
    """

    # Create a logger for the function:
    if "logger" in dir(): orig_logger = logger
    from AthenaCommon.Logging import logging
    logger = logging.getLogger( "xAODTrigSpacePointCountsCreator" )

    # Tell the user what's happening:
    logger.info( "Creating xAOD::TrigSpacePointCounts from TrigDec::TrigSpacePointCounts" )

    # Get the main sequence if necessary:
    if sequence == None:
        from AthenaCommon.AlgSequence import AlgSequence
        sequence = AlgSequence()
        pass

    # Access the stream if necessary:
    if stream == None:
        from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
        stream = MSMgr.GetStream( "StreamAOD" )
        pass

    # Schedule the converter algorithm:
    alg = xAODMaker__TrigSpacePointCountsCnvAlg()
    alg.AODKey = key
    alg.xAODKey = key
    sequence += alg

    # Add the created objects to the output:
    stream.AddItem( "xAOD::TrigSpacePointCountsContainer_v1#%s" % key )
    stream.AddItem( "xAOD::TrigSpacePointCountsAuxContainer_v1#%sAux." % key )

    # Reinstate the old logger if it existed:
    if "orig_logger" in dir(): logger = orig_logger

    pass
