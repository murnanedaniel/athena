# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

""" Tau slice specific flags  """


 
from AthenaCommon.Logging import logging
from AthenaCommon.JobProperties import JobProperty, JobPropertyContainer, jobproperties
from TriggerMenu.menu.CommonSliceHelper import CommonSliceHelper, AllowedList

__author__  = ''
__version__="$Revision: 1.34 $"
__doc__="Tau slice specific flags  "

_flags = [] 
class doTrackingApproach(JobProperty):
    """ Use tracking approach in Tau Trigger """
    statusOn=True
    allowedTypes=['bool']
    StoredValue=False
    
_flags.append(doTrackingApproach)

class signatures(JobProperty):
    """ signatures in Tau slice """
    statusOn=True
    allowedTypes = ['list', 'str']
    StoredValue = []

_flags.append(signatures)



# create container

class TauSlice(JobPropertyContainer, CommonSliceHelper):
    """ Tau Slice Flags """

from TriggerJobOpts.TriggerFlags import TriggerFlags
TriggerFlags.add_Container(TauSlice)

# add add common slice flags
TriggerFlags.TauSlice.import_JobProperties('TriggerMenu.menu.CommonSliceFlags')

for flag in _flags:
    TriggerFlags.TauSlice.add_JobProperty(flag)
del _flags

# make an alias
TauSliceFlags = TriggerFlags.TauSlice
