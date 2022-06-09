#                                                                
#Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
# 

#import sTGCRawMonLabels as labels
import MMRawDataMonitoring.MMRawMonLabels as labels

def getsTGCLabel(x,y):

    labelx = getattr(labels, x)
    labely = getattr(labels, y)
    return labelx,labely

def getsTGCLabelY(y):

    labely = getattr(labels, y)
    return labely

def getsTGCLabelX(x):

    labelx = getattr(labels, x)
    return labelx
