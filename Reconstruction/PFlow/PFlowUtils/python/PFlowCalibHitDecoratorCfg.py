# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def PFlowCalibHitDecoratorCfg(flags):
    result=ComponentAccumulator()

    from CaloCalibHitRec.CaloCalibHitDecoratorCfg import CaloCalibHitDecoratorCfg 
    result.merge(CaloCalibHitDecoratorCfg(flags))

    PFlowCalibPFODecoratorAlgorithm = CompFactory.PFlowCalibPFODecoratorAlgorithm()
    PFlowCalibPFODecoratorAlgorithm.TruthAttributerTool = CompFactory.CaloCalibClusterTruthAttributerTool("PFlowCalibPFOTruthAttributerTool")
    result.addEventAlgo(PFlowCalibPFODecoratorAlgorithm)

    return result
