#Matt A's custom python accumulator for R22 
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from pprint import pprint

#Trigger xAODType.ObjectType dict entry loading
import cppyy
try:
    cppyy.load_library('libxAODBaseObjectTypeDict')
except Exception:
    pass
from ROOT import xAODType
xAODType.ObjectType


def config_CHS_CSSK_merged(inputFlags,**kwargs):
    #
    import cppyy
    try:
        cppyy.load_library('libxAODBaseObjectTypeDict')
    except Exception:
        pass
    from ROOT import xAODType
    xAODType.ObjectType

    output_CA=ComponentAccumulator()
    StoreGateSvc=CompFactory.StoreGateSvc
    output_CA.addService(StoreGateSvc("DetectorStore"))

    from JetRecConfig.JetRecConfig import JetInputCfg
    from JetRecConfig.StandardJetConstits import stdConstitDic as cst
    

    output_CA.merge(JetInputCfg(inputFlags, cst.EMPFlowCSSK))
    return output_CA
    
    
                       



def UFOConfig(inputFlags, **kwargs):
    # Configure UFO reconstruction algorithm instead of TCC
    UFO_CA=ComponentAccumulator()
    #add Storegate
    StoreGateSvc=CompFactory.StoreGateSvc
    
    UFO_CA.addService(StoreGateSvc("DetectorStore"))
    from TrackCaloClusterRecTools.TrackCaloClusterConfig import runUFOReconstruction
    from JetRecConfig.StandardJetConstits import stdConstitDic as cst
    constituents=cst.EMPFlowCSSK
    inputFEcontainer="CSSKParticleFlowObjects"
    
    UFO_reco=runUFOReconstruction(configFlags=inputFlags,constits=constituents,inputFEcontainerkey=inputFEcontainer)
    UFO_CA.merge(UFO_reco)

    return UFO_CA
    


if __name__=="__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags as cfgFlags
    
    cfgFlags.Input.isMC=True
    
    #    cfgFlags.Input.Format="AOD"
    
    #    cfgFlags.Input.Files=["/scratch/anthony/TEST_AOD/TCC/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r10724/mc16_13TeV/AOD.14795494._007454.pool.root.1"]
    #    cfgFlags.Input.Files=["/scratch/anthony/TEST_AOD/TCC/valid1_test/valid1/AOD.24855256._000001.pool.root.1"]
    cfgFlags.Input.Files=["myAOD.root"]
    cfgFlags.Output.AODFileName="output_UFO_DAOD.root"
    cfgFlags.Exec.MaxEvents=20
    cfgFlags.Output.doWriteAOD=True
    #cfgFlags.Common.ProductionStep=0 # I think this is default
    cfgFlags.dump()
    cfgFlags.lock()
    

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg=MainServicesCfg(cfgFlags)

    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    inputList=["xAOD::TrackParticleContainer#InDetTrackParticles", "xAOD::TrackParticleAuxContainer#InDetTrackParticlesAux."]
    inputList.append("xAOD::FlowElementContainer#JetETMissChargedParticleFlowObjects")
    inputList.append("xAOD::FlowElementContainer#JetETMissNeutralParticleFlowObjects")
    inputList.append("xAOD::FlowElementAuxContainer#JetETMissNeutralParticleFlowObjectsAux.")
    inputList.append("xAOD::FlowElementAuxContainer#JetETMissChargedParticleFlowObjectsAux.")
    inputList.append("xAOD::FlowElementContainer#TrackCaloClustersCharged")
    inputList.append("xAOD::FlowElementAuxContainer#TrackCaloClustersChargedAux.")
    inputList.append("xAOD::FlowElementContainer#TrackCaloClustersNeutral")
    inputList.append("xAOD::FlowElementAuxContainer#TrackCaloClustersNeutralAux.")
    inputList.append("xAOD::FlowElementContainer#UFO*")
    inputList.append("xAOD::FlowElementAuxContainer#UFO*")
    inputList.append("xAOD::FlowElementContainer#CHS*")
    inputList.append("xAOD::FlowElementAuxContainer#CHS*")
    inputList.append("xAOD::FlowElementContainer#CSSK*")
    inputList.append("xAOD::FlowElementAuxContainer#CSSK*")
    inputList.append("xAOD::FlowElementContainer#*")
    inputList.append("xAOD::FlowElementAuxContainer#*")
    inputList.append("xAOD::MuonContainer#Muons")
    inputList.append("xAOD::MuonAuxContainer#*")
    inputList.append("xAOD::PhotonContainer#Photons")
    inputList.append("xAOD::PhotonAuxContainer#*")
    inputList.append("xAOD::ElectronContainer#Electrons")
    inputList.append("xAOD::ElectronAuxContainer#*")
    inputList.append("xAOD::TauJetContainer#*")
    inputList.append("xAOD::TauJetAuxContainer#*")
    
    
    cfg.merge(OutputStreamCfg(cfgFlags,"AOD",ItemList=inputList))
    
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg

    
    cfg.merge(PoolReadCfg(cfgFlags))

    cfg.merge(config_CHS_CSSK_merged(cfgFlags))

    cfg.merge(UFOConfig(cfgFlags)) 

    cfg.printConfig()

    cfg.run()

