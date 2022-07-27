# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from math import inf

# this function is only used below
def JetParticleAssociationCfg(ConfigFlags, jetCollName, partcollname, assocname, **options):

    acc=ComponentAccumulator()

    options["coneSizeFitPar1"] = +0.239
    options["coneSizeFitPar2"] = -1.220
    options["coneSizeFitPar3"] = -1.64e-5
    options["InputParticleContainer"] = partcollname
    options["OutputDecoration"] = assocname

    # -- create the association tool
    acc.setPrivateTools(
        CompFactory.JetParticleShrinkingConeAssociation(
            JetContainer=jetCollName, **options))

    return acc


def JetParticleAssociationAlgCfg(
        ConfigFlags,
        JetCollection,
        InputParticleCollection,
        OutputParticleDecoration,
        MinimumJetPt=None,
        MinimumJetPtFlag=None):

    acc=ComponentAccumulator()
    jetcol = JetCollection.replace("Track", "PV0Track")
    name=(jetcol + "_" + OutputParticleDecoration + "_assoc").lower()
    if MinimumJetPt is None:
        MinimumJetPt = ConfigFlags.BTagging.minimumJetPtForTrackAssociation
    if MinimumJetPt > 0.0 and MinimumJetPtFlag is None:
        ptflag = f'{OutputParticleDecoration}OverPtThreshold'
    elif MinimumJetPtFlag is not None:
        ptflag = MinimumJetPtFlag
    else:
        ptflag = ''

    # -- create the association algorithm
    acc.addEventAlgo(CompFactory.JetDecorationAlg(
        name=name,
        JetContainer=jetcol,
        Decorators=[
            acc.popToolsAndMerge(
                JetParticleAssociationCfg(
                    ConfigFlags,
                    jetcol,
                    InputParticleCollection,
                    OutputParticleDecoration,
                    MinimumJetPt=MinimumJetPt,
                    PassPtFlag=ptflag,
                ))
        ]
    ))

    return acc

def JetParticleFixedConeAssociationAlgCfg(ConfigFlags, fixedConeRadius, JetCollection="", InputParticleCollection="", OutputParticleDecoration="", **options):

    acc=ComponentAccumulator()
    jetcol = JetCollection.replace("Track", "PV0Track")

    options['JetContainer'] = jetcol
    options['Decorators'] = [CompFactory.JetParticleShrinkingConeAssociation(
                                                                             InputParticleContainer=InputParticleCollection,
                                                                             OutputDecoration=OutputParticleDecoration,
                                                                             coneSizeFitPar1=fixedConeRadius,
                                                                             coneSizeFitPar2=-inf,
                                                                             coneSizeFitPar3=0,
                                                                             **options)]
    options['name'] = (jetcol + "_" + OutputParticleDecoration + "_assoc").lower()

    # -- create the association algorithm
    acc.addEventAlgo(CompFactory.JetDecorationAlg(**options))

    return acc
