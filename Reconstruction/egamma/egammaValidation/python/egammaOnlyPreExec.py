# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from RecExConfig.RecFlags import rec
from RecExConfig.RecAlgsFlags import recAlgs
from ParticleBuilderOptions.AODFlags import AODFlags
from JetRec.JetRecFlags import jetFlags
from AthenaConfiguration.AllConfigFlags import ConfigFlags


def setRunEgammaOnlyRecoFlags():

    ConfigFlags.Detector.GeometryMuon = False
    ConfigFlags.Reco.EnablePFlow = False

    rec.doTrigger.set_Value_and_Lock(False)
    rec.doTau.set_Value_and_Lock(False)
    rec.doMuon.set_Value_and_Lock(False)
    rec.doBTagging.set_Value_and_Lock(False)
    rec.doCaloRinger.set_Value_and_Lock(False)
    rec.doAlfa.set_Value_and_Lock(False)
    rec.doAFP.set_Value_and_Lock(False)
    rec.doLucid.set_Value_and_Lock(False)

    recAlgs.doMuonSpShower.set_Value_and_Lock(False)
    recAlgs.doEFlow.set_Value_and_Lock(False)
    recAlgs.doEFlowJet.set_Value_and_Lock(False)
    recAlgs.doMissingET.set_Value_and_Lock(False)
    recAlgs.doMissingETSig.set_Value_and_Lock(False)

    AODFlags.ThinGeantTruth.set_Value_and_Lock(False)
    AODFlags.egammaTrackSlimmer.set_Value_and_Lock(False)
    AODFlags.ThinTRTStandaloneTracks.set_Value_and_Lock(False)
    AODFlags.ThinInDetForwardTrackParticles.set_Value_and_Lock(False)
    AODFlags.ThinNegativeEnergyNeutralPFOs.set_Value_and_Lock(False)

    jetFlags.Enabled = False

def setRunEgammaOnlyMergeFlags():
    rec.doTrigger.set_Value_and_Lock(False)
