#!/bin/sh

# art-include: master/AthDerivation
# art-include: master/Athena
# art-description: DAOD building EGAM2 EGAM3 EGAM4 data18
# art-type: grid
# art-output: *.pool.root

set -e

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data18_13TeV.00357772.physics_Main.recon.AOD.r13286/AOD.27654050._000557.pool.root.1 --outputDAODFile art.pool.root --reductionConf EGAM2 EGAM3 EGAM4 --maxEvents 500  --preExec 'from AthenaCommon.DetFlags import DetFlags; DetFlags.detdescr.all_setOff(); DetFlags.BField_setOn(); DetFlags.digitize.all_setOff(); DetFlags.detdescr.Calo_setOn(); DetFlags.simulate.all_setOff(); DetFlags.pileup.all_setOff(); DetFlags.overlay.all_setOff(); DetFlags.detdescr.Muon_setOn();' --postExec 'from DerivationFrameworkJetEtMiss.JetCommon import swapAlgsInSequence; swapAlgsInSequence(topSequence,"jetalg_ConstitModCorrectPFOCSSKCHS_GPFlowCSSK", "UFOInfoAlgCSSK" );'
