#!/bin/sh

# art-include: master/AthDerivation
# art-include: master/Athena
# art-description: DAOD building FTAG1 JETM1 mc20
# art-type: grid
# art-output: *.pool.root

set -e

Reco_tf.py --inputAODFile /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/DerivationFrameworkART/mc20_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3681_r13167/AOD.27162646._000001.pool.root.1 --outputDAODFile art.pool.root --reductionConf FTAG1 JETM1 --maxEvents 500  --preExec 'from AthenaCommon.DetFlags import DetFlags; DetFlags.detdescr.all_setOff(); DetFlags.BField_setOn(); DetFlags.digitize.all_setOff(); DetFlags.detdescr.Calo_setOn(); DetFlags.simulate.all_setOff(); DetFlags.pileup.all_setOff(); DetFlags.overlay.all_setOff(); DetFlags.detdescr.Muon_setOn();' --postExec 'from DerivationFrameworkJetEtMiss.JetCommon import swapAlgsInSequence; swapAlgsInSequence(topSequence,"jetalg_ConstitModCorrectPFOCSSKCHS_GPFlowCSSK", "UFOInfoAlgCSSK" );' --passThrough True
