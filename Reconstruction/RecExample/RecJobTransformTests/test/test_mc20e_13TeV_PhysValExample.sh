#!/bin/sh
#
# art-description: Run PhysicaValidation Example code on mc20e AOD. Uses AOD from test_mc20e_nopileup_13TeV.sh on 27 May 2021
# art-output: log.*
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena

export TRF_ECHO=True;
Reco_tf.py --inputAODFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/RecJobTransformTests/mc16_13TeV/valid1/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.recon.AOD.e4993_s3227_r12649/myAOD.pool.root  --outputNTUP_PHYSVALFile physval.root  --validationFlags doExample --conditionsTag="OFLCOND-MC16-SDR-RUN2-09"
RES=$?
echo "art-result: $RES Reco"
