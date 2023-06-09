#!/bin/sh
set -e

Reco_tf.py  --outputHISTFile=myMergedMonitoring_IDCosmic_0.root --maxEvents=500 --outputESDFile=myESD_IDCosmic_0.pool.root --outputAODFile=myAOD_IDCosmic_0.AOD.pool.root --outputTAGFile=myTAG_IDCosmic_0.root --ignoreErrors=False --inputBSFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data11_7TeV.00184074.physics_IDCosmic.merge.RAW._lb0100._SFO-ALL._0001.1,/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data11_7TeV.00184074.physics_IDCosmic.merge.RAW._lb0101._SFO-ALL._0001.1,/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data11_7TeV.00184074.physics_IDCosmic.merge.RAW._lb0102._SFO-ALL._0001.1 --preExec='rec.doTrigger=False;'

