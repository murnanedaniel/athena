#!/bin/sh
#
# art-description: MC16-style simulation of decaying staus using FullG4 (tests the Sleptons + Gauginos packages)
# art-architecture:  '#x86_64-intel'
# art-type: grid
# art-output: *.root
# art-output: PDGTABLE.*
# art-output: log.*
# art-include: 21.0/Athena
# art-include: 21.0/AthSimulation
# art-include: 21.3/Athena
# art-include: 21.9/Athena
# art-include: master/Athena
# art-include: master/AthSimulation

# MC16 setup
# ATLAS-R2-2016-01-00-01 and OFLCOND-MC16-SDR-14
Sim_tf.py \
    --conditionsTag 'default:OFLCOND-MC16-SDR-14' \
    --physicsList 'FTFP_BERT_ATL' \
    --truthStrategy 'MC15aPlusLLP' \
    --simulator 'FullG4' \
    --postExec 'EVNTtoHITS:ServiceMgr.ISF_InputConverter.GenParticleFilters["ISF_GenParticleInteractingFilter"].AdditionalNonInteractingParticleTypes=[-51,51,52,-53,53]' \
    --postInclude 'default:PyJobTransforms/UseFrontier.py' \
    --preInclude 'EVNTtoHITS:SimulationJobOptions/preInclude.BeamPipeKill.py,SimulationJobOptions/preInclude.FrozenShowersFCalOnly.py' \
    --preExec 'EVNTtoHITS:simFlags.TightMuonStepping=True; simFlags.SimBarcodeOffset.set_Value_and_Lock(200000); simFlags.TRTRangeCut=30.0;' \
    --DataRunNumber '284500' \
    --geometryVersion 'default:ATLAS-R2-2016-01-00-01' \
    --inputEVNTFile "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/mc16_13TeV.999999.DMsimp_s_spin1_SVJ.test.EVNT.pool.root" \
    --outputHITSFile "Hits.pool.root" \
    --maxEvents 10 \
    --imf False

rc=$?
echo  "art-result: $rc simulation"
rc2=-9999
if [ $rc -eq 0 ]
then
    ArtPackage=$1
    ArtJobName=$2
    art.py compare grid --entries 10 ${ArtPackage} ${ArtJobName} --mode=semi-detailed
    rc2=$?
fi
echo  "art-result: $rc2 regression"
