#!/bin/sh
#
# art-description: Run MT and ST simulation outside ISF, reading ttbar events, writing HITS, using MC16 geometry and conditions

# art-include: master/Athena
# art-type: grid
# art-architecture:  '#x86_64-intel'
# art-athena-mt: 8
# art-output: log.*
# art-output: test.MT.HITS.pool.root
# art-output: test.ST.HITS.pool.root

export ATHENA_CORE_NUMBER=8

Sim_tf.py \
--multithreaded \
--inputEVNTFile "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.evgen.EVNT.e4993.EVNT.08166201._000012.pool.root.1" \
--outputHITSFile "test.MT.HITS.pool.root" \
--maxEvents 10 \
--geometryVersion 'default:ATLAS-R2-2016-01-00-01' \
--conditionsTag 'default:OFLCOND-MC16-SDR-14' \
--simulator 'ATLFAST3MTEnergyOrdered' \
--postInclude 'default:PyJobTransforms/UseFrontier.py' \
--preInclude 'EVNTtoHITS:Campaigns/MC16Simulation.py' #\
#--postExec 'EVNTtoHITS:ServiceMgr.MessageSvc.enableSuppression=False;topSeq.ISF_Kernel_ATLFAST3MT.OutputLevel=VERBOSE;' \
#--imf False

rc=$?
mv log.EVNTtoHITS log.ATLFAST3MTAthenaMT
echo  "art-result: $rc ATLFAST3MTAthenaMT"
rc=0
rc2=-9999
if [ $rc -eq 0 ]
then
    ArtPackage=$1
    ArtJobName=$2
    unset ATHENA_CORE_NUMBER
    Sim_tf.py \
    --inputEVNTFile "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.evgen.EVNT.e4993.EVNT.08166201._000012.pool.root.1" \
  --outputHITSFile "test.ST.HITS.pool.root" \
  --maxEvents 10 \
  --geometryVersion 'default:ATLAS-R2-2016-01-00-01' \
  --conditionsTag 'default:OFLCOND-MC16-SDR-14' \
  --simulator 'ATLFAST3MTEnergyOrdered' \
  --postInclude 'default:PyJobTransforms/UseFrontier.py' \
  --preInclude 'EVNTtoHITS:Campaigns/MC16Simulation.py' #\
#  --postExec 'EVNTtoHITS:ServiceMgr.MessageSvc.enableSuppression=False;topSeq.ISF_Kernel_ATLFAST3MT.OutputLevel=VERBOSE;' \
#  --imf False
  rc2=$?
  mv log.EVNTtoHITS log.ATLFAST3MTAthena
fi
echo  "art-result: $rc2 ATLFAST3MTAthena"
rc3=-9999
if [ $rc -eq 0 ]
then
    ArtPackage=$1
    ArtJobName=$2
    unset ATHENA_CORE_NUMBER
    Sim_tf.py \
    --inputEVNTFile "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.evgen.EVNT.e4993.EVNT.08166201._000012.pool.root.1" \
  --outputHITSFile "test.ST.old.HITS.pool.root" \
  --maxEvents 10 \
  --geometryVersion 'default:ATLAS-R2-2016-01-00-01' \
  --conditionsTag 'default:OFLCOND-MC16-SDR-14' \
  --simulator 'ATLFAST3EnergyOrdered' \
  --postInclude 'default:PyJobTransforms/UseFrontier.py' \
  --preInclude 'EVNTtoHITS:Campaigns/MC16Simulation.py' #\
#  --postExec 'EVNTtoHITS:ServiceMgr.MessageSvc.enableSuppression=False;topSeq.ISF_Kernel_ATLFAST3.OutputLevel=VERBOSE;ServiceMgr.ISF_AFIIParticleBrokerSvc.OutputLevel=VERBOSE;' \
#  --imf False
  rc3=$?
  mv log.EVNTtoHITS log.ATLFAST3Athena
fi
echo  "art-result: $rc3 ATLFAST3Athena"
rc4=-9999
if [ $rc2 -eq 0 ]
then
    acmd.py diff-root test.MT.HITS.pool.root test.ST.HITS.pool.root --error-mode=resilient --mode=semi-detailed --order-trees --ignore-leaves RecoTimingObj_p1_EVNTtoHITS_timings index_ref
    rc4=$?
fi
echo  "art-result: $rc4 ATLFAST3MT_STvsMT"
rc5=-9999
if [ $rc3 -eq 0 ]
then
    acmd.py diff-root test.ST.HITS.pool.root test.ST.old.HITS.pool.root --error-mode=resilient --mode=semi-detailed --order-trees --ignore-leaves RecoTimingObj_p1_EVNTtoHITS_timings index_ref
    rc5=$?
fi
echo  "art-result: $rc5 ATLFAST3MTvsATLFAST3"
