#!/bin/bash

# art-description: DirectIOART AthenaMP inputFile:RDO protocol=DAVS
# art-type: grid
# art-output: *.pool.root
# art-include: 21.0/Athena

set -e

export ATHENA_PROC_NUMBER=2
Reco_tf.py --AMI q221 --inputRDOFile davs://lcg-lrz-http.grid.lrz.de:443/pnfs/lrz-muenchen.de/data/atlas/dq2/atlasdatadisk/rucio/mc16_13TeV/90/96/RDO.11373415._000001.pool.root.1 --outputRDO_TRIGFile art.pool.root

echo "art-result: $? DirectIOART_AthenaMP_RecoTF_inputRDO_protocol_DAVS"