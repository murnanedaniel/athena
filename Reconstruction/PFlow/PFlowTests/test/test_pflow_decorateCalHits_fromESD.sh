#!/bin/sh
#
# art-description: Athena runs pflow calibration hit decoration reconstruction
# art-type: grid
# art-include: master/Athena
# art-include: 23.0/Athena
# art-athena-mt: 8
# art-output: AOD.pool.root
# art-output: log_FE_validation.txt

export ATHENA_CORE_NUMBER=8 # set number of cores used in multithread to 8.

python $Athena_DIR/python/PFlowUtils/PFRunCalibHitDecorator_ESDtoAOD_mc20e.py
rc1=$?

echo "art-result: ${rc1}"
