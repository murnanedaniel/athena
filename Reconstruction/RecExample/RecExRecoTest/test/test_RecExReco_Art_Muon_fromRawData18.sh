#!/bin/sh
#
# art-description: Athena runs muon reconstruction from a RAW data18 file
# art-type: grid
# art-athena-mt: 8
# art-include: master/Athena
# art-include: 22.0/Athena
# art-output: *.log   

python $Athena_DIR/python/RecExRecoTest/MuonReco_RAWData18.py | tee temp.log
echo "art-result: ${PIPESTATUS[0]}"
test_postProcessing_Errors.sh temp.log

