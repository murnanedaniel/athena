#!/bin/bash

if [ -z ${TEST} ]; then
  export TEST="TrigP1Test"
fi

if [ -z ${NAME} ]; then
  export NAME="UNDEFINED"
fi

if [ -z ${JOB_LOG} ]; then
  export JOB_LOG="atn_test.log"
fi

if [ -z ${JOB_athenaHLT_LOG} ]; then
  export JOB_athenaHLT_LOG="atn_test_athenaHLT.log"
fi

if [ -z ${JOB_LOG_LVL} ]; then
  export JOB_LOG_LVL="INFO,ERROR"
  #export TDAQ_ERS_DEBUG_LEVEL=3
fi

# read test configurations from env variables

if [ -z "${ART_CMD}" ]; then
  echo "Please provide test command as \$ART_CMD" 2>&1 | tee -a ${JOB_LOG}
  echo "echo 'art-result: 999 ${NAME}.athena_mother'"  2>&1 | tee -a ${JOB_LOG}
  exit 1
fi

if [ -z ${ART_TIMEOUT} ]; then
  export ART_TIMEOUT="200m"
fi

if [ -z ${ART_FILE_NAME} ]; then
  export ART_FILE_NAME="/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/TrigP1Test/data18_13TeV.00360026.physics_EnhancedBias.MissingTowers._lb0151._SFO-6._0001.1.pool.root"
fi

ART_CMD="timeout ${ART_TIMEOUT} ${ART_CMD} --log-level ${JOB_LOG_LVL}"
ART_CMD=${ART_CMD/INPUT_FILE/$ART_FILE_NAME}
echo "----------------------------------------" 2>&1 | tee -a ${JOB_LOG}
echo "Running athenaHLT command:" 2>&1 | tee -a ${JOB_LOG}
echo ${ART_CMD} 2>&1 | tee -a ${JOB_LOG}
echo "#!/bin/bash" >> run_athenaHLT.sh
echo ${ART_CMD} >> run_athenaHLT.sh
echo "echo 'art-result: '\$?' ${NAME}.athena_mother'" >> run_athenaHLT.sh
printenv >> printenv.log
source run_athenaHLT.sh 2>&1 | tee -a ${JOB_athenaHLT_LOG}
echo "----------------------------------------" 2>&1 | tee -a ${JOB_LOG}

source exec_art_trigp1test_merge.sh 2>&1 | tee -a ${JOB_LOG}
source exec_art_trigp1test_post.sh 2>&1 | tee -a ${JOB_LOG}
source exec_art_trigp1test_summary.sh 2>&1 | tee -a ${JOB_LOG}

