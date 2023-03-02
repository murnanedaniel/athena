#!/bin/sh
#newgrp l2it
local_dir=$(pwd)

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
PS1="\[\033[1;32m\]\u\[\033[0m\]\[\033[32m\][\w]\[\033[0m\]> "
lsetup git
asetup 21.9,Athena,21.9.26
git clone https://gitlab.cern.ch/gnn4itkteam/athena.git
##     [This command works because I have created my own clone of the athena git repository. One does that once and only once.]

cd athena
git remote add upstream https://:@gitlab.cern.ch:8443/atlas/athena.git
git fetch upstream
git checkout 21.9.26-root-and-csv-files-from-RDO
#     [Now we have a local copy of the entire Athena code, version 21.9.26] that includes the dumper

cd ${local_dir}/athena
#     [i.e. return to the « head » directory of the local copy of athena]

mkdir build
cd build
cmake -DATLAS_PACKAGE_FILTER_FILE=${local_dir}/athena/Tracking/TrkDumpAlgs/package_filters.txt ${local_dir}/athena/Projects/WorkDir
#     [This configures a « sparse build » : only the packages that are listed in package_filter.txt will be recompiled, the rest is taken from the official ATLAS release.]

make
#     [This is the actual build]

source x86_64-centos7-gcc62-opt/setup.sh 
#     [activate our sparse build for execution]

cd ../../
sh athena/Tracking/TrkDumpAlgs/scripts/run_reco.sh

