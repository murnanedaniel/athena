export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup git
asetup 21.9,Athena,21.9.26
#mkdir build
cd athena
source build/x86_64-centos7-gcc62-opt/setup.sh
cd build
cmake -DATLAS_PACKAGE_FILTER_FILE=${local_dir}/athena/Tracking/TrkDumpAlgs/package_filters.txt ${local_dir}/athena/Projects/WorkDir
make
