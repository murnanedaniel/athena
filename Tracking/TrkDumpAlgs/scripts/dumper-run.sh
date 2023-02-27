local_dir=$(pwd)
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup git
asetup 21.9,Athena,21.9.26
source athena/build/x86_64-centos7-gcc62-opt/setup.sh
export PATH=$PATH:${local_dir}/athena/Reconstruction/RecJobTransforms/scripts/Reco_tf.py
bash ${local_dir}/athena/Tracking/TrkDumpAlgs/scripts/run_reco.sh
