#! /usr/bin/env bash


DCSC2_SYSTEMS="-sMDT -sTile -sTGC -sRPC -sTDQ -sMagnets -sGlobal -sTRT -sSCT -sLAr -sLucid -sTrigger -sAFP"
#DCSC2_SYSTEMS="-sMDT -sTile -sTGC -sRPC -sTDQ -sCSC -sMagnets -sGlobal -sPixels -sTRT -sLAr -sLucid" #outdated
#DEST_DB=COOLOFL_GLOBAL/COMP200 #outdated
DEST_DB=COOLOFL_GLOBAL/CONDBR2

RUN=$1
shift

if [ -z $RUN ]; then 
    echo "Usage: ./ExecuteDCSC2.sh [run_number]"
    exit 1
fi

echo "Running for $RUN"

export AtlasSetup=/afs/cern.ch/atlas/software/dist/AtlasSetup
pushd /afs/cern.ch/user/a/atlasdqm/ws/DCSCalc/prodarea > /dev/null
source $AtlasSetup/scripts/asetup.sh 22.0.72,Athena
source /afs/cern.ch/user/a/atlasdqm/DQCalculators/DCSCalc/prodarea/build/x86_64-centos7-gcc11-opt/setup.sh

CORAL_AUTH_PATH=/afs/cern.ch/user/a/atlasdqm/private
CORAL_DBLOOKUP_PATH=/afs/cern.ch/user/a/atlasdqm/private
#export FRONTIER_LOG_LEVEL=debug

#dcsc.py -h
dcsc.py $@ $DCSC2_SYSTEMS -r$RUN -o$DEST_DB --email-on-failure
#dcsc.py $@ $DCSC2_SYSTEMS -r$RUN -o$DEST_DB
