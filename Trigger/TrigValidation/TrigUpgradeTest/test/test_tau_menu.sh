#!/bin/sh
# art-description: athenaMT trigger test using the tau menu from TrigUpgradeTest job options
# art-type: build
# art-include: master/Athena
# art-output: *.log
# art-output: *.root
# art-output: *.regtest.new

export NAME="TrigUpgradeTest_tau_menu"
export SKIPEVENTS=10
export EVENTS=20
export THREADS=1
export SLOTS=1
export JOBOPTION="TrigUpgradeTest/full_menu.py"
export REGTESTEXP="TriggerSummaryStep.*HLT_.*|TriggerMonitorFinal.*HLT_.*|TrigSignatureMoniMT.*HLT_.*"
export EXTRA="doEmptyMenu=True;doTauSlice=True"

# Find the regtest reference installed with the release
export REGTESTREF=`find_data.py TrigUpgradeTest/tau_menu.ref`

source exec_TrigUpgradeTest_art_athenaMT.sh
source exec_TrigUpgradeTest_art_post.sh

