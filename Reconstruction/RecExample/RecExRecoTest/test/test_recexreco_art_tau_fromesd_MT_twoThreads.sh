#!/bin/sh
#
# art-description: Athena runs topoclustering from an ESD file
# art-type: grid
# art-include: master/Athena

athena --threads=2 tauRec/run_tau_standalone.py
echo "art-result: $?"
