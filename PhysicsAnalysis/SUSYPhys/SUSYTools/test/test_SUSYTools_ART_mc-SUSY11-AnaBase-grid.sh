#!/bin/sh

# art-description: SUSYTools ART test - TestSUSYToolsAlg.py
# art-type: grid
# art-include: 21.2/AnalysisBase
# art-output: hist-MC_DAOD_SUSY11.root
# art-output: dcube

# Create empty pool file
art.py createpoolfile

echo "Running SUSYTools test: \'TestSUSYToolsAlg.py -t MCe -m -1 -d 11\'"
TestSUSYToolsAlg.py -t MCe -m -1 -d 11
echo  "art-result: $? TEST"

mv submitDir/hist-ARTInput.root ./hist-MC_DAOD_SUSY11.root

echo "Running DCube post-processing"

tName="MCe-SUSY11"
dcubeRef=/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ART/References/hist-MC_DAOD_SUSY11.root
dcubeXml=/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ART/xml/dcube_config_master.xml

/cvmfs/atlas.cern.ch/repo/sw/art/dcube/bin/art-dcube ${tName} hist-MC_DAOD_SUSY11.root ${dcubeXml} ${dcubeRef}

echo "art-result: $? DCUBE"

echo "Done."
