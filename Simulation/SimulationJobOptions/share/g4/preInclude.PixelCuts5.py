#########################################################
#
# SimulationJobOptions/postInitOptions.PixelCuts5.py
# Zach Marshall
#
# For a special production to evaluate the effect of low
# energy secondaries on Pixel Clustering.
#
#########################################################

atlasG4log.info("G4 PIX Config: Setting PIX cut")

from AthenaCommon.CfgGetter import getPublicTool
pixelRegionTool = getPublicTool('PixelPhysicsRegionTool')
pixelRegionTool.ElectronCut = 0.005
pixelRegionTool.PositronCut = 0.005
pixelRegionTool.GammaCut = 0.005
