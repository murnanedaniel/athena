# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

"""
Configuration database for ISF_SimulationSelectors
Elmar Ritsch, 10/11/2014
"""

from AthenaCommon.CfgGetter import addTool, addToolClone, addService, addAlgorithm, \
     addTypesToExcludeIfDefaultValue, addNamesToExcludeIfDefaultValue, addFullNamesToExcludeIfDefaultValue, \
     addPropertiesToExcludeIfDefault, \
     addTypesToSkipIfNotAvailable, addNamesToSkipIfNotAvailable, addFullNamesToSkipIfNotAvailable, \
     addTypesOnlyToSkip

from AthenaCommon.Constants import *  # FATAL,ERROR etc.
import AthenaCommon.SystemOfUnits as Units


addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getElectronGeant4Selector"               , "ISF_ElectronGeant4Selector"              )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getNeutralGeant4Selector"                , "ISF_NeutralGeant4Selector"               )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getMuonGeant4Selector"                   , "ISF_MuonGeant4Selector"                  )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getMuonAFIIGeant4Selector"               , "ISF_MuonAFIIGeant4Selector"              )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getMuonFatrasSelector"                   , "ISF_MuonFatrasSelector"                  )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getPhotonConeFatrasSelector"             , "ISF_PhotonConeFatrasSelector"            )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getPhotonConeGeant4Selector"             , "ISF_PhotonConeGeant4Selector"            )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getBHadronProductsGeant4Selector"        , "ISF_BHadronProductsGeant4Selector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getBHadronProductsFatrasSelector"        , "ISF_BHadronProductsFatrasSelector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getTauProductsGeant4Selector"            , "ISF_TauProductsGeant4Selector"           )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getZProductsSimSelector"                 , "ISF_ZProductsSimSelector"                )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getZProductsGeant4Selector"              , "ISF_ZProductsGeant4Selector"             )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getHiggsLeptonsConeGeant4Selector"       , "ISF_HiggsLeptonsConeGeant4Selector"      )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getHiggsLeptonsConeGeant4CaloSelector"   , "ISF_HiggsLeptonsConeGeant4CaloSelector"  )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getWLeptonsConeGeant4Selector"           , "ISF_WLeptonsConeGeant4Selector"          )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getZLeptonsDirectionConeGeant4Selector"  , "ISF_ZLeptonsDirectionConeGeant4Selector" )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getZLeptonsPositionConeGeant4Selector"   , "ISF_ZLeptonsPositionConeGeant4Selector"  )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getJPsiLeptonsConeGeant4Selector"        , "ISF_JPsiLeptonsConeGeant4Selector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getWithinEta5FastCaloSimSelector"        , "ISF_WithinEta5FastCaloSimSelector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getEtaGreater5ParticleKillerSimSelector" , "ISF_EtaGreater5ParticleKillerSimSelector")
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getSubDetStickyGeant4SimSelector"        , "ISF_SubDetStickyGeant4SimSelector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getGlobalStickyGeant4SimSelector"        , "ISF_GlobalStickyGeant4SimSelector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultFastCaloSimSelector"           , "ISF_DefaultFastCaloSimSelector"          )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFastCaloSimPileupSelector"            , "ISF_FastCaloSimPileupSelector"           )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFastHitConvAlgFastCaloSimSelector"    , "ISF_FastHitConvAlgFastCaloSimSelector"   )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultLegacyAFIIFastCaloSimSelector" , "ISF_DefaultLegacyAFIIFastCaloSimSelector")
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFastHitConvAlgLegacyAFIIFastCaloSimSelector" , "ISF_FastHitConvAlgLegacyAFIIFastCaloSimSelector")
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultFatrasSelector"                , "ISF_DefaultFatrasSelector"               )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultFatrasNewExtrapolationSelector", "ISF_DefaultFatrasNewExtrapolationSelector")
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultParticleKillerSelector"        , "ISF_DefaultParticleKillerSelector"       )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultGeant4Selector"                , "ISF_DefaultGeant4Selector"               )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultAFIIGeant4Selector"            , "ISF_DefaultAFIIGeant4Selector"           )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultLongLivedGeant4Selector"       , "ISF_DefaultLongLivedGeant4Selector"      )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFullGeant4Selector"                   , "ISF_FullGeant4Selector"                  )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getPassBackGeant4Selector"               , "ISF_PassBackGeant4Selector"              )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFatrasPileupSelector"                 , "ISF_FatrasPileupSelector"                )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFatrasPileupSelector_noHits"          , "ISF_FatrasPileupSelector_noHits"         )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getFatrasRandomSelector"                 , "ISF_FatrasRandomSelector"                )
addTool("ISF_SimulationSelectors.ISF_SimulationSelectorsConfig.getDefaultParametricSimulationSelector"  , "ISF_DefaultParametricSimulationSelector" )
