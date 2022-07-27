# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from TriggerMenuMT.HLT.Config.Utility.ChainDefInMenu import ChainProp
from .SignatureDicts import ChainStore

def setupMenu():

    from AthenaCommon.Logging                 import logging
    log = logging.getLogger( __name__ )
    log.info('Executing menu....')

    chains = ChainStore()

    chains['Test'] = [
        # muons
        ChainProp(name='HLT_TestChain8_mv1step_L1MU5VF', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain8_mv1_L1MU8F', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain20_mv1_L1MU8F',stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain10_mv2_L1MU8F',stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain6_mEmpty2_L1MU5VF',  stream=['Main'], groups=['RATE:Test','BW:Other'] ),

        # egamma
        ChainProp(name='HLT_TestChain5_ev1_L1EM3', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain8_ev1_L1EM3', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain5_ev2_L1EM7', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain5_ev3_L1EM7', stream=['Main'], groups=['RATE:Test','BW:Other'] ),
        ChainProp(name='HLT_TestChain5_gv1_L1EM7', stream=['Main'], groups=['RATE:Test','BW:Other'] ),

        # combined
        ChainProp(name='HLT_TestChain6_mv1_TestChain10_ev1_L1EM7_MU8F', l1SeedThresholds=['MU8F','EM7'], stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel'), 
        ChainProp(name='HLT_TestChain6_mv2_TestChain8_ev2_L1EM7_MU8F',  l1SeedThresholds=['MU8F','EM7'], stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel'),
        ChainProp(name='HLT_TestChain5_ev1_TestChain8_ev1_L12EM3',      stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel'),
#        ChainProp(name='HLT_TestChain5_ev1_TestChain8_ev1_2TestChain6_mv1_L1EM3_L1EM3_L12MU5VF', stream=['Main'], groups=['RATE:Test','BW:Other']   ),
        ChainProp(name='HLT_2TestChain6_mv1_L12MU5VF',                   stream=['Main'], groups=['RATE:Test','BW:Other']   ),
        ChainProp(name='HLT_3TestChain6_mv1_L12MU5VF',                   stream=['Main'], groups=['RATE:Test','BW:Other']   ),
        ChainProp(name='HLT_TestChain6_mv1_TestChain10_mv1_L12MU5VF',   stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel' ),        
        ChainProp(name='HLT_2TestChain6_mEmpty1_L12MU5VF',               stream=['Main'], groups=['RATE:Test','BW:Other']), #may differ from manual
        ChainProp(name='HLT_TestChain6_mv1_TestChain5_ev1dr_L1EM7_MU8F', l1SeedThresholds=['MU8F','EM7'], stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel'), 
        ChainProp(name='HLT_2TestChain4_mv1dr_L12MU5VF', stream=['Main'], groups=['RATE:Test','BW:Other'] ),

        # FSNOSEED not implemented in emulation
        ChainProp(name='HLT_TestChain10_mEmpty1_TestChain6_mEmpty1_L12MU5VF',  l1SeedThresholds=['MU5VF','MU5VF'],  stream=['Main'], groups=['RATE:Test','BW:Other'], mergingStrategy = 'parallel'),
    ]

    return chains
