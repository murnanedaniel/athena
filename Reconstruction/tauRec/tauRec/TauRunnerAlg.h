/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TAUREC_TAURUNNERALG_H
#define TAUREC_TAURUNNERALG_H

#include "GaudiKernel/ToolHandle.h"
#include "AthenaBaseComps/AthAlgorithm.h"
#include "tauRecTools/ITauToolExecBase.h"
#include "tauRecTools/ITauToolBase.h"

#include "tauRecTools/TauEventData.h"

#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteHandle.h"

#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/PFOAuxContainer.h"

#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODCaloEvent/CaloClusterAuxContainer.h"


/**
 * @brief       Main class for tau candidate processing.
 */

class TauRunnerAlg: public AthAlgorithm
{
    public:
        //-----------------------------------------------------------------
        // Contructor and destructor
        //-----------------------------------------------------------------
        TauRunnerAlg( const std::string &name, ISvcLocator *pSvcLocator );
        ~TauRunnerAlg();

        //-----------------------------------------------------------------
        // Gaudi algorithm hooks
        //-----------------------------------------------------------------
        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private:
       
	ToolHandleArray<ITauToolBase>  m_tools;

        //ToolHandleArray<ITauToolExecBase>  m_tools;
	TauEventData m_data;

	SG::ReadHandleKey<xAOD::TauJetContainer> m_tauInputContainer{this,"Key_tauInputContainer","tmp_TauJets","input temp tau key"};
	SG::WriteHandleKey<xAOD::TauJetContainer> m_tauOutputContainer{this,"Key_tauOutputContainer","TauJets","output tau data key"};	
	
};

#endif // TAUREC_TAURUNNERALG_H
