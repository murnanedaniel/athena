/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef __TRIGEGAMMAMONTOOL_H
#define __TRIGEGAMMAMONTOOL_H
#include "TrigHLTMonitoring/IHLTMonTool.h"
class ITrigEgammaAnalysisBaseTool;
class TrigEgammaMonTool : public IHLTMonTool {

    public:
        TrigEgammaMonTool( const std::string & type, const std::string & name, const IInterface* parent);
        virtual ~TrigEgammaMonTool();
        virtual StatusCode init();
        virtual StatusCode book();
        virtual StatusCode fill();
        virtual StatusCode proc();
    private:
        ToolHandleArray< ITrigEgammaAnalysisBaseTool > m_asgtools;
        std::vector<std::string> m_asgToolNames;
        SG::ReadHandleKey<xAOD::JetContainer> m_jetContainerKey{ this, "JetContainerKey", "AntiKt4EMTopoJets", "jet container for OLR" };
    protected:

};
#endif //__TRIGEGAMMAANALYSIST0BASETOOL_H

