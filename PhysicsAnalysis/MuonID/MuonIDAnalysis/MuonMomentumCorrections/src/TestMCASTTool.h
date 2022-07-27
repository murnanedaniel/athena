/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MCAST_TOOLTESTER_H
#define MCAST_TOOLTESTER_H

// Gaudi/Athena include(s):
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// Local include(s):
#include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "PATInterfaces/SystematicsUtil.h"
#include "TrackInfo.h"

// Root include(s)
#include "TFile.h"
#include "TTree.h"

namespace CP {

    class TestMCASTTool : public AthAlgorithm {
    public:
        //::: Regular Algorithm constructor
        TestMCASTTool(const std::string& name, ISvcLocator* svcLoc);
        //::: Function initialising the algorithm
        virtual StatusCode initialize();
        //::: Function executing the algorithm
        virtual StatusCode execute();
        //::: Function finalizing the algoritm
        virtual StatusCode finalize();

    private:
        //::: Name of the output file
        std::string m_Output;
        //::: StoreGate key for the muon container to investigate
        std::string m_sgKey;
        //::: Connection to the smearing tool
        ToolHandle<CP::IMuonCalibrationAndSmearingTool> m_MCaSTool;

        std::vector<std::string> m_sysNames;
        std::vector<CP::SystematicSet> m_sysList;

        TFile* m_DebugFile;
        TTree* m_DebugTree;
        int m_SelCategoryRaw = 0, m_SelCategory = 0;
        std::unique_ptr<MMCTest::TrackInfo> m_Combined;
        std::unique_ptr<MMCTest::TrackInfo> m_InnerDet;
        std::unique_ptr<MMCTest::TrackInfo> m_MSExtr;
        std::unique_ptr<MMCTest::TrackInfo> m_MSOnlyExtr;

    };  // class TestMCASTTool

}  // namespace CP

#endif  // MCAST_TOOLTESTER_H
