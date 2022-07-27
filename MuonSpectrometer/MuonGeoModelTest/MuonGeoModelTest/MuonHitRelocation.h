/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUON_HIT_RELOCATION_H
#define MUON_HIT_RELOCATION_H

#include <string>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonSimEvent/CscHitIdHelper.h"
#include "MuonSimEvent/MdtHitIdHelper.h"
#include "MuonSimEvent/MicromegasHitIdHelper.h"
#include "MuonSimEvent/RpcHitIdHelper.h"
#include "MuonSimEvent/TgcHitIdHelper.h"
#include "MuonSimEvent/sTgcHitIdHelper.h"

// Histograms quantities related to G4 Hits in the Atlas
// Detector.  Used to validate GeoModel + much more.

class MuonHitRelocation : public AthAlgorithm {
public:
    MuonHitRelocation(const std::string& name, ISvcLocator* pSvcLocator);
    ~MuonHitRelocation();
    StatusCode initialize();
    StatusCode execute();
    long getIntStName(const std::string& stName) const;

private:
    bool m_checkMdt;
    bool m_checkRpc;
    bool m_checkTgc;
    bool m_checkCsc;
    bool m_checksTgc;
    bool m_checkMM;

    class Clockwork;
    Clockwork* m_c;
    const RpcHitIdHelper* m_rmuonHelper;
    const MdtHitIdHelper* m_mmuonHelper;
    const TgcHitIdHelper* m_tmuonHelper;
    const CscHitIdHelper* m_cmuonHelper;
    const sTgcHitIdHelper* m_stmuonHelper;
    const MicromegasHitIdHelper* m_mmmuonHelper;
    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc{this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};
};

#endif
