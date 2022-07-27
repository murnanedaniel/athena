/*
    Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// TauJets_LepRMParticleThinning.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef DERIVATIONFRAMEWORK_TAUJETS_LEPRMPARTICLETHINNING_H
#define DERIVATIONFRAMEWORK_TAUJETS_LEPRMPARTICLETHINNING_H

#include <string>
#include <atomic>

#include "AthenaBaseComps/AthAlgTool.h"
#include "DerivationFrameworkInterfaces/IThinningTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTau/TauJetContainer.h"

#include "StoreGate/ThinningHandleKey.h"
#include "StoreGate/ReadHandleKey.h"

#include "ExpressionEvaluation/ExpressionParserUser.h"

namespace DerivationFramework {

    class TauJets_LepRMParticleThinning : public extends<ExpressionParserUser<AthAlgTool>, IThinningTool> {
    public: 
        TauJets_LepRMParticleThinning(const std::string& t, const std::string& n, const IInterface* p);
        virtual StatusCode initialize() override;
        virtual StatusCode finalize() override;
        virtual StatusCode doThinning() const override;

    private:
        mutable std::atomic<unsigned int> m_ntot_taus {0};
        mutable std::atomic<unsigned int> m_ntot_trks {0};
        mutable std::atomic<unsigned int> m_ntot_ID_trks {0};
        mutable std::atomic<unsigned int> m_npass_taus {0};
        mutable std::atomic<unsigned int> m_npass_trks {0};
        mutable std::atomic<unsigned int> m_npass_ID_trks {0};
        
        StringProperty m_streamName{ this, "StreamName", "", "Name of the stream being thinned" };

        SG::ReadHandleKey<xAOD::TauJetContainer> m_originalTauKey{ this, "originalTauKey", "", ""};

        SG::ThinningHandleKey<xAOD::TauJetContainer> m_LepRMTauKey{ this, "LepRMTauKey", "TauJets_LepRM", "where Lep can be Muon or Elec" };

        SG::ThinningHandleKey<xAOD::TrackParticleContainer> m_inDetSGKey{ this, "InDetTrackParticlesKey", "InDetTrackParticles", "" };

        SG::ThinningHandleKey<xAOD::TauTrackContainer> m_tauTracksSGKey{ this, "TauTracksKey", "TauTracks_LepRM", "where Lep can be Muon or Elec" }; 

        Gaudi::Property<std::string> m_selectionString{ this, "SelectionString", "",""};
    };
}

#endif // DERIVATIONFRAMEWORK_TAUJETS_LEPRMPARTICLETHINNING_H
