#ifndef ISOLATIONSELECTION_TESTMARCOHELPERS_H
#define ISOLATIONSELECTION_TESTMARCOHELPERS_H

#include <AsgTools/ToolHandle.h>
#include <IsolationSelection/IsolationCloseByCorrectionTool.h>
#include <IsolationSelection/IsolationSelectionTool.h>

#include <xAODBase/IParticle.h>

#include <xAODBase/IParticleContainer.h>
#include <TTree.h>
#include <vector>
namespace CP{
        class IsolationCloseByCorrectionTool;
        class IsoCorrectionTestHelper {
            public:  
                IsoCorrectionTestHelper(TTree* outTree, const std::string& ContainerName);
                StatusCode Fill( xAOD::IParticleContainer* Particles);
                
               
                
            private:
                
                float Charge(const xAOD::IParticle* P) const;
                
                std::vector<float> m_pt;
                std::vector<float> m_eta;
                std::vector<float> m_phi;
                std::vector<float> m_e;
                std::vector<int> m_Q;
                
                std::vector<float> m_orig_TrackIsol;
                std::vector<float> m_corr_TrackIsol;
                
                std::vector<float> m_orig_CaloIsol;
                std::vector<float> m_corr_CaloIsol;
                
                std::vector<char> m_orig_passIso;
                std::vector<char> m_corr_passIso;
                
                
        
        };
        
    
}

#endif
