/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef jetsubstructuremomenttools_jetsubstructurebase_header
#define jetsubstructuremomenttools_jetsubstructurebase_header

#include "xAODCaloEvent/CaloCluster.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODBase/IParticle.h"

#include "JetRec/JetModifierBase.h"
#include <vector>

namespace fastjet {
  class PseudoJet;
}

class JetSubStructureMomentToolsBase :
  public JetModifierBase {
    public:
      // Constructor and destructor
      JetSubStructureMomentToolsBase(std::string name);

    protected:
      bool checkForConstituents(const xAOD::Jet &jet) const {
        if(jet.numConstituents() == 0) {
          ATH_MSG_WARNING("Attempting to use a substructure tool on a jet that has no constituent");
          return false;
        } else {
          return true;
        }
      }
      fastjet::PseudoJet buildPseudoJet (const xAOD::Jet & jet) const;
      fastjet::PseudoJet buildPseudoJet(const std::vector<const xAOD::IParticle*>& iparticles) const;
  };


#endif
