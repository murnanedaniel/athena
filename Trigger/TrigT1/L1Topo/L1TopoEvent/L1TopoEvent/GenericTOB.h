// Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#ifndef GENERICTOB_H
#define GENERICTOB_H

#include <iostream>
#include <string>

#include "L1TopoEvent/BaseTOB.h"
#include "L1TopoEvent/JetTOB.h"
#include "L1TopoEvent/gJetTOB.h"
#include "L1TopoEvent/jTauTOB.h"
#include "L1TopoEvent/eTauTOB.h"
#include "L1TopoEvent/jLJetTOB.h"
#include "L1TopoEvent/gLJetTOB.h"
#include "L1TopoEvent/jJetTOB.h"
#include "L1TopoEvent/ClusterTOB.h"
#include "L1TopoEvent/eEmTOB.h"
#include "L1TopoEvent/jEmTOB.h"
#include "L1TopoEvent/cTauTOB.h"
#include "L1TopoEvent/MuonTOB.h"
#include "L1TopoEvent/LateMuonTOB.h"
#include "L1TopoEvent/MuonNextBCTOB.h"
#include "L1TopoEvent/MetTOB.h"
#include "L1TopoEvent/jXETOB.h"
#include "L1TopoEvent/jTETOB.h"
#include "L1TopoEvent/gXETOB.h"
#include "L1TopoEvent/gTETOB.h"


// TODO implement sizecheck lile in ClusterTOB

namespace TCS {
   
   class GenericTOB : public BaseTOB {
   public:

      // default constructor
      GenericTOB(uint32_t roiWord = 0);

      // constructor from individual values
      GenericTOB(unsigned int Et, int eta, int phi, uint32_t roiWord = 0);

      // copy constructor
      GenericTOB(const GenericTOB & other);

      // constructor from jet
      GenericTOB(const JetTOB & jet, JetTOB::JetSize jetSize);

      // constructor from jFEX Tau
      GenericTOB(const jTauTOB & tau);

      // constructor from jFEX Em
      GenericTOB(const jEmTOB & jem);

      // constructor from jFEX LJet
      GenericTOB(const jLJetTOB & jet);

      // constructor from gFEX LJet
      GenericTOB(const gLJetTOB & jet);

      // constructor from jFEX Jet
      GenericTOB(const jJetTOB & jet);

      // constructor from gFEX Jet
      GenericTOB(const gJetTOB & jet);

      // constructor from cluster
      GenericTOB(const ClusterTOB & cluster);

      // constructor from eFEX Em
      GenericTOB(const eEmTOB & eem);

      // constructor from eFEX Tau
      GenericTOB(const eTauTOB & etau);

      // constructor from cTau
      GenericTOB(const cTauTOB & ctau);

      // constructor from muon
      GenericTOB(const MuonTOB & muon);
      
      // constructor from lateMuon
      GenericTOB(const LateMuonTOB & lateMuon);

      // constructor from muonNextBC
      GenericTOB(const MuonNextBCTOB & muonNextBC);

      // constructor from met
      GenericTOB(const MetTOB & met);

      // constructor from jFEX XE
      GenericTOB(const jXETOB & jxe);

      // constructor from jFEX TE
      GenericTOB(const jTETOB & jte);

      // constructor from gFEX XE
      GenericTOB(const gXETOB & gxe);

      // constructor from gFEX TE
      GenericTOB(const gTETOB & gte);

      // destructor
      ~GenericTOB();

      static GenericTOB* createOnHeap(const GenericTOB &);
      static void clearHeap();

      static const Heap<TCS::GenericTOB>& heap() { return fg_heap; }

   public:
      unsigned int Et() const { return m_Et; }
      unsigned int EtWide() const { return m_EtWide; }
      unsigned int EtNarrow() const { return m_EtNarrow; }

      int Ex() const { return m_Ex; }
      int Ey() const { return m_Ey; }
      unsigned int Et2() const { return m_Et2; }
      unsigned int sumEt() const { return m_sumEt; }

      int eta() const { return m_eta; }
      int phi() const { return m_phi; }

      // See definitions at TrigT1Interfaces/MuCTPIL1TopoCandidate.h 
      int bw2or3() const { return m_bw2or3; }
      int innerCoin() const { return m_innerCoin; }
      int goodMF() const { return m_goodMF; }
      int charge() const { return m_charge; }
      int is2cand() const { return m_is2cand; }
      
      double EtDouble() const { return m_EtDouble; }
      double etaDouble() const { return m_etaDouble; }
      double phiDouble() const { return m_phiDouble; }

      double ExDouble() const { return m_ExDouble; }
      double EyDouble() const { return m_EyDouble; }
      double sumEtDouble() const { return m_sumEtDouble; }

      virtual void print(std::ostream &o) const;

      void setTobType(inputTOBType_t tobType) { m_tobType = tobType; }

      inputTOBType_t tobType() const { return m_tobType; }

   private:
      unsigned int m_Et { 0 };
      unsigned int m_EtNarrow { 0 };
      unsigned int m_EtWide { 0 };

      int m_Ex { 0 };
      int m_Ey { 0 };
      unsigned int m_Et2 { 0 };
      unsigned int m_sumEt { 0 };

      int m_eta { 0 };
      int m_phi { 0 };

      int m_bw2or3 { 0 };
      int m_innerCoin { 0 };
      int m_goodMF { 0 };
      int m_charge { 0 };
      int m_is2cand { 0 };

      double m_EtDouble { 0 };
      double m_etaDouble { 0 };
      double m_phiDouble { 0 };

      double m_ExDouble { 0 };
      double m_EyDouble { 0 };
      double m_sumEtDouble { 0 };

      inputTOBType_t   m_tobType { NONE };

      static thread_local  Heap<TCS::GenericTOB> fg_heap;
   };  
}

#endif
