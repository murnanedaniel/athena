//  LateMuonTOB.h 
//  TopoCore
//  Based on MuonTOB.h Created by Joerg Stelzer on 11/10/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.

#ifndef __TopoCore__LateMuonTOB__
#define __TopoCore__LateMuonTOB__

#include <iostream>

#include "L1TopoEvent/BaseTOB.h"
#include "L1TopoEvent/Heap.h"

namespace TCS {
    
   class LateMuonTOB : public BaseTOB {
   public:
      
      static unsigned int nBitsEt() { return g_nBitsEt; }
      static unsigned int nBitsIsolation() { return g_nBitsIsolation; }
      static unsigned int nBitsEta() { return g_nBitsEta; }
      static unsigned int nBitsPhi() { return g_nBitsPhi; }


      // default constructor
      LateMuonTOB(uint32_t roiWord = 0);
      
      // constructor with individual values
      LateMuonTOB(unsigned int et, unsigned int isolation, int eta, int phi, uint32_t roiWord = 0);

      // constructor with initial values
      LateMuonTOB(const LateMuonTOB & latemuon);

      // destructor
      virtual ~LateMuonTOB();

      // accessors
      unsigned int Et() const { return m_Et; }
      unsigned int isolation() const { return m_isolation; }
      int eta() const { return m_eta; }
      int phi() const { return m_phi; }

      double EtDouble() const { return m_EtDouble; }
      double EtaDouble() const { return m_etaDouble; }
      double PhiDouble() const { return m_phiDouble; }
      
      // setters
      void setEt(unsigned int et) { m_Et = sizeCheck(et, nBitsEt()); }
      void setIsolation(unsigned int et) { m_isolation = sizeCheck(et, nBitsIsolation()); }
      void setEta(int eta) { m_eta = sizeCheck(eta, nBitsEta()); }
      void setPhi(int phi) { m_phi = sizeCheck(phi, nBitsPhi()); }
      
      void setEtDouble(double et) { m_EtDouble = et; }
      void setEtaDouble(double eta) { m_etaDouble = eta; }
      void setPhiDouble(double phi) { m_phiDouble = phi; }

      inputTOBType_t tobType() const { return LATEMUON; }
      
      // memory management
      static LateMuonTOB* createOnHeap(const LateMuonTOB& cl);
      static void clearHeap();
      static const Heap<TCS::LateMuonTOB>& heap() { return fg_heap; }

      static unsigned int instances() { return fg_instances; }

      virtual void print(std::ostream &o) const;

   private:
      static unsigned int g_nBitsEt;
      static unsigned int g_nBitsIsolation;
      static unsigned int g_nBitsEta;
      static unsigned int g_nBitsPhi;
      
      unsigned int m_Et;
      unsigned int m_isolation;
      int m_eta;
      int m_phi;

      double m_EtDouble;
      double m_etaDouble;
      double m_phiDouble;

      static unsigned int fg_instances;

      static Heap<TCS::LateMuonTOB> fg_heap;
   };
   
} // end of namespace TCS

#endif /* defined(__TopoCore__LateMuonTOB__) */
