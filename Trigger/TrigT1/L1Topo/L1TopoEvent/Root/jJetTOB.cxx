//  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#include "L1TopoEvent/jJetTOB.h"

thread_local TCS::Heap<TCS::jJetTOB> TCS::jJetTOB::fg_heap("Jet");

// constructors
// default constructor
TCS::jJetTOB::jJetTOB(uint32_t roiWord, const std::string& tobName) :
   BaseTOB( roiWord,tobName )
{}

// constructor with initial values
TCS::jJetTOB::jJetTOB(unsigned int Et, int eta, unsigned phi, uint32_t roiWord, const std::string& tobName) :
   BaseTOB( roiWord,tobName )
   , m_Et(Et)
   , m_eta(eta)
   , m_phi(phi)
{}

// copy constructor
TCS::jJetTOB::jJetTOB(const jJetTOB & jet) = default;

TCS::jJetTOB::~jJetTOB() = default;

TCS::jJetTOB*
TCS::jJetTOB::createOnHeap(const jJetTOB& jet) {
   return fg_heap.create(jet);
}

void
TCS::jJetTOB::clearHeap() {
   return fg_heap.clear();
}

void
TCS::jJetTOB::print(std::ostream &o) const {
   o << "jJet energy: " << Et() << ", eta: " << eta() << ", phi: " << phi();
}
