/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/
//  eEmSort.cxx
//  TopoCore
//  algorithm to create sorted lists for eEMs, et order applied
//
#include "L1TopoAlgorithms/eEmSort.h"
#include "L1TopoEvent/TOBArray.h"
#include "L1TopoEvent/ClusterTOBArray.h"
#include "L1TopoEvent/GenericTOB.h"
#include <algorithm>

REGISTER_ALG_TCS(eEmSort)

bool SortByEtLargesteEm(TCS::GenericTOB* tob1, TCS::GenericTOB* tob2)
{
   return tob1->Et() > tob2->Et();
}

// constructor
TCS::eEmSort::eEmSort(const std::string & name) : SortingAlg(name) {
   defineParameter( "InputWidth", 120 ); // for FW
   defineParameter( "InputWidth1stStage", 30 ); // for FW
   defineParameter( "OutputWidth", 6);
   defineParameter( "IsoMask", 0);
   defineParameter( "MinEta", 0);
   defineParameter( "MaxEta", 63);
   defineParameter( "DoIsoCut", 1);
}


// destructor
TCS::eEmSort::~eEmSort() {}


TCS::StatusCode
TCS::eEmSort::initialize() {
   m_numberOfClusters = parameter("OutputWidth").value();
   m_iso = parameter("IsoMask").value();
   m_minEta = parameter("MinEta").value();
   m_maxEta = parameter("MaxEta").value();
   m_doIsoCut = parameter( "DoIsoCut").value();
   return TCS::StatusCode::SUCCESS;
}


TCS::StatusCode
TCS::eEmSort::sort(const InputTOBArray & input, TOBArray & output) {

   const ClusterTOBArray & clusters = dynamic_cast<const ClusterTOBArray&>(input);

   // fill output array with GenericTOB buildt from clusters
   for(ClusterTOBArray::const_iterator cl = clusters.begin(); cl!= clusters.end(); ++cl ) {
      const GenericTOB gtob(**cl);

      if (parType_t(std::abs((*cl)-> eta())) < m_minEta) continue; 
      if (parType_t(std::abs((*cl)-> eta())) > m_maxEta) continue;
      // isolation cut
      if (m_iso != 0 ) {
          unsigned int isobit(0x1 << (m_iso-1));
          if(m_doIsoCut && ((parType_t((*cl)->isolation()) & isobit) != isobit)) continue;
      }
      
      output.push_back( gtob );
   }

   // sort
   output.sort(SortByEtLargesteEm);
   

   // keep only max number of clusters
   int par = m_numberOfClusters;
   unsigned int maxNumberOfClusters = (unsigned int)(par<0?0:par);
   if(maxNumberOfClusters>0) {
      while( output.size()> maxNumberOfClusters ) {
         output.pop_back();
      }
   }
   return TCS::StatusCode::SUCCESS;
}

