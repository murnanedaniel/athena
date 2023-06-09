// emacs: this is -*- c++ -*-
/** @file Associator_DeltaRCosmics.h */

#ifndef TrigInDetAnalysisUtils_Associator_DeltaRCosmics_H
#define TrigInDetAnalysisUtils_Associator_DeltaRCosmics_H


#include <iostream>
#include <string>
#include <cmath>
#include <map>

#include "TrigInDetAnalysis/TrackAssociator.h"
#include "TrigInDetAnalysis/Track.h"


class Associator_DeltaRCosmics : public TrackAssociator {

public:
  
  Associator_DeltaRCosmics(const std::string& name, double deltaR) : TrackAssociator(name), m_deltaR2(deltaR*deltaR) {} 
  
  ~Associator_DeltaRCosmics() { } 
  
  virtual void match(const std::vector<TIDA::Track*>& referenceTracks, 
		     const std::vector<TIDA::Track*>& testTracks) {

    // Clear previously filled association map
    clear();
    
    // Loop over reference tracks
    std::vector<TIDA::Track*>::const_iterator reference, referenceEnd=referenceTracks.end();
    for(reference=referenceTracks.begin(); reference!=referenceEnd; reference++) {
      
      // Loop over test tracks and find the closest
      TIDA::Track* bestMatch = NULL;
      double bestDeltaR=1000;
      
      std::vector<TIDA::Track*>::const_iterator test, testEnd=testTracks.end();
      for(test=testTracks.begin(); test!=testEnd; test++) {

	// Evaluate distance between reference and test tracks
	double deta = (*reference)->eta() + (*test)->eta();
	double dphi = (*reference)->phi() - (*test)->phi() - M_PI;
	if(dphi>M_PI)  dphi-=2*M_PI;
	if(dphi<-M_PI) dphi+=2*M_PI;
	double deltaR =  deta*deta+dphi*dphi;

	// Check if this is the best match so far
	if(bestMatch==NULL || deltaR<bestDeltaR) { 
	  bestDeltaR = deltaR;
	  bestMatch  = (*test);
	} 
      }

      // Check if the best match is within delta R specifications
      if(bestMatch && bestDeltaR<m_deltaR2) { 
	// Create reference->test and test->reference associations
	mmatched.insert(map_type::value_type(*reference, bestMatch));
	mrevmatched.insert( map_type::value_type(bestMatch, *reference));
      }
    }
  }
  
private:
  
  double m_deltaR2;
  
};


#endif  // TrigInDetAnalysisUtils_Associator_DeltaRCosmics_H
