/* emacs: this is -*- c++ -*- */
/**
 **     @file    Associator_TruthMatch.h
 **
 **     @author  mark sutton
 **     @date    Fri 11 Jan 2019 07:06:39 CET 
 **
 **     Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 **/


#ifndef TIDAUTILS_ASSOCIATOR_TRUTHMATCH_H
#define TIDAUTILS_ASSOCIATOR_TRUTHMATCH_H

#include <set>
#include <map>


#include "TrigInDetAnalysis/TrackAssociator.h"
#include "TrigInDetAnalysisUtils/Associator_BestMatch.h"
#include "TrigInDetAnalysis/Track.h"



#if 0
class Associator_TruthMatcher : public Associator_BestMatcher {

public:
  Associator_TruthMatcher() :
    Associator_MatcherBase( "Truth" )
  { }

  virtual double distance( TIDA::Track* refTrack, TIDA::Track* testTrack ) {
    if      (testTrack->match_barcode() == -1)                  return 1;
    else if (testTrack->match_barcode() == refTrack->barcode()) return 0;
    else                                                        return 1;
  }
  
};
#endif


// class Associator_TruthMatcher : public Associator_MatcherBase{
class Associator_TruthMatcher : public TrackAssociator {
  
public:

  Associator_TruthMatcher() : TrackAssociator("Truth")  { }
  
  virtual ~Associator_TruthMatcher() { }

  virtual TrackAssociator* clone() override { return new Associator_TruthMatcher(*this); }

  //Fill reference tracks in matching step
  virtual void match( const std::vector<TIDA::Track*>& refTracks, 
		      const std::vector<TIDA::Track*>& testTracks) override {

    //std::cout<<"refTracks.size() "<<refTracks.size()<<" \t testTracks.size() "<<testTracks.size()<<std::endl;
    
    for (unsigned int i = 0; i < refTracks.size(); i++) {
      
      //        std::cout<<refTracks[i]->author() <<std::endl;
      
      for (unsigned int j = 0; j < testTracks.size(); j++) {

	
	if ( distance( refTracks[i], testTracks[j] ) < 1. ) {
	  //              std::cout<<"MATCHED"<<std::endl;
	  mmatched.insert(    map_type::value_type(refTracks[i],testTracks[j]));
	  mrevmatched.insert( map_type::value_type(testTracks[j],refTracks[i]));
	}
      }
    }
    
    return;
  }
  
  
  virtual double distance( TIDA::Track* refTrack, TIDA::Track* testTrack ) {
      if      (testTrack->match_barcode() == -1)                   return 1;
      else if (testTrack->match_barcode() == refTrack->barcode())  return 0;
      else                                                         return 1;
  }
  
private:

  // double md;

};


#endif //  TIDAUTILS_ASSOCIATOR_TRUTHMATCH_H
                                                        
