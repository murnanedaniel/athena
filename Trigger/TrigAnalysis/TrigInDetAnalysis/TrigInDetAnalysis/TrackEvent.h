// emacs: this is -*- c++ -*-
//
//   @file    TrackEvent.h        
//
//            Basic event class to contain a vector of
//            chains for trigger analysis                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: TrackEvent.h, v0.0   Mon  1 Feb 2010 11:43:51 GMT sutt $


#ifndef __TRACKEVENT_H
#define __TRACKEVENT_H

#include <iostream>
#include <vector>
#include <string>

#include "TrigInDetAnalysis/TrackChain.h"
#include "TrigInDetAnalysis/TrackVertex.h"
//#include "TrigInDetTruthEvent/TrigInDetTrackTruthMap.h"


#include "TObject.h"



class TrackEvent : public TObject {

public:

  TrackEvent();

  virtual ~TrackEvent();

  /// accessors
  void run_number(unsigned  r)    { m_run_number     = r; }  
  void event_number(unsigned e)   { m_event_number   = e; }  
  void lumi_block(unsigned lb)    { m_lumi_block     = lb; }  
  void time_stamp(unsigned t)     { m_time_stamp     = t; }  
  void bunch_crossing_id(unsigned b) { m_bunch_crossing_id = b; }  
	void mu(double m) { m_mu = m;}

  unsigned run_number()     const { return m_run_number;   } 
  unsigned event_number()   const { return m_event_number; }  
  unsigned lumi_block()     const { return m_lumi_block;   }  
  unsigned time_stamp()     const { return m_time_stamp;   }  
  unsigned bunch_crossing_id() const { return m_bunch_crossing_id;  }  

  /// FIXME: what is this ? need a comment describing any not 
  ///        descriptive variable name
  double mu() const { return m_mu; } /// vertex multiplicity ?   

  /// NB all these could be avoided simply be inheriting 
  ///    from and std::vector<TrackChain> rather than 
  ///    a member variable 

  /// number of chains added to this event
  unsigned size() const { return m_chains.size(); } 

  /// methods to add and access chains 
  void addChain(const std::string& chainname) { 
     m_chains.push_back(TrackChain(chainname));
  }

  void addVertex(const TrackVertex& v) { 
     m_vertices.push_back(v);
  }


  const std::vector<TrackChain>& chains() const { return m_chains; };
  std::vector<TrackChain>&       chains()       { return m_chains; };
  
  //void setTruthMap(TrigInDetTrackTruthMap truthmap) {
  //	m_truthmap = truthmap;
  //}
  
  /// clear the event
  void clear() { m_chains.clear(); m_vertices.clear(); } 
 
  /// get the last chain from the vector 
  TrackChain& back() { return m_chains.back(); }
  
  // iterators
  std::vector<TrackChain>::iterator begin() { return m_chains.begin(); }
  std::vector<TrackChain>::iterator end()   { return m_chains.end(); }

  std::vector<TrackChain>::const_iterator begin() const { return m_chains.begin(); }
  std::vector<TrackChain>::const_iterator end()   const { return m_chains.end(); }
  
  /// vector operator
  TrackChain& operator[](int i) { return m_chains.at(i); }

  const std::vector<TrackVertex> vertices() const { return m_vertices; }

  std::vector<std::string> chainnames() const;

  void erase( const std::string& name );


private:
 
  unsigned m_run_number; 
  unsigned m_event_number;  
  unsigned m_lumi_block;
  unsigned m_time_stamp;

  unsigned m_bunch_crossing_id;
  double   m_mu;  /// vertex multiplicity ?   

  /// trigger chain information
  std::vector<TrackChain> m_chains;

  std::vector<TrackVertex> m_vertices;
  
  ClassDef(TrackEvent,3)

};


inline std::ostream& operator<<( std::ostream& s, const TrackEvent& t ) { 
  s << "Event run: " << t.run_number() 
    << "\tevent: "   << t.event_number() 
    << "\tlb: "      << t.lumi_block() 
    << "\tbc: "      << t.bunch_crossing_id()
    << "\ttime: "    << t.time_stamp();
  for ( unsigned i=0 ; i<t.vertices().size() ; i++ ) s << "\n" << t.vertices()[i];
  for ( unsigned i=0 ; i<t.chains().size()   ; i++ ) s << "\n" << t.chains()[i];
  return s;
}


#endif  // __TRACKEVENT_H 










