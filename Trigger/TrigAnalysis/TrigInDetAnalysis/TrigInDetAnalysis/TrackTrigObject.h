/* emacs: this is -*- c++ -*- */
/**
 **     @file    TrackTrigObject.h
 **
 **     @author  mark sutton
 **     @date    Mon  5 Nov 2012 00:07:15 GMT 
 **
 **     Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 **/


#ifndef  TRACKTRIGOBJECT_H
#define  TRACKTRIGOBJECT_H

#include <iostream>

#include "TrigInDetAnalysis/Track.h"

#include "TObject.h"

class TrackTrigObject : public TObject {

public:

  TrackTrigObject() { } /// just for root, don't want this really
  TrackTrigObject(double eta, double phi, double pt, double z0, 
		  int type, unsigned long id=0);

  TrackTrigObject(const TIDA::Track& t, 
		  int type, unsigned long id=0);

  virtual ~TrackTrigObject() { } 
  
  double eta() const { return m_eta; }
  double phi() const { return m_phi; }
  double pt()  const { return m_pt; }
  double z0()  const { return m_z0; }
  
  int  type()  const { return m_type; }

  unsigned long id() const { return m_id; }
  
  void addChild(unsigned long id) { m_children.push_back(id); }
  
  const std::vector<unsigned long>&  children() const { return m_children; }
  std::vector<unsigned long>         children()       { return m_children; }
  
private:
  
  /// parameters
  double m_eta;
  double m_phi;
  double m_pt;
  double m_z0;

  /// generic type information
  int    m_type;

  /// identifier
  unsigned long m_id;

  /// identifiers of children
  std::vector<unsigned long> m_children; 

  ClassDef(TrackTrigObject,1);

};


inline std::ostream& operator<<( std::ostream& s, const TrackTrigObject& t ) { 
  s << "Object [" 
    << "\ttype=" << t.type()
    << ":\teta=" << t.eta()
    << "\tphi="  << t.phi()
    << "\tpt="   << t.pt()*0.001 << " GeV"
    << " ] \t(";
  
  std::string n = "";
  if ( t.children().size()>3 ) n = "\n\t"; 
  for ( size_t i=0 ; i<t.children().size() ; i++ ) s << n << "  0x" << std::hex << t.children()[i] << std::dec;
  s << "  )" ;

  return s;
}


#endif  // TRACKTRIGOBJECT_H 










