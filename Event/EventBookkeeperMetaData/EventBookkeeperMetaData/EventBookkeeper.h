/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////// 
// Header file for class EventBookkeeper
// Author: David Cote, September 2008. <david.cote@cern.ch>
/////////////////////////////////////////////////////////////////// 
#ifndef EVENTBOOKKEEPER_H 
#define EVENTBOOKKEEPER_H 
 
#include "AthenaKernel/CLASS_DEF.h"
#include <string>
#include <iosfwd>
#include <stdint.h>
#include <vector>

#include "GaudiKernel/MsgStream.h"

//fwd declare
class EventBookkeeperCollection;

class EventBookkeeper
{
 public:
  //Constructors
  EventBookkeeper();
  EventBookkeeper(const std::string &name);
  EventBookkeeper(const std::string &name, 
                  const std::string &description, 
                  const std::string &logic);
  // Copy constructors: 
  EventBookkeeper( const EventBookkeeper& rhs );
  EventBookkeeper& operator=(const EventBookkeeper& rhs);

  //  Destructor
  ~EventBookkeeper();

  void Print(const std::string &indent="", std::ostream& out = std::cout);
  void PrintToMsg(MsgStream &, const std::string &indent="");
  void PrintFamily(const std::string &indent="");

  // get() and set() methods
  const std::string& getName() const { return m_name; }
  void setName( const std::string& name );

  const std::string& getInputStream() const { return m_inputstream; }
  void setInputStream( const std::string& inputstream );

  const std::string& getOutputStream() const { return m_outputstream; }
  void setOutputStream( const std::string& outputstream );
  void setOutputStreamOfFamily( const std::string &outputstream );

  const std::string& getDescription() const { return m_description; }
  void setDescription( const std::string &description );

  const std::string& getLogic() const { return m_logic; }
  void setLogic( const std::string& logic );

  uint64_t getNAcceptedEvents() const { return m_nAcceptedEvents; }
  void setNAcceptedEvents( uint64_t nEvents );
  void addNAcceptedEvents( uint64_t nEvents );

  double getNWeightedAcceptedEvents() const { return m_nWeightedAcceptedEvents; }
  void setNWeightedAcceptedEvents( double nWeightedEvents );
  void addNWeightedAcceptedEvents( double nWeightedEvents );

  void updateAcceptedEventsIncludingFamily(const EventBookkeeper* eb);

  int getCycle() const { return m_cycle; }
  void setCycle( int cycle );

  const std::vector<EventBookkeeper*>* getChildrenEventBookkeepers() const { return m_childrenEB; }
  void fillWithWholeFamily( EventBookkeeperCollection* family );
  void setChildrenEventBookkeepers(std::vector<EventBookkeeper*>*  childrenEB );
  void AddChild(EventBookkeeper* eb);
  void AddChildren( std::vector<EventBookkeeper*>* children );
  EventBookkeeper* AddNewChild(const std::string& name,
                               const std::string& description);

  bool isEqualTo( const EventBookkeeper *eb );

 private:
  void SetDefaultDataMemberValues();

  std::string m_name;
  std::string m_description;
  std::string m_inputstream;
  std::string m_outputstream;
  std::string m_logic;
  uint64_t m_nAcceptedEvents;
  double m_nWeightedAcceptedEvents;
  int m_cycle;
  std::vector<EventBookkeeper*>* m_childrenEB;

  //Additional functions and data for EventBookkeeperCollection with flat structure 
  //This special mode is only foreseen for dumping in a TTree, please don't use it otherwise
  //In this case, the m_childrenEB are not usable, being replaced m_childrenIndices
  friend class CutFlowSvc;
  friend class EventBookkeeperCollection;
  EventBookkeeper* DeepCopyForFlatStructure( EventBookkeeperCollection* collFLAT );
  int m_parentIndex;
  std::vector<unsigned int>* m_childrenIndices; 

  //Helper data members for CutFlowSvc (transient-only)
  bool m_declaredChildFilter;
  bool m_declaredTopFilter;
};


//this is automatically generated by: 'clid -m EventBookkeeper'
CLASS_DEF( EventBookkeeper , 241669778 , 1 )

#endif //> EVENTBOOKKEEPER_H
