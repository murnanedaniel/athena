/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigConfL1Data/DiffStruct.h"
#include "TrigConfL1Data/HelperFunctions.h"
#include "TrigConfL1Data/L1DataDef.h"
#include "TrigConfL1Data/TriggerItem.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

#include "boost/lexical_cast.hpp"

using namespace std;
using namespace TrigConf;

TriggerItem::TriggerItem() : TrigConfData()
{}

TriggerItem::~TriggerItem() {
   delete m_TopNode;
}


uint16_t
TrigConf::TriggerItem::bunchgroupMask() const {
   uint16_t bgmask(0);
   m_TopNode->getBunchGroupsMask(bgmask);
   return bgmask;
}


TrigConf::DiffStruct*
TrigConf::TriggerItem::compareTo(const TriggerItem* o) const {
   DiffStruct * ds = new DiffStruct("TriggerItem");
   ds->check("ctpid", ctpId(), o->ctpId());
   ds->check("definition", definition(), o->definition());
   ds->check("complex_deadtime", complex_deadtime(), o->complex_deadtime());
   ds->check("triggertype", triggerType(), o->triggerType());
   if(!ds->empty()) {
      ds->name = name();
   } else {
      delete ds; ds=0;
   }
   return ds;
}


void
TriggerItem::writeXML(std::ostream & xmlfile, int indentLevel, int indentWidth) const {

   string final_def = m_Definition;
   bool end=false;
   string::size_type index=0;
   while(!end){
      std::string::size_type pos=final_def.find("&",index);
      if(pos!= std::string::npos){
         end=false;
         index = pos+1;
         //do the replace
         final_def.replace(pos, 1, "&amp;");
      }else{
         end=true;
         index=0;
      }
   }

   indent(xmlfile, indentLevel, indentWidth)
      << "<TriggerItem ctpid=\"" << m_CtpId
      << "\" partition=\"" << partition()
      << "\" name=\"" << name() 
      << "\" complex_deadtime=\"" << complex_deadtime()
      << "\" definition=\"" << final_def
      << "\" trigger_type=\"" << TrigConf::uint2bin(m_TriggerType, partition()==1 ? 8 : 4)
      << "\" version=\"" << version() << "\">" << endl;

   if (m_TopNode)
      m_TopNode->writeXML(xmlfile, indentLevel+1, indentWidth);

   indent(xmlfile, indentLevel, indentWidth) << "</TriggerItem>" << endl;
}


void
TriggerItem::buildLogic(std::string & logic, std::vector<std::string> & conditionList) const {
   if (m_TopNode)
      m_TopNode->buildLogic(conditionList, logic);
}


void
TriggerItem::setCondition(const string& logic,
                          const vector<string>& conditions,
                          const vector<TriggerThreshold*>& thrs) {
   // logic is a string like "!(1|2)&(3|4)"
   // vector of conditions of length equal the largest number in logic string (e.g 4)
   // a condition can be either a single name for internal triggers or a komma separated list of multiplicity, conditions name and threshold name
   // thrs is a list of thresholds
  
   m_Definition = "("+logic+")";
   
   try {
      m_TopNode = TrigConf::parse(logic, conditions, thrs);
   }
   catch(exception e) {
      cout << "WARNING: Could not set condition of triggeritem " << name() << " because: " << e.what() << endl;
   }
}


void
TriggerItem::print(const std::string& indent, unsigned int detail) const {
   if(detail>=3) {
      cout << indent << "==================================" << endl; 
      cout << indent << "TriggerItem:  "; 
      printNameIdV(indent);
      cout << indent << "Complex Deadtime: " << m_ComplexDeadtime << endl;
      cout << indent << "Definition:       " << m_Definition << endl;
      cout << indent << "CTP_id  :         " << m_CtpId << endl;
      cout << indent << "TriggerType(hex): 0x" << std::hex << m_TriggerType << std::dec << endl;
      cout << indent << "Bunchgroup mask : 0x" << uint2bin(bunchgroupMask(), 16) << endl;
      if (m_TopNode) m_TopNode->print(indent);
      else cout << indent << "TopNode is NULL pointer" << endl;
      cout << indent << "==================================" << endl; 
   }
}


std::ostream &
TrigConf::operator<<(std::ostream & o, const TrigConf::TriggerItem & item) {
   o << "TriggerItem: " << (const TrigConfData&)item;
   o << "   Complex Deadtime: " << item.m_ComplexDeadtime    << endl;
   o << "   Definition:       " << item.m_Definition  << endl;
   o << "   CTP_id:           " << item.m_CtpId       << endl;
   o << "   TriggerType(hex): " << hex << item.m_TriggerType << dec << endl;
   if (item.topNode()) o << *item.topNode();
   else o << "   No TopNode" << endl;
   return o;
}


string
TriggerItem::__str__() const {
   stringstream s;
   s << *this;
   return s.str();
}
