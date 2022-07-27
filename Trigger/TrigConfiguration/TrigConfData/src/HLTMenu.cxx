/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigConfData/HLTMenu.h"

using TV = boost::property_tree::ptree::value_type;  // tree-value type
using namespace std;

TrigConf::HLTMenu::HLTMenu()
{}

TrigConf::HLTMenu::HLTMenu(const boost::property_tree::ptree & data) 
   : DataStructure(data)
{
   load();
}

void
TrigConf::HLTMenu::load()
{
   if(! isInitialized() || empty() ) {
      return;
   }
   m_name = getAttribute("name");
}

void
TrigConf::HLTMenu::clear()
{
   m_smk = 0;
}

std::size_t 
TrigConf::HLTMenu::size() const
{
   return data().get_child("chains").size();
}

unsigned int
TrigConf::HLTMenu::smk() const {
   return m_smk;
}

void
TrigConf::HLTMenu::setSMK(unsigned int smk) {
   m_smk = smk;
}

TrigConf::HLTMenu::const_iterator
TrigConf::HLTMenu::begin() const
{
   return {data().get_child("chains"), 0,  [](auto & x){auto chain = Chain(x.first, x.second); return chain; }};
}

TrigConf::HLTMenu::const_iterator
TrigConf::HLTMenu::end() const
{
   auto & pt = data().get_child("chains");
   return { pt, pt.size(), [](auto & x){return Chain(x.second);} };
}


std::vector<TrigConf::DataStructure>
TrigConf::HLTMenu::streams() const
{
   std::vector<DataStructure> strlist;
   auto streams = data().get_child_optional("streams");
   if(streams) {
      strlist.reserve(streams->size());
      for( auto & strData : *streams ) {
         strlist.emplace_back( strData.second );
      }
   }
   return strlist;
}


std::map<std::string, std::vector<std::string>>
TrigConf::HLTMenu::sequencers() const
{
   std::map<std::string, std::vector<std::string>> result;
   const auto & sequencers = getObject("sequencers");

   for( auto & sequence : sequencers.getKeys() ) {
      for( auto & alg : sequencers.getList(sequence) ) {
         result[sequence].emplace_back(alg.getValue<std::string>());
      }
   }

   return result;
}


void
TrigConf::HLTMenu::printMenu(bool full) const
{
   cout << "HLT menu '" << name() << "'" << endl;
   cout << "Streams: " << streams().size() << endl;
   cout << "Chains: " << size() << endl;
   if(full) {
      int c(0);
      for( auto & chain : *this ) {
         cout << "  " << c++ << ": " << chain.name() << endl;
      }
   }
}
