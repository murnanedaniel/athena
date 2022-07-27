/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "PersistentDataModelTPCnv/DataHeader_p5.h"

#include "CxxUtils/MD5.h"

#include <uuid/uuid.h>
#include <sstream>

DataHeaderElement_p5::DataHeaderElement_p5() : m_token(), m_oid2(0U) {}
DataHeaderElement_p5::DataHeaderElement_p5(const DataHeaderElement_p5& rhs) : m_token(rhs.m_token),
	m_oid2(rhs.m_oid2) {}
DataHeaderElement_p5::~DataHeaderElement_p5() {}

DataHeaderElement_p5& DataHeaderElement_p5::operator=(const DataHeaderElement_p5& rhs) {
   if (this != &rhs) {
      m_token = rhs.m_token;
      m_oid2 = rhs.m_oid2;
   }
   return(*this);
}

const std::string& DataHeaderElement_p5::token() const {
   return(m_token);
}

long long int DataHeaderElement_p5::oid2() const {
   return(m_oid2);
}


DataHeaderForm_p5::DataHeaderForm_p5() : m_map(), m_uints() {}
DataHeaderForm_p5::DataHeaderForm_p5(const DataHeaderForm_p5& rhs) : m_map(rhs.m_map), m_uints(rhs.m_uints) {}
DataHeaderForm_p5::~DataHeaderForm_p5() {}
DataHeaderForm_p5& DataHeaderForm_p5::operator=(const DataHeaderForm_p5& rhs) {
   if (&rhs != this) {
      m_map = rhs.m_map;
      m_uints = rhs.m_uints;
   }
   return(*this);
}

const std::vector<std::string>& DataHeaderForm_p5::map() const {
   return(m_map);
}

void DataHeaderForm_p5::insertMap(const std::string& element) {
   m_map.push_back(element);
}

const std::vector<unsigned int>& DataHeaderForm_p5::params(unsigned int entry) const {
   return(m_uints[entry - 1]);
}

void DataHeaderForm_p5::insertParam(unsigned int param, unsigned int entry) {
   m_uints[entry - 1].push_back(param);
}

unsigned int DataHeaderForm_p5::size() const {
   return(m_uints.size());
}

void DataHeaderForm_p5::resize(unsigned int size) {
   m_uints.resize(size);
}


DataHeader_p5::DataHeader_p5() : m_dataHeader(), m_dhFormToken(), m_dhFormMdx() {}
DataHeader_p5::DataHeader_p5(const DataHeader_p5& rhs) : m_dataHeader(rhs.m_dataHeader),
	m_dhFormToken(rhs.m_dhFormToken),
	m_dhFormMdx(rhs.m_dhFormMdx) {}
DataHeader_p5::~DataHeader_p5() {
}

DataHeader_p5& DataHeader_p5::operator=(const DataHeader_p5& rhs) {
   if (this != &rhs) {
      m_dataHeader = rhs.m_dataHeader;
      m_dhFormToken = rhs.m_dhFormToken;
      m_dhFormMdx = rhs.m_dhFormMdx;
   }
   return(*this);
}

const std::vector<DataHeaderElement_p5>& DataHeader_p5::elements() const {
   return(m_dataHeader);
}

const std::string& DataHeader_p5::dhFormToken() const {
   return(m_dhFormToken);
}

void DataHeader_p5::setDhFormToken(const std::string& formToken,
                                   const DataHeaderForm_p5& dhForm)
{
  m_dhFormToken = formToken;
  std::ostringstream stream;
  for (const std::string& s : dhForm.map()) {
    stream << s << "\n";
  }
  for (unsigned int entry = 1; entry <= dhForm.size(); ++entry) {
    for (unsigned int x : dhForm.params(entry)) {
      stream << x << ",";
    }
    stream << "\n";
  }
  MD5 checkSum((unsigned char*)stream.str().c_str(), stream.str().size());
  uuid_t checkSumUuid;
  checkSum.raw_digest((unsigned char*)(&checkSumUuid));
  char text[37];
  uuid_unparse_upper(checkSumUuid, text);
  m_dhFormMdx = text;
}
const std::string& DataHeader_p5::dhFormMdx() const {
   return(m_dhFormMdx);
}
