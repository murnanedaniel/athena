/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// Automatically generated from scripts in /home/hessey/prog/gmx2geo/makeshape
//    Handle tube shape elements.
//
#ifndef GEO_MODEL_XML_MAKE_CONS_H
#define GEO_MODEL_XML_MAKE_CONS_H

#include "GeoModelXml/Element2GeoItem.h"

class MakeCons: public Element2GeoItem {
public:
    MakeCons();
    const RCBase * make(const xercesc::DOMElement *element, GmxUtil &gmxUtil) const;
};

#endif // GEO_MODEL_XML_MAKE_CONS_H
