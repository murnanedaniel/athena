/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef Mdt_H
#define Mdt_H

#include "MuonGeoModel/DetectorElement.h"

#include <string>
#include <vector>

class GeoFullPhysVol;

namespace MuonGM {

class MultiLayer;
class Cutout;
class Component;
class MdtComponent;

class Mdt: public DetectorElement {

public:
   double width;
   double length;
   double thickness;
   double longWidth;	// for trapezoidal layers
   int index;
   double tubelenStepSize;
   double tubePitch;
	
   Mdt(Component* s1, std::string s2);
   ~Mdt();
   MultiLayer* layer;
   GeoFullPhysVol* build();
   GeoFullPhysVol* build(std::vector<Cutout*>);
   void print();

private:
   MdtComponent* m_component;
   Mdt & operator=(const Mdt &right);
   Mdt(const Mdt&);

};

} // namespace MuonGM


#endif
