/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef sTGC_Technology_H
#define sTGC_Technology_H

#include "AGDDKernel/AGDDTechnology.h"
#include <vector>
namespace MuonGM {

// Description class to build sTGC chambers

class sTGC_Technology: public AGDDTechnology {
public:
	double thickness;
	int nlayers;
	double gasThickness;
	double pcbThickness; // Included for Backwards Compatibility
	double pcbThickness150;
	double pcbThickness200;
	double coverThickness;
	double f4Thickness;
  double f5Thickness;
  double f6Thickness;



	int geoLevel;

	// inner structure parameters (to be defined)

	// constructor
	inline sTGC_Technology(const std::string& s,
                               AGDDDetectorStore& ds);
	inline double Thickness() const;
};

sTGC_Technology::sTGC_Technology(const std::string& s,
                                 AGDDDetectorStore& ds): AGDDTechnology(s, ds),
  geoLevel(0)
{
}

double sTGC_Technology::Thickness() const
{
	//thickness=nlayers*(gasThickness+pcbThickness) + pcbThickness;
	return thickness;
}

} // namespace MuonGM

#endif
