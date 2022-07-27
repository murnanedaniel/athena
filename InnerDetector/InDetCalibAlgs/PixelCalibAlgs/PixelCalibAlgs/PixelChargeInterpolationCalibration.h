/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef PixelChargeInterpolationCalibration_h
#define PixelChargeInterpolationCalibration_h

#include "CxxUtils/checker_macros.h"

class string;
template < class T, class Allocator > class vector;
class TDirectory;

namespace PixelCalib{

class PixelChargeInterpolationHistograms;
class PixelChargeInterpolationPlot;
class PixelChargeInterpolationParameters;

class PixelChargeInterpolationCalibration{

public :
	// Modified stuff
	PixelChargeInterpolationCalibration(const std::string& tag,
			PixelChargeInterpolationParameters &parametersModel);
	virtual ~PixelChargeInterpolationCalibration();
	void Fill(int DetType, double GeVTrkPt,
		double alpha, double DeltaRow, double digresphi,
		double resphi, double OmegaPhi,
		double TrkEta, double DeltaCol, double digreseta,
		double reseta, double OmegaEta);
	void Read(TDirectory *readdir = 0);
	void Write(TDirectory *writedir = 0);
	
	int Analyze( const std::string& output, std::vector<std::string> &reference_names);

private:

  PixelChargeInterpolationCalibration(const PixelChargeInterpolationCalibration &);
  PixelChargeInterpolationCalibration &operator=(const PixelChargeInterpolationCalibration&);

	PixelChargeInterpolationHistograms *m_DigitalCalibration;
	PixelChargeInterpolationHistograms *m_AnalogCalibration;
	PixelChargeInterpolationPlot *m_plots;
};

}
#endif
