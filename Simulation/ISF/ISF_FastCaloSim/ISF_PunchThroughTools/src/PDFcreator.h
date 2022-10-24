/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_PUNCHTHROUGHTOOLS_SRC_PDFCREATOR_H
#define ISF_PUNCHTHROUGHTOOLS_SRC_PDFCREATOR_H

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include <map>

//ISF includes
#include "ISF_FastCaloSimEvent/TFCS1DFunction.h"

namespace CLHEP {
  class HepRandomEngine;
}

namespace ISF
{
  /** @class PDFcreator

      Creates random numbers with a distribution given as ROOT TF1.
      The TF1 function parameters will be retrieved from a histogram given by addPar.

      @author  Elmar Ritsch <Elmar.Ritsch@cern.ch>
      @maintainer/updater Thomas Carter <thomas.michael.carter@cern.ch>
  */

  class PDFcreator
  {

  public:
    /** construct the class with a given TF1 and a random engine */
    PDFcreator() {} ;

    ~PDFcreator() { };

    /** all following is used to set up the class */
    void setName( std::string PDFname ) { m_name = PDFname; }; // set the pdf's name
    void addToEnergyEtaRangeHist1DMap(double energy, std::vector<double> etaMinEtaMax, TFCS1DFunction *hist); //add entry to map linking energy, eta window and histogram
    void addToEnergyEtaRangeHist2DMap(double energy, std::vector<double> etaMinEtaMax, std::map< double , TFCS1DFunction* > *hist); //add entry to map linking energy, eta window and histogram

    /** get the random value with this method, by providing the input parameters */
    double getRand(CLHEP::HepRandomEngine* rndmEngine, const std::vector<double>& inputPar, const double& outEnergy = 0., const double& randMin = 0., const double& randMax = 0.) const;
    std::string getName() const {return m_name;};
    static bool compareEnergy1D(const std::pair< double , std::map< std::vector<double>, TFCS1DFunction*> > map, const double value){ return map.first < value; };
    static bool compareEnergy2D(const std::pair< double , std::map< std::vector<double>, std::map< double , TFCS1DFunction* >* > > map, const double value){ return map.first < value; };
    static bool compareEtaMax1D(const std::pair< std::vector<double>, TFCS1DFunction*> map, const double value){ return map.first.at(1) < value; };
    static bool compareEtaMax2D(const std::pair< std::vector<double>, std::map< double , TFCS1DFunction* >* > map, const double value){ return map.first.at(1) < value; };

  private:
    std::string                         m_name;               //!< Give pdf a name for debug purposes
    std::map< double , std::map< std::vector<double>, TFCS1DFunction*> > m_energy_etaRange_hists1D; //!< map of energies to map of eta ranges to 1D histograms
    std::map< double , std::map< std::vector<double>, std::map< double , TFCS1DFunction* >* > > m_energy_etaRange_hists2D; //!< map of energies to map of eta ranges to 2D histograms
    constexpr static double s_sqrtOf2 = M_SQRT2;

  };
}

#endif
