/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JETCALIBTOOLS_JETSMEARINGCORRECTION_H
#define JETCALIBTOOLS_JETSMEARINGCORRECTION_H


#include "JetCalibTools/JetCalibrationStep.h"

#include <memory>
#include <vector>

#include "boost/thread/tss.hpp"

#include "TRandom3.h"
#include "TEnv.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

class JetSmearingCorrection
    : virtual public JetCalibrationStep
{

    public:
        JetSmearingCorrection();
        JetSmearingCorrection(const std::string& name, TEnv* config, TString jetAlgo, TString calibAreaTag, bool dev);
        virtual ~JetSmearingCorrection();

        virtual StatusCode initialize() override;
        virtual StatusCode calibrate(xAOD::Jet& jet, JetEventInfo&) const override;

        virtual StatusCode getNominalResolutionData(const xAOD::Jet& jet, double& resolution) const override;
        virtual StatusCode getNominalResolutionMC(  const xAOD::Jet& jet, double& resolution) const override;

    private:
        // Helper methods
        StatusCode getSigmaSmear(xAOD::Jet& jet, double& sigmaSmear) const;
        StatusCode getNominalResolution(const xAOD::Jet& jet, TH1* histo, const std::vector< std::unique_ptr<TH1> >& projections, double& resolution) const;
        StatusCode readHisto(double& returnValue, TH1* histo, double x) const;
        StatusCode readHisto(double& returnValue, TH1* histo, const std::vector< std::unique_ptr<TH1> >& projections, double x, double y) const;
        StatusCode cacheProjections(TH1* fullHistogram, std::vector< std::unique_ptr<TH1> >& cacheLocation, const std::string& type);
        TRandom3* getTLSRandomGen(unsigned long seed) const;

        // Private enums
        enum class SmearType
        {
            UNKNOWN=0,  // Unknown/unset/etc
            Pt,         // pT smearing (JpTR)
            Mass,       // Mass smearing (JMR)
            FourVec     // Four-vector smearing (JER)
        };
        enum class HistType
        {
            UNKNOWN=0,  // Unknown/unset/etc
            Pt,         // 1D histogram, binned in pT
            PtEta,      // 2D histogram, binned in x=pT y=eta
            PtAbsEta    // 2D histogram, binned in x=pT y=|eta|
        };
        enum class InterpType
        {
            UNKNOWN=0,  // Unknown/unset/etc
            Full,       // Full interpolation
            None,       // No interpolation
            OnlyX,      // Interpolate only in the x dimension
            OnlyY       // Interpolate only in the y dimension
        };

        // Class variables from constructor
        TEnv* m_config;
        const TString m_jetAlgo;
        const TString m_calibAreaTag;
        const bool m_dev;
        mutable boost::thread_specific_ptr<TRandom3> m_rand_tls; // thread-specific random number generator
        
        // Class variables read in from the config file
        TString m_jetOutScale;
        SmearType  m_smearType;
        HistType   m_histType;
        InterpType m_interpType;
        std::unique_ptr<TH1> m_smearResolutionMC;
        std::unique_ptr<TH1> m_smearResolutionData;

        // Variables to cache projections in case of 1-D interpolaton in 2-D or 3-D histograms
        std::vector< std::unique_ptr<TH1> > m_cachedProjResMC;
        std::vector< std::unique_ptr<TH1> > m_cachedProjResData;
};


#endif

