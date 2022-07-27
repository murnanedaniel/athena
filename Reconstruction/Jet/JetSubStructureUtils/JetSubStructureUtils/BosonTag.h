/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef jetsubstructureutils_bosontag_header
#define jetsubstructureutils_bosontag_header

// @Author: Giordon Stark
// @Email: gstark@cern.ch

// c++ includes
#include <set>
#include <string>

// EDM includes
#include <xAODJet/Jet.h>

// forward-declare the ROOT includes
class TFile;
class TH2;

namespace JetSubStructureUtils {
  class BosonTag {
    public:
      // struct holding the configuration details read in from the recommendations file
      struct CONFIG {
        // vector of coefficients for calculating the mean mass value
        std::vector<float> m_mass_params;
        // variation to apply to mean mass value to calculate mass window
        float m_mass_window;
        // vector of coefficients for calculating the D2 cut value
        std::vector<float> m_d2_params;
        // direction of how to apply D2 cut value (LEFT: cut < val; RIGHT: val < cut)
        std::string m_d2_cut_direction;
        // internal use mostly to signify that the object was successfully configured
        bool m_isConfig;

        // initialize to default values
        CONFIG();
        // load in all values from the file and store in object
        bool setConfigs(const std::vector<float> mass_params, float mass_window, const std::vector<float> d2_params, const std::string d2_cut_direction);
      };

      // standard tool constructor, will read in recommendations file if tagger_alg = "smooth"
      BosonTag( std::string working_point           = "medium",
                std::string tagger_alg              = "smooth",
#ifdef ROOTCORE
                std::string recommendations_file    = "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Wtagging_MC15_Prerecommendations_20150809.dat",
#else
                std::string recommendations_file    = "JetSubStructureUtils/data/config_13TeV_Wtagging_MC15_Prerecommendations_20150809.dat",
#endif
                bool debug                          = false,
                bool verbose                        = false);

      // this is recommended usage, pass in jet, get true/false
      int result(const xAOD::Jet& jet) const;
      // sometimes you don't have certain properties set so pass them in
      //    to select the appropriate tagging recommendation
      int result(const xAOD::Jet& jet, std::string algorithm_name) const;

      // given the jet and configurations, return the corresponding CONFIG object
      std::pair<bool, CONFIG> get_configuration(std::string algorithm_name) const;

      // given the jet and configurations, return the string representation of the jet
      //        eg: AK10LCTRIMF5R20, CA10LCPRUNR50Z15, CA12LCBDRSM100R30Y15
      std::pair<bool, std::string> get_algorithm_name(const xAOD::Jet& jet,
                                     const xAOD::JetAlgorithmType::ID jet_algorithm,
                                     const float size_parameter,
                                     const xAOD::JetInput::Type jet_input,
                                     const xAOD::JetTransform::Type jet_transform) const;

    private:
      std::string m_working_point,
                  m_tagger_alg,
                  m_recommendations_file;
      bool m_debug,
           m_verbose;

      // this is so we don't error out in general, esp. for athena jobs
      bool m_bad_configuration;

      /* map<workingPoint,
             map<tagger,
                 map<algorithm, CONFIG>
                >
            >
      */
      // use (fastjet::JetAlgorithm) jet->getAlgorithmType() later
      // #include <fastjet/JetDefinition.hh>
      std::map<std::string, std::map<std::string, CONFIG>> m_configurations;

      // main 4 details for classifying a jet
      static const SG::AuxElement::ConstAccessor<int> s_AlgorithmType;
      static const SG::AuxElement::ConstAccessor<float> s_SizeParameter;
      static const SG::AuxElement::ConstAccessor<int> s_InputType;
      static const SG::AuxElement::ConstAccessor<int> s_TransformType;

      // for trimming
      static const SG::AuxElement::ConstAccessor<float> s_RClus;
      static const SG::AuxElement::ConstAccessor<float> s_PtFrac;

      // for pruning
      static const SG::AuxElement::ConstAccessor<float> s_RCut;
      static const SG::AuxElement::ConstAccessor<float> s_ZCut;

      // for splitting
      // static const SG::AuxElement::ConstAccessor<int> NSubjetMax ("NSubjetMax");
      static const SG::AuxElement::ConstAccessor<char> s_BDRS;
      /* MuMax, YMin, RClus */
      // static const SG::AuxElement::ConstAccessor<float> RClus ("RClus"); // defined above for trimming
      static const SG::AuxElement::ConstAccessor<float> s_YMin;
      static const SG::AuxElement::ConstAccessor<float> s_MuMax;
      //    under splitting; for Run-1 tagger
      static const SG::AuxElement::ConstAccessor<float> s_YFilt;

      // for D2
      static const SG::AuxElement::ConstAccessor<float> s_D2;
      static const SG::AuxElement::ConstAccessor<float> s_ECF1;
      static const SG::AuxElement::ConstAccessor<float> s_ECF2;
      static const SG::AuxElement::ConstAccessor<float> s_ECF3;

  };
}

#endif
