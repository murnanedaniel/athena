/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TOPEVENTSELECTIONTOOLS_JETFTAGEFFPLOTS_H_
#define TOPEVENTSELECTIONTOOLS_JETFTAGEFFPLOTS_H_

#include <string>

#include "TopEventSelectionTools/EventSelectorBase.h"
#include "TopEventSelectionTools/PlotManager.h"
#include "PMGTools/PMGTruthWeightTool.h"
#include <AsgTools/AnaToolHandle.h>
#include "FTagAnalysisInterfaces/IBTaggingSelectionTool.h"
#include "TopCorrections/ScaleFactorRetriever.h"


class TFile;

namespace EL {
class Worker;
}

namespace top {
class TopConfig;

/**
 * @brief An example of how to quickly make some plots at a certain point in
 * the cutflow.
 */
class JetFtagEffPlots : public EventSelectorBase {
 public:
  
  JetFtagEffPlots(const std::string& name, TFile* outputFile, const std::string& params, std::shared_ptr<top::TopConfig> config,
                 EL::Worker* wk = nullptr);


  bool apply(const top::Event& event) const override;

  
  std::string name() const override;

 private:
  // File units are MeV and normally people like plots in GeV.
  static const double toGeV;

  // Easy access to histograms.
  std::shared_ptr<PlotManager> m_hists               = nullptr;


  // Nominal hash value
  std::size_t m_nominalHashValue;

  std::string m_CDIfile;
  
  bool m_fill_total_hists;

  bool m_use_event_weight;
  
  //optional suffix you can add to your histogram
  std::string m_histogram_suffix;

  // pT and eta bin edges
  std::string m_ptBins;
  std::string m_etaBins;
  
  float m_max_pT;
  float m_min_pT;
  int m_N_pT_bins;

  float m_max_Eta;
  float m_min_Eta;
  int m_N_Eta_bins;

  // name of the jet collection - this is needed for the names of the histograms
  std::string m_jetCollection;
  std::string m_WP;
  std::string m_tagger;
  // shared pointed to instance of TopConfig
  std::shared_ptr<top::TopConfig> m_config;


  asg::AnaToolHandle<IBTaggingSelectionTool> m_selection_tool;
 
  top::ScaleFactorRetriever* m_sfRetriever;

  // function to translate the binnings into vector of bin edges
  void formatBinning(const std::string& str, std::vector<double>& binEdges );

  //helper function to fill histograms
  void FillHistograms(std::shared_ptr<PlotManager> h_ptr, double w_event, const top::Event& event) const;

  
};

}  // namespace top

#endif  // TOPEVENTSELECTIONTOOLS_JETFTAGEFFPLOTS_H_
