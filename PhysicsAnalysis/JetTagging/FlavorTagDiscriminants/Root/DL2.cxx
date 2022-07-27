/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "FlavorTagDiscriminants/DL2.h"
#include "FlavorTagDiscriminants/BTagTrackIpAccessor.h"
#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/NanReplacer.hh"

#include "xAODBTagging/BTaggingUtilities.h"

namespace {
  const std::string jetLinkName = "jetLink";
}


namespace FlavorTagDiscriminants {

  // DL2
  //
  // TODO: make this work with more input nodes
  DL2::DL2(const lwt::GraphConfig& graph_config,
           const std::vector<FTagInputConfig>& inputs,
           const std::vector<FTagTrackSequenceConfig>& track_sequences,
           const FTagOptions& options):
    m_jetLink(jetLinkName),
    m_input_node_name(""),
    m_graph(new lwt::LightweightGraph(graph_config,graph_config.outputs.begin()->first)),
    m_variable_cleaner(nullptr),
    m_defaultValue(options.default_output_value)
  {
    // set up inputs
    if (graph_config.inputs.size() > 1) {
      throw std::logic_error("We don't currently support graphs with "
                             "more than one input");
    } else if (graph_config.inputs.size() == 1){
      m_input_node_name = graph_config.inputs.at(0).name;
      m_variable_cleaner.reset(new lwt::NanReplacer(
                                 graph_config.inputs.at(0).defaults,
                                 lwt::rep::all));
    }

    auto [vb, vj, ds] = dataprep::createBvarGetters(inputs);
    m_varsFromBTag = vb;
    m_varsFromJet = vj;
    m_dataDependencyNames += ds;

    auto [tsb, td] = dataprep::createTrackGetters(
      track_sequences, options);
    m_dataDependencyNames += td;
    m_trackSequenceBuilders = tsb;

    auto [decorators, dd] = dataprep::createDecorators(
      graph_config, options);
    m_dataDependencyNames += dd;
    m_decorators = decorators;

    auto [track_validity, is_defaults, ipdd] = dataprep::createIpChecker(
      graph_config, options);
    m_invalid_track_checker = track_validity;
    m_is_defaults = is_defaults;
    m_dataDependencyNames += ipdd;
  }

  void DL2::decorate(const xAOD::BTagging& btag) const {
    auto jetLink = m_jetLink(btag);
    if (!jetLink.isValid()) {
      throw std::runtime_error("invalid jetLink");
    }
    const xAOD::Jet& jet = **jetLink;
    decorate(jet, btag);
  }
  void DL2::decorate(const xAOD::Jet& jet) const {
    decorate(jet, jet);
  }
  void DL2::decorateWithDefaults(const SG::AuxElement& jet) const {
    // save out things
    for (const auto& dec: m_decorators) {
      for (const auto& node: dec.second) {
        // save something that is clearly wrong
        node.second(jet) = m_defaultValue;
      }
    }
  }

  void DL2::decorate(const xAOD::Jet& jet, const SG::AuxElement& btag) const {
    using namespace internal;
    std::vector<NamedVar> vvec;
    for (const auto& getter: m_varsFromBTag) {
      vvec.push_back(getter(btag));
    }
    for (const auto& getter: m_varsFromJet) {
      vvec.push_back(getter(jet));
    }
    std::map<std::string, std::map<std::string, double> > nodes;
    if (m_variable_cleaner) {
      std::map<std::string, double> variables(vvec.begin(), vvec.end());
      auto cleaned = m_variable_cleaner->replace(variables);

      // Note, you can hack in more variables to `cleaned` here.

      // put the cleaned inputs into the node structure
      nodes[m_input_node_name] =  cleaned;
    }

    // add track sequences, check if any are invalid
    char invalid = 0;
    std::map<std::string, std::map<std::string, std::vector<double>>> seqs;
    for (const auto& builder: m_trackSequenceBuilders) {

      Tracks sorted_tracks = builder.tracksFromJet(jet, btag);
      if (m_invalid_track_checker(sorted_tracks)) invalid = 1;
      Tracks flipped_tracks = builder.flipFilter(sorted_tracks, jet);

      for (const auto& seq_builder: builder.sequencesFromTracks) {
        seqs[builder.name].insert(seq_builder(jet, flipped_tracks));
      }
    }

    for (const auto& def: m_is_defaults) {
      def(btag) = invalid;
    }
    if (invalid) {
      decorateWithDefaults(btag);
      return;
    }

    // save out things
    for (const auto& dec: m_decorators) {
      // the second argument to compute(...) is for sequences
      auto out_vals = m_graph->compute(nodes, seqs, dec.first);
      for (const auto& node: dec.second) {
        node.second(btag) = out_vals.at(node.first);
      }
    }
  }

  FTagDataDependencyNames DL2::getDataDependencyNames() const {
    return m_dataDependencyNames;
  }
}
