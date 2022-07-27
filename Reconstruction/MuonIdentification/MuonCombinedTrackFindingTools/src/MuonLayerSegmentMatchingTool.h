/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUON_MUONLAYERSEGMENTMATCHINGTOOL_H
#define MUON_MUONLAYERSEGMENTMATCHINGTOOL_H

#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "MuonCombinedToolInterfaces/IMuonLayerSegmentMatchingTool.h"
#include "MuonLayerEvent/MuonSystemExtension.h"
#include "MuonRecHelperTools/MuonEDMPrinterTool.h"
#include "MuonSegmentMakerToolInterfaces/IMuonSegmentSelectionTool.h"
#include "MuonSegmentTaggerToolInterfaces/IMuTagMatchingTool.h"
#include "TrkExInterfaces/IExtrapolator.h"

namespace Muon {

    class MuonSegment;

    class MuonLayerSegmentMatchingTool : virtual public Muon::IMuonLayerSegmentMatchingTool, public AthAlgTool {
    public:
        /** Default AlgTool functions */
        MuonLayerSegmentMatchingTool(const std::string& type, const std::string& name, const IInterface* parent);
        virtual ~MuonLayerSegmentMatchingTool() = default;
        StatusCode initialize() override;

        /**IMuonLayerSegmentMatchingTool interface: select */
        void select(const EventContext& ctx, const MuonSystemExtension::Intersection& intersection,
                    const std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments,
                    std::vector<std::shared_ptr<const Muon::MuonSegment> >& selectedSegments) const override;

    private:
        /** match segment to intersection */
        bool match(const EventContext& ctx, const MuonSystemExtension::Intersection& intersection, const MuonSegment& segment) const;

        /// Helper tool for debugging purposes
        PublicToolHandle<MuonEDMPrinterTool> m_printer{this, "MuonEDMPrinterTool", "Muon::MuonEDMPrinterTool/MuonEDMPrinterTool"};

        ToolHandle<Trk::IExtrapolator> m_extrapolator{this, "Extrapolator", "Trk::Extrapolation/AtlasExtrapolator"};

        ToolHandle<IMuTagMatchingTool> m_matchingTool{this, "MatchTool", "MuTagMatchingTool/MuTagMatchingTool"};
    };
}  // namespace Muon

#endif
