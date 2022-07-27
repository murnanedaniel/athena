/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUON_MUONLAYERSEGMENTFINDERTOOL_H
#define MUON_MUONLAYERSEGMENTFINDERTOOL_H

#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgTool.h"
#include "CscSegmentMakers/ICscSegmentFinder.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "MuonDetDescrUtils/MuonSectorMapping.h"
#include "MuonHoughPatternTools/MuonLayerHoughTool.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonLayerEvent/MuonSystemExtension.h"
#include "MuonRecHelperTools/MuonEDMPrinterTool.h"
#include "MuonRecToolInterfaces/IMuonPRDSelectionTool.h"
#include "MuonRecToolInterfaces/IMuonSegmentMaker.h"
#include "MuonSegmentMakerToolInterfaces/IMuonClusterSegmentFinder.h"
#include "MuonSegmentMakerToolInterfaces/IMuonClusterSegmentFinderTool.h"
#include "MuonSegmentMakerToolInterfaces/IMuonLayerSegmentFinderTool.h"
#include "MuonRecToolInterfaces/HoughDataPerSec.h"

namespace Muon {

    class MuonSegment;
    struct MuonLayerPrepRawData;
    class MuonLayerROTs;
    class MdtDriftCircleOnTrack;
    class MuonClusterOnTrack;
    class MuonLayerSegmentFinderTool : virtual public IMuonLayerSegmentFinderTool, public AthAlgTool {
    public:
        /** Default AlgTool functions */
        MuonLayerSegmentFinderTool(const std::string& type, const std::string& name, const IInterface* parent);
        virtual ~MuonLayerSegmentFinderTool() = default;
        StatusCode initialize() override;

        /**IMuonLayerSegmentFinderTool interface: find */
        void find(const EventContext& ctx,
                  const MuonSystemExtension::Intersection& intersection, 
                  const  MuonLayerPrepRawData& layerPrepRawData,
                  std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const override;

        void findMdtSegmentsFromHough(const EventContext& ctx,
                                      const MuonSystemExtension::Intersection& intersection, 
                                      std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const override;
    private:
        /** find segments from PRD clusters */
        void findClusterSegments(const EventContext& ctx, const MuonSystemExtension::Intersection& intersection,
                                 const MuonLayerPrepRawData& layerPrepRawData,
                                 std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const;

        /** find csc segments */
        void findCscSegments(const EventContext& ctx, const MuonLayerPrepRawData& layerPrepRawData,
                             std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const;

        /** find mdt segments from hits in the layer */
        void findMdtSegments(const MuonSystemExtension::Intersection& intersection, const MuonLayerPrepRawData& layerPrepRawData,
                             std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const;

        /** find mdt segments main routine */
        void findMdtSegments(const MuonSystemExtension::Intersection& intersection, const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                             const std::vector<const MuonClusterOnTrack*>& clusters,
                             std::vector<std::shared_ptr<const Muon::MuonSegment> >& segments) const;

        ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc{
            this,
            "MuonIdHelperSvc",
            "Muon::MuonIdHelperSvc/MuonIdHelperSvc",
        };

        PublicToolHandle<MuonEDMPrinterTool> m_printer{
            this,
            "MuonEDMPrinterTool",
            "Muon::MuonEDMPrinterTool/MuonEDMPrinterTool",
        };
        ToolHandle<IMuonPRDSelectionTool> m_muonPRDSelectionTool{
            this,
            "MuonPRDSelectionTool",
            "Muon::MuonPRDSelectionTool/MuonPRDSelectionTool",
        };
        ToolHandle<IMuonSegmentMaker> m_segmentMaker{
            this,
            "SegmentMaker",
            "Muon::DCMathSegmentMaker/DCMathSegmentMaker",
        };
        ToolHandle<ICscSegmentFinder> m_csc2dSegmentFinder{
            this,
            "Csc2DSegmentMaker",
            "Csc2dSegmentMaker/Csc2dSegmentMaker",
        };
        ToolHandle<ICscSegmentFinder> m_csc4dSegmentFinder{
            this,
            "Csc4DSegmentMaker",
            "Csc4dSegmentMaker/Csc4dSegmentMaker",
        };
        ToolHandle<IMuonClusterSegmentFinderTool> m_clusterSegMakerNSW{
            this,
            "NSWMuonClusterSegmentFinderTool",
            "Muon::MuonClusterSegmentFinderTool/MuonClusterSegmentFinderTool",
        };
        /// Use the hough data to find sectors in the speectrometer traversed by a muon.
        SG::ReadHandleKey<Muon::HoughDataPerSectorVec> m_houghDataPerSectorVecKey{this, "HoughKey",
                                                                                   "", "HoughDataPerSectorVec key"};

        const Muon::MuonSectorMapping m_muonSectorMapping{};

    };
}  // namespace Muon

#endif
