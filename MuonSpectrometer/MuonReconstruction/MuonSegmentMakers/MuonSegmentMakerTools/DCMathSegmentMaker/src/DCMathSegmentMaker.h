/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DCMATHSEGMENTMAKER_H
#define DCMATHSEGMENTMAKER_H

#include <list>
#include <set>
#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgTool.h"
#include "EventPrimitives/EventPrimitives.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GeoPrimitives/GeoPrimitives.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonPrepRawData/RpcPrepDataContainer.h"
#include "MuonRIO_OnTrack/MuonClusterOnTrack.h"
#include "MuonRecHelperTools/IMuonEDMHelperSvc.h"
#include "MuonRecHelperTools/MuonEDMPrinterTool.h"
#include "MuonRecToolInterfaces/IMdtDriftCircleOnTrackCreator.h"
#include "MuonRecToolInterfaces/IMdtSegmentFinder.h"
#include "MuonRecToolInterfaces/IMuonClusterOnTrackCreator.h"
#include "MuonRecToolInterfaces/IMuonCompetingClustersOnTrackCreator.h"
#include "MuonRecToolInterfaces/IMuonSegmentFittingTool.h"
#include "MuonRecToolInterfaces/IMuonSegmentMaker.h"
#include "MuonSegment/MuonSegment.h"
#include "MuonSegmentMakerInterfaces/IDCSLFitProvider.h"
#include "MuonSegmentMakerToolInterfaces/IMuonSegmentSelectionTool.h"
#include "MuonSegmentMakerToolInterfaces/IMuonSegmentTriggerHitAssociator.h"
#include "TrkDriftCircleMath/Cluster.h"
#include "TrkDriftCircleMath/DCSLFitter.h"
#include "TrkDriftCircleMath/DCSLHitSelector.h"
#include "TrkDriftCircleMath/DCStatistics.h"
#include "TrkDriftCircleMath/DriftCircle.h"
#include "TrkDriftCircleMath/MdtChamberGeometry.h"
#include "TrkFitterInterfaces/ITrackFitter.h"
#include "TrkSurfaces/Surface.h"
#include "MuonStationIntersectCond/MuonIntersectGeoData.h"
namespace Trk {
    class RIO_OnTrack;
    class PlaneSurface;
}  // namespace Trk

namespace TrkDriftCircleMath {
    class MdtMultiChamberGeometry;
    class Segment;
}  // namespace TrkDriftCircleMath

namespace MuonGM {
    class MuonDetectorManager;
}

namespace Muon {
    class MdtPrepData;
}

namespace Muon {

    class MdtDriftCircleOnTrack;

    /**
        Function object to sort MuonClusterOnTrack pointers by Identier.
        Usage:

        std::vector<const MuonClusterOnTrack*> myClusters;

        ...

        std::sort( myClusters.begin(), myClusters.end(), SortClustersById() );
    */
    struct SortClustersById {
        bool operator()(const MuonClusterOnTrack* cl1, const MuonClusterOnTrack* cl2) { return cl1->identify() < cl2->identify(); }
    };

    /**
        Function object to sort pairs containing a double and a pointer to a MuonClusterOnTrack.
        The object is used to sort MuonClusterOnTrack objects along the trajectory of a give segment.
        Usage:

        std::vector<std::pair<double,const MuonClusterOnTrack*> > myDistClusters;

        ...

        std::sort( myDistClusters.begin(), myDistClusters.end(), SortByDistanceToSegment() );
    */
    struct SortByDistanceToSegment {
        bool operator()(const std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>& hit1,
                        const std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>& hit2) {
            return hit1.first < hit2.first;
        }
    };

    /**
        @class DCMathSegmentMaker

        Implementation of a IMuonSegmentMaker.

        For more details look at the mainpage of this package.
    */
    class DCMathSegmentMaker : virtual public IMuonSegmentMaker, public IMuonSegmentTriggerHitAssociator, public AthAlgTool {
    public:
        // pair of eta-phi hits in the same gasgap
        typedef std::pair<std::vector<const MuonClusterOnTrack*>, std::vector<const MuonClusterOnTrack*> > EtaPhiHitsPair;
        // map to sort hit per gasgap
        typedef std::map<Identifier, EtaPhiHitsPair> IdHitMap;
        typedef IdHitMap::iterator IdHitMapIt;
        typedef std::map<Identifier, IdHitMap> ChIdHitMap;
        typedef ChIdHitMap::iterator ChIdHitMapIt;

        typedef std::vector<const MuonClusterOnTrack*>::const_iterator ROTCit;

        struct HitInXZ {
            HitInXZ(Identifier i, bool isM, bool measP, double lx, double lz, double lxmin, double lxmax, double phmin, double phmax) :
                id(i), isMdt(isM), measPhi(measP), x(lx), z(lz), xmin(lxmin), xmax(lxmax), phimin(phmin), phimax(phmax) {}
            Identifier id;
            bool isMdt;
            bool measPhi;
            double x;
            double z;
            double xmin;
            double xmax;
            double phimin;
            double phimax;
        };

        struct Cluster2D {
            /** constructor taking a single phi hit */
            Cluster2D(const Identifier elId, const Identifier ggId, const Amg::Vector2D& lp, double err, const MuonClusterOnTrack* ecl,
                      const MuonClusterOnTrack* pcl) :
                detElId(elId), gasGapId(ggId), locPos(lp), error(err), etaHit(ecl), phiHit(pcl) {
                if (ecl || pcl) { surface().localToGlobal(locPos, Amg::Vector3D::UnitZ(), globalPos); }
                if (pcl) phiHits.push_back(pcl);
            }
            /** constructor taking a vector of phi hits */
            Cluster2D(const Identifier elId, const Identifier ggId, const Amg::Vector2D& lp, double err, const MuonClusterOnTrack* ecl,
                      const std::vector<const MuonClusterOnTrack*>& phs) :
                detElId(elId), gasGapId(ggId), locPos(lp), error(err), etaHit(ecl), phiHits(phs) {
                // if phiHits to empty point phiHit to first hit in PhiHits
                phiHit = phiHits.empty() ? 0 : phiHits.front();
                if (ecl || phiHit) { surface().localToGlobal(locPos, Amg::Vector3D::UnitZ(), globalPos); }
            }
            Identifier detElId;
            Identifier gasGapId;
            Amg::Vector2D locPos;
            double error;  // assume same error for eta and phi
            const MuonClusterOnTrack* etaHit;
            const MuonClusterOnTrack* phiHit;
            std::vector<const MuonClusterOnTrack*> phiHits;
            const Trk::Surface& surface() const {
                if (etaHit)
                    return etaHit->associatedSurface();
                else
                    return phiHit->associatedSurface();
            }
            Identifier identify() const {
                if (etaHit)
                    return etaHit->identify();
                else
                    return phiHit->identify();
            }
            Amg::Vector3D globalPos;
            bool is2D() const { return etaHit && phiHit; }
            bool corrupt() const { return (!etaHit && !phiHit) || error < 0.01; }
        };
        typedef std::vector<Cluster2D> ClusterVec;
        typedef ClusterVec::iterator ClusterIt;
        typedef ClusterVec::const_iterator ClusterCit;
        typedef std::pair<ClusterVec, ClusterVec> ClusterVecPair;

        struct TubeEnds {
            TubeEnds() = default;
            double lxmin{0};
            double lxmax{0};
            double phimin{0};
            double phimax{};
        };

        struct segmentCreationInfo {  // miscellaneous objects needed for segment creation
            segmentCreationInfo(ClusterVecPair& spVecs, TrkDriftCircleMath::MdtMultiChamberGeometry* multiGeo, Amg::Transform3D gToStation,
                                Amg::Transform3D amdbToGlobal, double pmin, double pmax) :
                clusters(spVecs.first, spVecs.second),
                geom(multiGeo),
                globalTrans(gToStation),
                amdbTrans(amdbToGlobal),
                phimin(pmin),
                phimax(pmax) {}
            ClusterVecPair clusters;
            TrkDriftCircleMath::MdtMultiChamberGeometry* geom;
            Amg::Transform3D globalTrans;
            Amg::Transform3D amdbTrans;
            double phimin;
            double phimax;
        };

    public:
        DCMathSegmentMaker(const std::string&, const std::string&, const IInterface*);

        virtual ~DCMathSegmentMaker() = default;

        virtual StatusCode initialize();

        /** find segments starting from a list of RIO_OnTrack objects, implementation of IMuonSegmentMaker interface
          routine.

          Will call:

           std::vector<const MuonSegment*>* find( const Amg::Vector3D& gpos, const Amg::Vector3D& gdir,
                                                    const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                                  const std::vector<const MuonClusterOnTrack*>&  clusters, bool
          hasPhiMeasurements);
        */
        void find(const std::vector<const Trk::RIO_OnTrack*>& rios, Trk::SegmentCollection* segColl = nullptr) const;

        /** find segments starting from a list of RIO_OnTrack objects in multiple chambers, implementation of
           IMuonSegmentMaker interface routine Will call:

           std::vector<const MuonSegment*>* find( const Amg::Vector3D& gpos, const Amg::Vector3D& gdir,
                                                    const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                                  const std::vector<const MuonClusterOnTrack*>&  clusters, bool
           hasPhiMeasurements);
         */
        void find(const std::vector<const Trk::RIO_OnTrack*>& rios1, const std::vector<const Trk::RIO_OnTrack*>& rios2) const;

        /** find segments starting from:
             - a list of MdtDriftCircleOnTrack objects
             - a list of MuonClusterOnTrack objects

             Implementation of IMuonSegmentMaker interface routine

             Will call:

             std::vector<const MuonSegment*>* find( const Amg::Vector3D& gpos, const Amg::Vector3D& gdir,
                                                    const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                                    const std::vector<const MuonClusterOnTrack*>&  clusters,
                                                    bool hasPhiMeasurements, double momentum );
         */
        void find(const std::vector<const MdtDriftCircleOnTrack*>& mdts, const std::vector<const MuonClusterOnTrack*>& clusters,
                  Trk::SegmentCollection* segColl = nullptr) const;

        /** find segments starting from:
            - an estimate of the global position and direction of the particle in the chamber
            - a list of MdtDriftCircleOnTrack
            - a list of MuonClusterOnTrack
            - a boolean to indicate whether the external prediction should be used to set the
              @f$ \phi @f$-direction of the segment
            - an estimate of the momentum of the particle

            The global direction is used to perform a seeded search for segments.
        */
        void find(const Amg::Vector3D& gpos, const Amg::Vector3D& gdir, const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                  const std::vector<const MuonClusterOnTrack*>& clusters, bool hasPhiMeasurements = false,
                  Trk::SegmentCollection* segColl = nullptr, double momentum = 1e9, double sinAngleCut = 0) const;

        /** find segments starting from:
            - a track prediction
            - a list of MdtDriftCircleOnTrack objects in multiple chambers, sorted by chamber
            - a list of MuonClusterOnTrack objects in multiple chambers, sorted by chamber

            Implementation of IMuonSegmentMaker interface routine

            Will call:

            std::vector<const MuonSegment*>* find( const Amg::Vector3D& gpos, const Amg::Vector3D& gdir,
                                                   const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                                   const std::vector<const MuonClusterOnTrack*>&  clusters,
                                                   bool hasPhiMeasurements, double momentum );
        */
        void find(const Trk::TrackRoad& road, const std::vector<std::vector<const MdtDriftCircleOnTrack*> >& mdts,
                  const std::vector<std::vector<const MuonClusterOnTrack*> >& clusters, Trk::SegmentCollection* segColl,
                  bool hasPhiMeasurements = false, double momentum = 1e9) const;

        /** associate trigger hits to an existing segment
            - a segment
            - a list of MuonClusterOnTrack objects
            - a flag indicating whether the eta hits should be added

            returns the new segment, owner ship passed to the caller.

            Implementation of IMuonSegmentTriggerHitAssociator interface routine
        */
       MuonSegment* associateTriggerHits(const MuonSegment& seg, const std::vector<const MuonClusterOnTrack*>& clus,
                                                bool includeEtaHits) const;

    private:
        /** helper routine to print Identifers */
        std::string printSP(std::pair<double, double> resPull, const Cluster2D& spacePoint) const;

        /** apply error scaling for low mometum tracks */
        bool errorScalingRegion(const Identifier& id) const;

        /** calculate error scaling factor */
        double errorScaleFactor(const Identifier& id, double curvature, bool hasPhiMeasurements) const;

        std::vector<Identifier> calculateHoles(const EventContext& ctx, Identifier chid, const Amg::Vector3D& gpos, const Amg::Vector3D& gdir, bool hasMeasuredCoordinate,
                                               std::set<Identifier>& deltaVec, std::set<Identifier>& outoftimeVec,
                                               const std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const;

        TrkDriftCircleMath::DCVec createDCVec(const std::vector<const MdtDriftCircleOnTrack*>& mdts, double errorScale,
                                              std::set<Identifier>& chamberSet, double& phimin, double& phimax,
                                              TrkDriftCircleMath::DCStatistics& dcStatistics, const Amg::Transform3D& gToStation,
                                              const Amg::Transform3D& amdbToGlobal) const;
        ClusterVecPair create1DClusters(const std::vector<const MuonClusterOnTrack*>& clusters) const;
        ClusterVecPair create2DClusters(const std::vector<const MuonClusterOnTrack*>& clusters) const;

        void handleChamber(IdHitMap& gasGapHitMap) const;
        ClusterVecPair createSpacePoints(ChIdHitMap& chIdHitMap) const;
        ClusterVecPair createSpacePoints(IdHitMap& gasGapHitMap) const;
        Cluster2D createSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit, const MuonClusterOnTrack* phiHit) const;
        Cluster2D createRpcSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit,
                                      const std::vector<const MuonClusterOnTrack*>& phiHits) const;
        Cluster2D createTgcSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit, const MuonClusterOnTrack* phiHit) const;
        TrkDriftCircleMath::CLVec createClusterVec(const Identifier& chid, ClusterVec& spVec, const Amg::Transform3D& gToStation) const;

       void associateMDTsToSegment(
            const Amg::Vector3D& gdir, TrkDriftCircleMath::Segment& segment, const std::vector<const MdtDriftCircleOnTrack*>& mdts,
            TrkDriftCircleMath::MdtMultiChamberGeometry* multiGeo, const Amg::Transform3D& gToStation, const Amg::Transform3D& amdbToGlobal,
            std::set<Identifier>& deltaVec, std::set<Identifier>& outoftimeVec,
            std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const;
        std::pair<std::pair<int, int>, bool> associateClustersToSegment(
            const TrkDriftCircleMath::Segment& segment, const Identifier& chid, const Amg::Transform3D& gToStation, ClusterVecPair& spVecs,
            double phimin, double phimax, std::vector<std::pair<double, std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const;
        
        std::unique_ptr<DataVector<const Trk::MeasurementBase>> createROTVec(
            std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const;

        double distanceToSegment(const TrkDriftCircleMath::Segment& segment, const Amg::Vector3D& hitPos,
                                 const Amg::Transform3D& gToStation) const;
        std::pair<double, double> residualAndPullWithSegment(const TrkDriftCircleMath::Segment& segment, const Cluster2D& spacePoint,
                                                             const Amg::Transform3D& gToStation) const;

        TrkDriftCircleMath::MdtChamberGeometry createChamberGeometry(const Identifier& chid, const Amg::Transform3D& gToStation) const;

        const MdtDriftCircleOnTrack* findFirstRotInChamberWithMostHits(const std::vector<const MdtDriftCircleOnTrack*>& mdts) const;

        bool updateSegmentPhi(const Amg::Vector3D& gpos, const Amg::Vector3D& gdir, Amg::Vector2D& segLocPos,
                              Trk::LocalDirection& segLocDir, Trk::PlaneSurface& surf, const std::vector<const Trk::MeasurementBase*>& rots,
                              double phimin, double phimax) const;

        /** rotate the angle XZ until all hits are in bounds.
            returns a pair of bool, double, the bool is false if the code did not update the dXdZ. the double is the updated
           dXdZ */
        std::pair<bool, double> rotateLocalAngleXZIntoBounds(double xline, double zline, double dXdZ, std::vector<HitInXZ>& hits) const;

        /** check whether all hits are in bounds in the XZ plane */
        bool checkBoundsInXZ(double xline, double zline, double dXdZ, std::vector<HitInXZ>& hits) const;

        /** calculate positions of tube ends */
        TubeEnds localTubeEnds(const MdtDriftCircleOnTrack& mdt, const Amg::Transform3D& gToSegment,
                               const Amg::Transform3D& segmentToG) const;

        /** update phi ranges */
        void updatePhiRanges(double phiminNew, double phimaxNew, double& phiminRef, double& phimaxRef) const;

        /** check whether phi is consistent with segment phi */
        bool checkPhiConsistency(double phi, double phimin, double phimax) const;

        /** update the global direction, keeping the phi of the input road direction but using the local angle YZ */
        Amg::Vector3D updateDirection(double linephi, const Trk::PlaneSurface& surf, const Amg::Vector3D& roaddir,
                                      bool isCurvedSegment) const;

        void addEtaHits(std::vector<const MuonClusterOnTrack*>& clusters,
                        std::vector<std::unique_ptr<const Trk::MeasurementBase> >& tash_bin, bool isEndcap) const;

        std::unique_ptr<MuonSegment> createSegment(const EventContext& ctx, TrkDriftCircleMath::Segment& segment, const Identifier& chid, const Amg::Vector3D& roadpos,
                                   const Amg::Vector3D& roaddir2, const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                   bool hasPhiMeasurements, segmentCreationInfo& sInfo) const;

        const MdtPrepData* findMdt(const Identifier& id) const;

        /** pointers to IdHelpers */
        SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey{
            this,
            "DetectorManagerKey",
            "MuonDetectorManager",
            "Key of input MuonDetectorManager condition data",
        };

        ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc{
            this,
            "MuonIdHelperSvc",
            "Muon::MuonIdHelperSvc/MuonIdHelperSvc",
        };
        ServiceHandle<IMuonEDMHelperSvc> m_edmHelperSvc{
            this,
            "edmHelper",
            "Muon::MuonEDMHelperSvc/MuonEDMHelperSvc",
            "Handle to the service providing the IMuonEDMHelperSvc interface",
        };  //<! edm helper tool

        ToolHandle<IMdtDriftCircleOnTrackCreator> m_mdtCreator{
            this,
            "MdtCreator",
            "Muon::MdtDriftCircleOnTrackCreator/MdtDriftCircleOnTrackCreator",
        };  //<! mdt rio ontrack creator
        ToolHandle<IMdtDriftCircleOnTrackCreator> m_mdtCreatorT0{
            this,
            "MdtCreatorT0",
            "Muon::MdtDriftCircleOnTrackCreator/MdtDriftCircleOnTrackCreator",
        };  //<! mdt rio ontrack creator
        ToolHandle<IMuonClusterOnTrackCreator> m_clusterCreator{
            this,
            "MuonClusterCreator",
            "Muon::MuonClusterOnTrackCreator/MuonClusterOnTrackCreator",
        };  //<! cluster rio ontrack creator
        ToolHandle<IMuonCompetingClustersOnTrackCreator> m_compClusterCreator{
            this,
            "MuonCompetingClustersCreator",
            "Muon::TriggerChamberClusterOnTrackCreator/TriggerChamberClusterOnTrackCreator",
        };  //<! competing clusters rio ontrack creator
        PublicToolHandle<MuonEDMPrinterTool> m_printer{
            this,
            "EDMPrinter",
            "Muon::MuonEDMPrinterTool/MuonEDMPrinterTool",
        };  //<! printer helper tool
        ToolHandle<IMdtSegmentFinder> m_segmentFinder{
            this,
            "MdtSegmentFinder",
            "Muon::MdtMathSegmentFinder/MdtMathSegmentFinder",
        };  //<! segment finder tool
        ToolHandle<IMuonSegmentFittingTool> m_segmentFitter{
            this,
            "SegmentFitter",
            "Muon::MuonSegmentFittingTool/MuonSegmentFittingTool",
        };  //<! segment fitting tool
        ToolHandle<IMuonSegmentSelectionTool> m_segmentSelectionTool{
            this,
            "SegmentSelector",
            "Muon::MuonSegmentSelectionTool/MuonSegmentSelectionTool",
        };  //<! segment selection tool
        ToolHandle<IDCSLFitProvider> m_dcslFitProvider{
            this,
            "DCFitProvider",
            "",
        };

        Gaudi::Property<double> m_sinAngleCut{this, "SinAngleCut", 0.2};                        //<! cut on the angle between the segment and the prediction
        Gaudi::Property<bool> m_doGeometry{this, "DoGeometry", true};                           //<! use chamber geometry in segment finding
        Gaudi::Property<bool> m_curvedErrorScaling{this, "CurvedErrorScaling", true};                   //<! rescale errors for low momenta
        Gaudi::Property<bool> m_doSpacePoints{this, "UseTriggerSpacePoints", true};                        //<! use cluster space points for association
        Gaudi::Property<bool> m_createCompetingROTsEta{this, "CreateCompetingROTsEta", true};               //<! create competing ROTs for the clusters
        Gaudi::Property<bool> m_createCompetingROTsPhi{this, "CreateCompetingROTsPhi", true};               //<! create competing ROTs for the clusters
        Gaudi::Property<bool> m_refitParameters{this, "RefitSegment", false};                      //<! refit segment if there are sufficient phi hits and update the segment parameters
        Gaudi::Property<bool> m_addUnassociatedPhiHits{this, "AddUnassociatedPhiHits", false};               //<! if there are phi hits without associated eta hit add them to segment
        Gaudi::Property<bool> m_strictRoadDirectionConsistencyCheck{this, "StrictRoadDirectionConsistencyCheck", true};  //<! check if direction of road is consistent with IP (default: true),
                                                     // should be off for cosmics
        Gaudi::Property<double> m_maxAssociateClusterDistance{this, "MaxAssociateClusterDistance", 3000.};        //<! maximum distance for clusters to be associated to segment (default: 3000
                                                     //(mm))
        Gaudi::Property<bool> m_allMdtHoles{this, "AllMdtHoles", false};                          //<! add all mdt holes without bound checks / flag to decide whether to apply bound checks during the hole search
        Gaudi::Property<bool> m_removeDeltas{this, "RemoveDeltasFromSegmentQuality", true};                         //<! do not add delta electrons to MuonSegmentQuality::holes
        Gaudi::Property<bool> m_reject1DTgcSpacePoints{this,"Reject1DTgcSpacePoints", true };               //<! remove 1D tgc space points / reject TGC eta hits that are not associated with a phi hit in the same gas gap
        Gaudi::Property<bool> m_usePreciseError{this, "UsePreciseError", false};
        Gaudi::Property<bool> m_outputFittedT0{this, "OutputFittedT0", false};
        Gaudi::Property<double> m_preciseErrorScale{this, "PreciseErrorScale", 2.};
        Gaudi::Property<bool> m_doTimeOutChecks{this, "UseTimeOutGard", false};
        
        Gaudi::Property<bool> m_recoverBadRpcCabling{this, "RecoverBadRpcCabling", false};
        Gaudi::Property<bool> m_updatePhiUsingPhiHits{this, "UpdatePhiUsingPhiHits", false};
        Gaudi::Property<bool> m_assumePointingPhi{this, "AssumePointingPhi", false };
        Gaudi::Property<bool> m_redo2DFit{this, "Redo2DFit", true};




        SG::ReadHandleKey<Muon::RpcPrepDataContainer> m_rpcKey{this, "RpcPrepDataContainer", "RPC_Measurements"};
        SG::ReadHandleKey<Muon::TgcPrepDataContainer> m_tgcKey{this, "TgcPrepDataContainer", "TGC_Measurements"};
        SG::ReadHandleKey<Muon::MdtPrepDataContainer> m_mdtKey{this, "MdtPrepDataContainer", "MDT_DriftCircles"};
        
        SG::ReadCondHandleKey<Muon::MuonIntersectGeoData> m_chamberGeoKey{this, "ChamberGeoKey", "MuonStationIntersects", "Pointer to hole search service"};
   
        
    };

}  // namespace Muon
#endif
