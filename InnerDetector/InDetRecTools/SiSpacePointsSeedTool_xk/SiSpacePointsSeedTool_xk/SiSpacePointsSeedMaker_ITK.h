// -*- C++ -*-

/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////////////////////
//  Header file for class SiSpacePointsSeedMaker_ITK
/////////////////////////////////////////////////////////////////////////////////
// Version 1.0 3/10/2004 I.Gavrilenko
/////////////////////////////////////////////////////////////////////////////////

#ifndef SiSpacePointsSeedMaker_ITK_H
#define SiSpacePointsSeedMaker_ITK_H

#include "InDetRecToolInterfaces/ISiSpacePointsSeedMaker.h"
#include "AthenaBaseComps/AthAlgTool.h"

#include "BeamSpotConditionsData/BeamSpotData.h"
#include "SiSPSeededTrackFinderData/SiSpacePointForSeedITK.h"
#include "SiSPSeededTrackFinderData/SiSpacePointsSeedMakerEventData.h"
#include "TrkSpacePoint/SpacePointContainer.h" 
#include "TrkSpacePoint/SpacePointOverlapCollection.h"
#include "TrkEventUtils/PRDtoTrackMap.h"

//for validation
#include "GaudiKernel/ITHistSvc.h"
#include "TFile.h"
#include "TTree.h"




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MagField cache
#include "MagFieldConditions/AtlasFieldCacheCondObj.h"
#include "MagFieldElements/AtlasFieldCache.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <list>
#include <vector>

class MsgStream;


namespace InDet {

  ///////////////////////////////////////////////////////////////////
  // Object function for ordering space point in R coordinate order
  ///////////////////////////////////////////////////////////////////
  
  class SiSpacePointsITKComparison_R {
 
  public: 
    bool operator () (InDet::SiSpacePointForSeedITK* s1,InDet::SiSpacePointForSeedITK* s2) {
      return((*s1).radius() < (*s2).radius());
    }
  };


  using EventData = SiSpacePointsSeedMakerEventData;

  /**
   * @class SiSpacePointsSeedMaker_ATLxk
   * Class for track candidates generation using space points information
   * for standard Atlas geometry
   *
   * In AthenaMT, event dependent cache inside SiSpacePointsSeedMaker_ITK
   * is not preferred. SiSpacePointsSeedMakerEventData = EventData class
   * holds event dependent data for SiSpacePointsSeedMaker_ITK.
   * Its object is instantiated in SiSPSeededTrackFinder::execute.
   */

  class SiSpacePointsSeedMaker_ITK : 
    public extends<AthAlgTool, ISiSpacePointsSeedMaker>
  {
    ///////////////////////////////////////////////////////////////////
    // Public methods:
    ///////////////////////////////////////////////////////////////////
      
  public:
      
    ///////////////////////////////////////////////////////////////////
    /// @name Standard tool methods
    ///////////////////////////////////////////////////////////////////
    //@{
    SiSpacePointsSeedMaker_ITK
    (const std::string&,const std::string&,const IInterface*);
    virtual ~SiSpacePointsSeedMaker_ITK() = default;
    virtual StatusCode initialize() override;
    virtual StatusCode finalize() override;
    //@}

    ///////////////////////////////////////////////////////////////////
    /// @name Methods to initialize tool for new event or region
    ///////////////////////////////////////////////////////////////////
    //@{
    virtual void newEvent (const EventContext& ctx, EventData& data, int iteration) const override;
    virtual void newRegion(const EventContext& ctx, EventData& data,
                           const std::vector<IdentifierHash>& vPixel, const std::vector<IdentifierHash>& vSCT) const override;
    virtual void newRegion(const EventContext& ctx, EventData& data,
                           const std::vector<IdentifierHash>& vPixel, const std::vector<IdentifierHash>& vSCT,
                           const IRoiDescriptor& iRD) const override;
    //@}

    ///////////////////////////////////////////////////////////////////
    /// @name Methods to initilize different strategies of seeds production
    ///////////////////////////////////////////////////////////////////
    //@{

    /// with two space points with or without vertex constraint
    virtual void find2Sp(EventData& data, const std::list<Trk::Vertex>& lv) const override;

    /// with three space points with or without vertex constraint
    virtual void find3Sp(const EventContext& ctx, EventData& data, const std::list<Trk::Vertex>& lv) const override;

    /// with three space points with or without vertex constraint
    /// with information about min and max Z of the vertex
    virtual void find3Sp(const EventContext& ctx, EventData& data, const std::list<Trk::Vertex>& lv, const double* zVertex) const override;

    /// with variable number space points with or without vertex constraint
    /// Variable means (2,3,4,....) any number space points
    virtual void findVSp(const EventContext& ctx, EventData& data, const std::list<Trk::Vertex>& lv) const override;
    //@}
      
    ///////////////////////////////////////////////////////////////////
    /// @name Iterator through seeds pseudo collection
    /// produced accordingly methods find    
    ///////////////////////////////////////////////////////////////////
    //@{
    virtual const SiSpacePointsSeed* next(const EventContext& ctx, EventData& data) const override;
    //@}

    virtual void writeNtuple(const SiSpacePointsSeed* seed, const Trk::Track* track, int seedType, long eventNumber) const override;
    virtual bool getWriteNtupleBoolProperty() const override;

    ///////////////////////////////////////////////////////////////////
    /// @name Print internal tool parameters and status
    ///////////////////////////////////////////////////////////////////
    //@{
    virtual MsgStream& dump(EventData& data, MsgStream& out) const override;
    //@}

  private:
    /// enum for array sizes
    /// Note that this stores the maximum capacities, the actual binnings 
    /// do not always use the full size. See data members below for the 
    /// actual binning paramaters, which are determined in buildFramework. 
    enum Size {arraySizePhi=200,     ///< capacity of the 1D phi arrays 
                arraySizeZ=11,       ///< capacity of the 1D z arrays
                arraySizePhiZ=arraySizePhi*arraySizeZ,   ///< capacity for the 2D phi-z arrays 
                arraySizeNeighbourBins=9,  ///< array size to store neighbouring phi-z-regions in the seed finding
                arraySizePhiV=100,         ///< array size in phi for vertexing 
                arraySizeZV=3,             ///< array size in z for vertexing
                arraySizePhiZV=arraySizePhiV*arraySizeZV,      ///< array size in phi-Z 2D for the vertexing
                arraySizeNeighbourBinsVertex=6};       ///< array size to store neighbouring phi-z regions for the vertexing


    ///////////////////////////////////////////////////////////////////
    // Private data and methods
    ///////////////////////////////////////////////////////////////////

    /// @name Data handles
    //@{
    SG::ReadHandleKey<SpacePointContainer> m_spacepointsSCT{this, "SpacePointsSCTName", "SCT_SpacePoints", "SCT space points container"};
    SG::ReadHandleKey<SpacePointContainer> m_spacepointsPixel{this, "SpacePointsPixelName", "PixelSpacePoints", "Pixel space points container"};
    SG::ReadHandleKey<SpacePointOverlapCollection> m_spacepointsOverlap{this, "SpacePointsOverlapName", "OverlapSpacePoints"};
    SG::ReadHandleKey<Trk::PRDtoTrackMap> m_prdToTrackMap{this,"PRDtoTrackMap","","option PRD-to-track association"};
    SG::ReadCondHandleKey<InDet::BeamSpotData> m_beamSpotKey{this, "BeamSpotKey", "BeamSpotData", "SG key for beam spot"};
    SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_fieldCondObjInputKey {this, "AtlasFieldCacheCondObj", "fieldCondObj",
                                                                           "Name of the Magnetic Field conditions object key"};
    //@}

    /// @name Properties, which will not be changed after construction
    //@{
    BooleanProperty m_pixel{this, "usePixel", true};
    BooleanProperty m_sct{this, "useSCT", true};
    BooleanProperty m_useOverlap{this, "useOverlapSpCollection", true};
    IntegerProperty m_maxsize{this, "maxSize", 50000};
    IntegerProperty m_maxsizeSP{this, "maxSizeSP", 5000};
    IntegerProperty m_maxOneSize{this, "maxSeedsForSpacePoint", 5};
    FloatProperty m_etamax{this, "etaMax", 2.7};
    FloatProperty m_r1minv{this, "minVRadius1", 0.};
    FloatProperty m_r1maxv{this, "maxVRadius1", 60.};
    FloatProperty m_r2minv{this, "minVRadius2", 70.};
    FloatProperty m_r2maxv{this, "maxVRadius2", 200.};
    FloatProperty m_drmin{this, "mindRadius", 5.};
    FloatProperty m_drmax{this, "maxdRadius", 300.};
    FloatProperty m_zmin{this, "minZ", -250.};
    FloatProperty m_zmax{this , "maxZ", +250.};
    FloatProperty m_r_rmin{this, "radMin", 0.};
    FloatProperty m_binSizeR{this, "radStep", 2.};
    FloatProperty m_dzver{this, "maxdZver", 5.};
    FloatProperty m_dzdrver{this, "maxdZdRver", 0.02};
    FloatProperty m_maxdImpact{this, "maxdImpact", 10.};
    FloatProperty m_maxdImpactSSS{this, "maxdImpactSSS", 20.};
    FloatProperty m_divermax{this, "maxdImpactForDecays", 20.};
    FloatProperty m_dzmaxPPP{this, "dZmaxForPPPSeeds", 600.};

    FloatProperty m_maxScore{this, "maximumAcceptedSeedScore", 100.};
    IntegerProperty m_maxOneSizeSSS{this, "maxSeedsForSpacePointStrips", 5};
    IntegerProperty m_maxOneSizePPP{this, "maxSeedsForSpacePointPixels", 5};
    BooleanProperty m_alwaysKeepConfirmedPixelSeeds{this, "alwaysKeepConfirmedPixelSeeds", false};
    BooleanProperty m_alwaysKeepConfirmedStripSeeds{this, "alwaysKeepConfirmedStripSeeds", false};
    BooleanProperty m_fastTracking{this, "useFastTracking", false};
    FloatProperty m_seedScoreBonusPPP{this, "seedScoreBonusPPP", -200.};
    FloatProperty m_seedScoreBonusSSS{this, "seedScoreBonusSSS", -400.};
    BooleanProperty m_optimisePhiBinning{this, "optimisePhiBinning", true};
    FloatProperty m_rminSSS{this, "radMinSSS", 400.};
    FloatProperty m_rmaxSSS{this, "radMaxSSS", 1000.};
    BooleanProperty m_isLRT{this, "isLRT", false};
    FloatProperty m_drminPPP{this, "mindRadiusPPP", 6.};
    FloatProperty m_drmaxPPP{this, "maxdRadiusPPP", 140.};
    FloatProperty m_drminSSS{this, "mindRadiusSSS", 20.};
    FloatProperty m_drmaxSSS{this, "maxdRadiusSSS", 3000.};
    FloatProperty m_dImpactCutSlopeUnconfirmedSSS{this, "dImpactCutSlopeUnconfirmedSSS", 1.0};
    FloatProperty m_dImpactCutSlopeUnconfirmedPPP{this, "dImpactCutSlopeUnconfirmedPPP", 0.};
    FloatProperty m_seedScoreBonusConfirmationSeed{this, "seedScoreBonusConfirmationSeed", -200.};
    BooleanProperty m_useSeedConfirmation{this, "useSeedConfirmation", false};
    FloatProperty m_rminPPPFast{this, "m_rminPPPFast", 70.};
    //@}

    /// @name Properties, which will be updated in initialize
    //@{
    FloatProperty m_etamin{this, "etaMin", 0.};
    FloatProperty m_r_rmax{this, "radMax", 1100.};
    FloatProperty m_ptmin{this, "pTmin", 500.};
    FloatProperty m_umax{this, "minSeedsQuality", 0.};
    //@}

    /// @name Properties, which will be updated in event methods, checketa is prepared in EventData.
    //@{
    BooleanProperty m_checketa{this, "checkEta", false};
    //@}

    /// @name Properties, which are not used in this implementation of SiSpacePointsSeedMaker_ITK class
    //@{
    BooleanProperty m_dbm{this, "useDBM", false};
    UnsignedIntegerProperty m_maxNumberVertices{this, "maxNumberVertices", 99};
    FloatProperty m_r1min{this, "minRadius1", 0.};
    FloatProperty m_r1max{this, "maxRadius1", 600.};
    FloatProperty m_r2min{this, "minRadius2", 0.};
    FloatProperty m_r2max{this, "maxRadius2", 600.};
    FloatProperty m_r3min{this, "minRadius3", 0.};
    FloatProperty m_r3max{this, "maxRadius3", 600.};
    FloatProperty m_rapcut{this, "RapidityCut", 2.7};
    FloatProperty m_diverpps{this, "maxdImpactPPS", 1.7};
    //@}

    /// @name Data member, which is not updated at all
    //@{
    float m_drminv{20.};
    //@}

    /// @name Data members, which are updated only in initialize
    //@{
    bool m_initialized{false};
    int m_outputlevel{0};
    
    int m_fNmax{0};
    int m_fvNmax{0};
    int m_rfz_b[arraySizePhiZ]{};
    int m_rfz_t[arraySizePhiZ]{};
    int m_rfz_ib[arraySizePhiZ][arraySizeNeighbourBins]{};
    int m_rfz_it[arraySizePhiZ][arraySizeNeighbourBins]{};
    int m_rfzv_n[arraySizePhiZV]{};
    int m_rfzv_i[arraySizePhiZV][arraySizeNeighbourBinsVertex]{};
    float m_dzdrmin0{0.};
    float m_dzdrmax0{0.};
    float m_ipt{0.};
    float m_ipt2{0.};
    float m_COF{0.};
    float m_sF{0.};
    float m_sFv{0.};
    float m_dzMaxFast   {200.};
    float m_R2MaxFast   {2500.};        
    float m_zmaxPPP     {2700.};
    float m_rmaxPPP     {140.};   
    float m_zmaxSSS     {2700.};
    float m_dzmaxSSS    {900.};    
    float m_drminSeedConf{5.};
    //@}

     /// @name Properties of writeNTuple for validation purpose
    //@{
    Gaudi::Property<bool> m_writeNtuple {this, "WriteNtuple", false, "Flag to write Validation Ntuples"};
    ///Flag to write validation ntuples. Turned off by default
    ITHistSvc* m_thistSvc;
    TTree* m_outputTree;
    mutable std::mutex m_mutex;
    mutable std::string          m_treeName               ATLAS_THREAD_SAFE;
    mutable TString              m_treeFolder             ATLAS_THREAD_SAFE;
    mutable float                  m_d0                   ATLAS_THREAD_SAFE = 0;
    mutable float                  m_z0                   ATLAS_THREAD_SAFE = 0;
    mutable float                  m_pt                   ATLAS_THREAD_SAFE = 0;
    mutable float                  m_eta                  ATLAS_THREAD_SAFE = 0;
    mutable double                 m_x1                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_x2                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_x3                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_y1                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_y2                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_y3                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_z1                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_z2                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_z3                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_r1                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_r2                   ATLAS_THREAD_SAFE = 0;
    mutable double                 m_r3                   ATLAS_THREAD_SAFE = 0;
    mutable float                  m_quality              ATLAS_THREAD_SAFE = 0;
    mutable int                    m_type                 ATLAS_THREAD_SAFE = 0;
    mutable double                 m_dzdr_t               ATLAS_THREAD_SAFE = 0;
    mutable double                 m_dzdr_b               ATLAS_THREAD_SAFE = 0;
    mutable bool                   m_givesTrack           ATLAS_THREAD_SAFE = 0;
    mutable float                  m_trackPt              ATLAS_THREAD_SAFE = 0;
    mutable float                  m_trackEta             ATLAS_THREAD_SAFE = 0;
    mutable long                   m_eventNumber          ATLAS_THREAD_SAFE = 0;
    //@}


    /// @name Binning parameters 
    int m_nBinsR{0};              ///<  number of bins in the radial coordinate 
    int m_maxPhiBinPPP{0};           ///<  number of bins in phi 
    float m_inverseBinSizePhiPPP{0};   ///<  cache the inverse bin size in phi which we use - needed to evaluate phi bin locations
    int m_maxPhiBinSSS{0};
    float m_inverseBinSizePhiSSS{0};

    
    /** Seed score thresholds defined based on the modifiers defined 
      * as configurables above.
      * These allow to categorise the seeds based on their quality score.
      * The modifiers above are much larger than the range of the raw 
      * unmodified scores would be, resulting in a grouping of the scores
      * based on the modifiers applied. The thresholds are just below 
      * the value reachable without the additional modifier
    **/   
    float m_seedScoreThresholdPPPConfirmationSeed{0.};    ///< max (score is assigned negative sign) score for PPP seeds with confirmation seed requirement. 
    float m_seedScoreThresholdSSSConfirmationSeed{0.};    ///< max (score is assigned negative sign) score for SSS seeds with confirmation seed requirement. 


    /// arrays associating bins to each other for SP formation
    std::array<int,arraySizePhiZ> m_nNeighbourCellsBottomPPP;  ///< number of neighbouring phi-z bins to consider when looking for "bottom SP" candidates for each phi-z bin
    std::array<int,arraySizePhiZ> m_nNeighbourCellsTopPPP;  ///< number of neighbouring phi-z bins to consider when looking for "top SP" candidates for each phi-z bin
    std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ> m_neighbourCellsBottomPPP; ///< mapping of neighbour cells in the 2D phi-z binning to consider  for the "bottom SP" search for central SPs in each phi-z bin. Number of valid entries stored in m_nNeighboursPhiZbottom
    std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ> m_neighbourCellsTopPPP; ///< mapping of neighbour cells in the 2D phi-z binning to consider  for the "top SP" search for central SPs in each phi-z bin. Number of valid entries stored in m_nNeighboursPhiZtop

    std::array<int,arraySizePhiZ> m_nNeighbourCellsBottomSSS;
    std::array<int,arraySizePhiZ> m_nNeighbourCellsTopSSS;
    std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ> m_neighbourCellsBottomSSS;
    std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ> m_neighbourCellsTopSSS;


    ///////////////////////////////////////////////////////////////////
    // Private methods
    ///////////////////////////////////////////////////////////////////
    /// @name Disallow default instantiation, copy, assignment
    //@{
    SiSpacePointsSeedMaker_ITK() = delete;
    SiSpacePointsSeedMaker_ITK(const SiSpacePointsSeedMaker_ITK&) = delete;
    SiSpacePointsSeedMaker_ITK &operator=(const SiSpacePointsSeedMaker_ITK&) = delete;
    //@}

    MsgStream& dumpConditions(EventData& data, MsgStream& out) const;
    MsgStream& dumpEvent(EventData& data, MsgStream& out) const;

    void buildFrameWork();

    void buildConnectionMaps(std::array<int, arraySizePhiZ>& nNeighbourCellsBottom,
			       std::array<int, arraySizePhiZ>& nNeighbourCellsTop,
			       std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ>& neighbourCellsBottom,
			       std::array<std::array<int, arraySizeNeighbourBins>, arraySizePhiZ>& neighbourCellsTop,
			       int maxPhiBin, bool isSSS);

    void buildBeamFrameWork(EventData& data) const;

    /** Determine the expected azimuthal trajectory displacement
      *  in phi in presence of the magnetic field for a particle 
      *  with momentum pTmin and impact parameter maxd0, 
      *  moving from a radial coordinate Rmin outward to Rmax.
      *  This method is used to determine the optimal binning 
      *  of the phi-z regions we consider in the seed making,
      *  to ensure we contain the hits from our softest tracks
      *  in a set of consecutive bins.
      *  @param[in] pTmin: minimum pt cut applied in MeV
      *  @param[in] maxD0: maximum d0 allowed 
      *  @param[in] Rmin: starting radius for trajectory displacement
      *  @param[in] Rmax: end radius for trajectory displacement
    **/
    float azimuthalStep(const float pTmin,const float maxd0,const float Rmin,const float Rmax) const;


    /** Create a SiSpacePointForSeed from the space point. 
      * This will also add the point to the data object's
      * l_spforseed list and update its i_spforseed iterator 
      * to point to the entry after the new SP 
      * for further additions.
      * Returns a nullptr if the SP fails the eta cut, 
      * should we apply one 
      * @param[in,out] data: Provides beam spot location, receives updates to the l_spforseed and i_spforseed members 
      * @param[in] sp: Input space point. 
    **/
    SiSpacePointForSeedITK* newSpacePoint(EventData& data, const Trk::SpacePoint*const& sp) const;
    SiSpacePointForSeedITK* newSpacePoint(EventData& data, const Trk::SpacePoint*const& sp, float* r, bool usePixSctInform=false) const;

    void newSeed
    (EventData& data,
     SiSpacePointForSeedITK*&,SiSpacePointForSeedITK*&,float) const;

    void newOneSeed
    (EventData& data, 
     SiSpacePointForSeedITK*&,SiSpacePointForSeedITK*&,
     SiSpacePointForSeedITK*&,float,float) const;

    void newOneSeedWithCurvaturesComparison
    (EventData& data,
     SiSpacePointForSeedITK*&,SiSpacePointForSeedITK*&,float) const;

    void fillSeeds(EventData& data) const;
    void fillLists(EventData& data) const;
    void fillListsFast(const EventContext& ctx, EventData& data) const;
    void pixInform(const Trk::SpacePoint* sp, float* r) const;
    void sctInform(EventData& data,const Trk::SpacePoint* sp, float* r) const;
    void erase(EventData& data) const;
    void production2Sp(EventData& data) const;
    void production3Sp(EventData& data) const;

    /** \brief: Seed production from space points. 
       * 
       * This method will try to find 3-SP combinations within a 
       * local phi-z region in the detector. 
       * 
       * The central SP of the seed will be taken from this region
       * (technically via the first entry of the bottom candidate array, 
       * which always points to the phi-z bin of interest itself). 
       * 
       * The top SP is allowed to come from the same or one of several close-by
       * phi-Z bins, as is the bottom SP. 
       * 
       * All SP collections are expected to be internally sorted in the radial coordinate.
       * 
       * @param[in,out] data: Event data
       * @param[in,out] iter_bottomCands: collection of iterators over SP collections for up to 9 phi-z cells to consider for the bottom space-point search 
       * @param[in,out] iter_endBottomCands: collection of end-iterators over the 
       * SP collections  for up to 9 phi-z cells to consider for the bottom space-point search 
       * @param[in,out] iter_topCands: collection of iterators over SP collections for up to 9 phi-z cells to consider for the top space-point search 
       * @param[in,out] iter_endTopCands: collection of end-iterators over the 
       * SP collections  for up to 9 phi-z cells to consider for the top space-point search 
       * @param[in] numberBottomCells: Number of bottom cells to consider. Determines how many entries in iter_(end)bottomCands are expected to be valid. 
       * @param[in] numberTopCells: Number of top cells to consider.Determines how many entries in iter_(end)topCands are expected to be valid. 
       * @param[out] nseed: Number of seeds found 
       **/ 
      void production3SpSSS
      (EventData& data,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_bottomCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_endBottomCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_topCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_endTopCands,
      const int numberBottomCells, const int numberTopCells, int& nseed) const;

      void production3SpPPP
      (EventData& data,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_bottomCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_endBottomCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_topCands,
      std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & iter_endTopCands,
      const int numberBottomCells, const int numberTopCells, int& nseed) const;

      /// as above, but for the trigger 
      void production3SpTrigger
      (EventData& /*data*/,
       std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & /*rb*/,
       std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & /*rbe*/,
       std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & /*rt*/,
       std::array<std::list<InDet::SiSpacePointForSeedITK*>::iterator, arraySizeNeighbourBins> & /*rte*/,
       const int /*numberBottomCells*/, const int /*numberTopCells*/, int& /*nseed*/) const;

    /** This creates all possible seeds with the passed central and bottom SP, using all top SP 
       * candidates which are stored in the data.CmSp member.  Seeds are scored by a quality score 
       * seeded by abs(d0), and modified if there is a second-seed confirmation or in case of PPP/SSS 
       * topologies. Then, they are written out via the newOneSeed method. 
       * @param[in] SPb Bottom Space point for the seed creation
       * @param[in] SP0 Central Space point for the seed creation
       * @param[in] Zob z0 estimate 
    **/ 
    void newOneSeedWithCurvaturesComparisonSSS
      (EventData& data, SiSpacePointForSeedITK*& SPb, SiSpacePointForSeedITK*& SP0, float Zob) const;
    void newOneSeedWithCurvaturesComparisonPPP
      (EventData& data, SiSpacePointForSeedITK*& SPb, SiSpacePointForSeedITK*& SP0, float Zob) const;
    void newOneSeedWithCurvaturesComparisonSeedConfirmation
      (EventData& data, SiSpacePointForSeedITK*& SPb, SiSpacePointForSeedITK*& SP0, float Zob) const;

 
     /** Helper method to determine if a seed 
       * is 'confirmed' - this means that a second
       * seed exists with compatible curvature, 
       * the same bottom and central SP, but a 
       * different third SP. This information 
       * is stored in a modification of the 
       * seed quality, which we check here. 
       * @param[in] bottomSP: bottom space point
       * @param[in] topSP: top space point 
       * @param[in] quality: seed quality
       * @return true if the seed is confirmed, false otherwise 
    **/ 
    bool isConfirmedSeed(const InDet::SiSpacePointForSeedITK* bottomSP, const InDet::SiSpacePointForSeedITK* topSP, float quality) const;


    void sort(std::vector<FloatInt>& s, int start, int size) const;
    bool newVertices(EventData& data, const std::list<Trk::Vertex>&) const;
    void findNext(EventData& data) const;
    bool isZCompatible(EventData& data, float&,float&,float&) const;
    void convertToBeamFrameWork(EventData& data, const Trk::SpacePoint*const&,float*) const;
    bool isUsed(const Trk::SpacePoint*, const Trk::PRDtoTrackMap &prd_to_track_map) const;

    void initializeEventData(EventData& data) const;
  };

} // end of name space

///////////////////////////////////////////////////////////////////
// Object-function for curvature seeds comparison
///////////////////////////////////////////////////////////////////

class comCurvatureITK {
public:
  bool operator ()
  (const std::pair<float,InDet::SiSpacePointForSeedITK*>& i1, 
   const std::pair<float,InDet::SiSpacePointForSeedITK*>& i2)
  {
    return i1.first < i2.first;
  }
};


///////////////////////////////////////////////////////////////////
// Test is space point used
///////////////////////////////////////////////////////////////////

namespace InDet {
  inline
  bool SiSpacePointsSeedMaker_ITK::isUsed(const Trk::SpacePoint* sp, const Trk::PRDtoTrackMap &prd_to_track_map) const
  {
    const Trk::PrepRawData* d = sp->clusterList().first;
    if (!d || !prd_to_track_map.isUsed(*d)) return false;

    d = sp->clusterList().second;
    if (!d || prd_to_track_map.isUsed(*d)) return true;

    return false;
  }

///////////////////////////////////////////////////////////////////
  // The procedure sorts the elements into ascending order.
  ///////////////////////////////////////////////////////////////////
  
  inline
  void SiSpacePointsSeedMaker_ITK::sort(std::vector<FloatInt>& s, int start, int size) const
  {
    //QuickSort for fast tracking currently buggy
    //TBC if really faster than std::sort
    //Using std::sort in all cases for now
    std::sort(s.begin()+start,s.begin()+start+size,[](const FloatInt a,const FloatInt b)->bool {return a.Fl < b.Fl;});
  }

}

#endif // SiSpacePointsSeedMaker_ITK_H
