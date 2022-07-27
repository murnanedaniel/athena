/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


/** @file InDetGlobalTrackMonTool.h
 * Implementation of inner detector global track monitoring tool
 * based on Segment tool
 *
 *@author
 * Heidi Sandaker <Heidi.Sandaker@cern.ch> @n

 * Arshak Tonoyan <Arshak.Tonyoan@cern.ch> @n
 * Thomas Burgess <Thomas.Burgess@cern.ch>
 * Gaetano Barone <Gaetano.Barone@cern.ch> @n 
 * 
 * $Id: InDetGlobalTrackMonTool.h,v 1.24 2009-04-06 17:18:19 kastanas Exp $
 *
 ****************************************************************************/

#ifndef InDetGlobalTrackMonTool_H
#define InDetGlobalTrackMonTool_H

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "TrkToolInterfaces/IResidualPullCalculator.h"
#include "TrkToolInterfaces/ITrackHoleSearchTool.h"
#include "TrkVertexFitterInterfaces/ITrackToVertexIPEstimator.h"
#include "TrkToolInterfaces/IUpdator.h"

#include "PixelGeoModel/IBLParameterSvc.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/VertexContainer.h"

//Detector Managers
#include "AtlasDetDescr/AtlasDetectorID.h"
#include "InDetIdentifier/TRT_ID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetIdentifier/PixelID.h"

//Local
//Framework
#include "StoreGate/ReadHandleKey.h"

//Standard c++
#include <string>
#include <map>
#include <vector>
#include "GaudiKernel/ToolHandle.h"
#include "AthenaBaseComps/AthAlgorithm.h"

//Predeclarations
class IInterface;
class StatusCode;
class TH1I;
class TH1F;
class TH2F;
class TH1I_LW;
class TH1F_LW;
class TH2F_LW;
class TProfile_LW;
class TProfile;
class TProfile2D;

/// Monitoring tool derived from InDetGlobalMotherMonTool
/// Contains the global track information for the ID

class InDetGlobalTrackMonTool : public ManagedMonitorToolBase//public InDetGlobalMotherTrigMonTool
{
public:
    ///Default constructor
    InDetGlobalTrackMonTool( 
	const std::string & type,
	const std::string & name,
	const IInterface* parent);
    
    ///Virtual destructor
    virtual ~InDetGlobalTrackMonTool() {}
    
    virtual StatusCode initialize();

    ///@name histos Book, fill and proc histograms
    ///@{

    ///@copydoc InDetGlobalMotherMonTool::bookHistograms()
    virtual StatusCode bookHistograms();
    virtual StatusCode bookHistogramsRecurrent();

    ///@copydoc InDetGlobalMotherMonTool::fillHistograms()
    virtual StatusCode fillHistograms();

    ///@copydoc InDetGlobalMotherMonTool::procHistograms()
    virtual StatusCode procHistograms();

    /// Functions to fill individual sets of histograms
    void FillForwardTracks( const xAOD::TrackParticle *trackPart);
    void FillEtaPhi( const xAOD::TrackParticle *trackPart );
    void FillHits( const xAOD::TrackParticle *trackPart );
    void FillTIDE();
    void FillHoles( const xAOD::TrackParticle *trackPart );
    void FillHitMaps( const xAOD::TrackParticle *trackPart );
    void FillHoleMaps( const xAOD::TrackParticle *trackPart );
    ///@} 
 
private:
    ToolHandle< InDet::IInDetTrackSelectionTool > m_baseline_selTool;
    ToolHandle< InDet::IInDetTrackSelectionTool > m_tight_selTool;

    /// Switch for hole searching
    bool m_doHolePlots;
    bool m_DoHoles_Search;
    /// Switch for hitmaps
    bool m_doHitMaps;

    bool m_doTide;
    bool m_doTideResiduals;
    bool m_doForwardTracks;
    bool m_doIBL;

    unsigned int m_nBinsEta;
    unsigned int m_nBinsPhi;
    unsigned int m_trackBin;
    unsigned int m_trackMax{};

    /// Contants for various histogram properties
    const float m_c_etaRange;
    const float m_c_etaTrackletsMin;
    const float m_c_etaTrackletsMax;
    const float m_c_etaRangeTRT;

    const float m_c_range_LB;

    const std::array<std::string,4> m_c_detector_labels;

    ServiceHandle <IBLParameterSvc> m_IBLParameterSvc;
    ToolHandle <Trk::ITrackHoleSearchTool> m_holes_search_tool;
    ToolHandle<Trk::IResidualPullCalculator> m_residualPullCalculator;
    PublicToolHandle< Trk::ITrackToVertexIPEstimator >  m_trackToVertexIPEstimator
       {this,"TrackToVertexIPEstimator","Trk::TrackToVertexIPEstimator",""};
    ToolHandle<Trk::IUpdator>             m_iUpdator;
    
    SG::ReadHandleKey<xAOD::TrackParticleContainer> m_TrackParticleName{this,"TrackParticleContainerName","InDetTrackParticles","TrackParticle Collection for Global Monitoring"};
    SG::ReadHandleKey<xAOD::TrackParticleContainer> m_ForwardTrackParticleName{this,"ForwardTrackParticleContainerName","InDetForwardTrackParticles","Forward TrackParticle Collection for Global Monitoring"};

    SG::ReadHandleKey<xAOD::JetContainer> m_JetsName{this,"JetCollection","AntiKt4EMTopoJets","Jet Collection for Global Track Monitoring"};
    SG::ReadHandleKey<xAOD::VertexContainer> m_vertexKey { this, "VertexContainer", "PrimaryVertices", "primary vertex container" };

    //--- Shift histograms ----------------------------------------
    
    // Gaetano Holes : 
    // SCT Holes : 
    TH1F * m_sct_holes;
    // TRT Holes : 
    TH1F * m_trt_holes;
    // Pixel Holes : 
    TH1F * m_pixel_holes;
    // Combined Holes : 
    TH1F * m_comb_holes;
    // Silicon  vs TRT holes 
    TH2F  *m_silicon_vs_trt;
    TH2F  *m_sct_vs_pixels;
    
    // Combined number of holes vs track quality
    TH2F *m_holes_quality;
    TProfile *m_holes_quality_profile;
    
    TH1I * m_Trk_Base{};
    ///Distribution of eta vs phi for combined tracks
    TH2F  * m_Trk_eta_phi_Base;
    TH2F  * m_Trk_eta_phi_Tight;
    TProfile2D  * m_Trk_eta_phi_Tight_ratio;
    TProfile2D  * m_Trk_eta_phi_noIBLhit_ratio;
    TProfile2D  * m_Trk_eta_phi_noBLhit_ratio;
    TProfile2D  * m_Trk_eta_phi_noTRText_ratio;
    
    // Total number of tracks per LB for each subdetector 
    TProfile * m_Trk_nBase_LB;
    TProfile * m_Trk_nTight_LB;

    TProfile * m_Trk_noIBLhits_LB;
    TProfile * m_Trk_noBLhits_LB;
    TProfile * m_Trk_noTRText_LB;

    TProfile * m_Trk_noIBLhits_frac_LB;
    TProfile * m_Trk_noBLhits_frac_LB;
    TProfile * m_Trk_noTRText_frac_LB;

    ///@name Detector managers
    ///{@

    /// the TRT ID helper
    const TRT_ID *m_trtID;
    
    /// the SCT ID helper 
    const SCT_ID *m_sctID;  

    /// the Pixel ID helper 
    const PixelID *m_pixelID;
    
    ///@}

    // changes Holes mapping 
    TH2F     *m_holes_eta_phi;
    TProfile2D  *m_holes_eta_pt;
    TProfile2D  *m_holes_phi_pt;
    TProfile2D  *m_holes_eta_phi_n;
    TProfile *m_holes_hits;
    TH2F     *m_holesvshits;
    TH2F     *m_holesvshits_ECA;
    TH2F     *m_holesvshits_ECC;
    TH2F     *m_holesvshits_BA;

    TH2F     *m_ID_hitmap_x_y;
    TH2F     *m_ID_hitmap_x_y_eca;
    TH2F     *m_ID_hitmap_x_y_ecc;
    TH2F     *m_HolesMAP_XY;
    TH2F     *m_HolesMAP_ZX;
    TH2F     *m_HolesMAP_ZR;


    std::array<TProfile2D *, 4> m_trk_hits_eta_phi;
    std::array<TProfile2D *, 4> m_trk_disabled_eta_phi;
    std::array<TProfile *,4> m_trk_hits_LB;
    
    TProfile2D * m_trk_shared_pix_eta_phi;
    TProfile2D * m_trk_split_pix_eta_phi;
    TProfile2D * m_trk_shared_sct_eta_phi;

    TProfile2D * m_trk_holes_pix_eta_phi;
    TProfile2D * m_trk_holes_sct_eta_phi;

    TProfile * m_trk_jetassoc_d0_reso_dr{}; 
    TProfile * m_trk_jetassoc_z0_reso_dr{}; 
    TProfile * m_trk_jetassoc_split_pix_dr{}; 
    TProfile * m_trk_jetassoc_shared_pix_dr{}; 

    TProfile * m_trk_jetassoc_res_pix_l0_x_dr{};
    TProfile * m_trk_jetassoc_res_pix_l1_x_dr{};
    TProfile * m_trk_jetassoc_res_pix_l2_x_dr{};
    TProfile * m_trk_jetassoc_res_pix_l3_x_dr{};

    TProfile * m_trk_jetassoc_res_pix_l0_y_dr{};
    TProfile * m_trk_jetassoc_res_pix_l1_y_dr{};
    TProfile * m_trk_jetassoc_res_pix_l2_y_dr{};
    TProfile * m_trk_jetassoc_res_pix_l3_y_dr{};

    TProfile * m_trk_jetassoc_res_pix_eca_x_dr{};
    TProfile * m_trk_jetassoc_res_pix_eca_y_dr{};

    TProfile * m_trk_jetassoc_res_pix_ecc_x_dr{};
    TProfile * m_trk_jetassoc_res_pix_ecc_y_dr{};

    TProfile * m_trk_jetassoc_ip_reso_lb{}; 
    TProfile * m_trk_jetassoc_split_pix_lb{}; 
    TProfile * m_trk_jetassoc_shared_pix_lb{}; 
   
    //--- Combined tracks debug histograms-----------------------------------
    TH2F * m_Trk_FORW_FA_eta_phi;
    TH2F * m_Trk_FORW_FC_eta_phi;

    TH1F * m_Trk_FORW_qoverp;
    TH1F * m_Trk_FORW_chi2;

    ///Number of PIX hits per track
    TH1I * m_Trk_FORW_FA_nPIXhits;
    TH1I * m_Trk_FORW_FC_nPIXhits;
    
    ///@}

    template <class histClass>
    StatusCode registerManHist ( 
	histClass * & target, const std::string & path, Interval_t interval,
	const std::string & name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	const std::string & xlabel, const std::string & ylabel = "" )
	{
	    target = new histClass( name.c_str(), title.c_str(), 
				    nbinsx, xlow, xhi );
	    target->GetXaxis()->SetTitle( xlabel.c_str() );
	    target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return regHist( target, path, interval );
	}

    template <class histClass>
    StatusCode registerManHist ( 
	histClass * & target, const std::string & path, Interval_t interval,
	const std::string & name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	int nbinsy, double ylow, double yhi,
	const std::string & xlabel, const std::string & ylabel = "" )
	{
	    target = new histClass( name.c_str(), title.c_str(), 
				    nbinsx, xlow, xhi,
				    nbinsy, ylow, yhi );
	    target->GetXaxis()->SetTitle( xlabel.c_str() );
	    target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return regHist( target, path, interval );
	}
    
    template <class histClass>
    StatusCode registerHistI ( 
	MonGroup & theGroup, histClass * & target, const std::string & name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	const std::string & xlabel = "", const std::string & ylabel = "" )
	{
	    target = histClass::create( name.c_str(), title.c_str(), nbinsx, xlow, xhi );
	    
	    if ( xlabel != "" )
		target->GetXaxis()->SetTitle( xlabel.c_str() );

	    if ( ylabel != "" )
		target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return theGroup.regHist( target );
	}

    template <class histClass>
    StatusCode registerHistI ( 
	MonGroup & theGroup, histClass * & target, const std::string & name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	int nbinsy, double ylow, double yhi,
	const std::string & xlabel = "", const std::string & ylabel = "" )
	{
	    target = histClass::create( name.c_str(), title.c_str(), 
					nbinsx, xlow, xhi,
					nbinsy, ylow, yhi );
	    
	    if ( xlabel != "" )
		target->GetXaxis()->SetTitle( xlabel.c_str() );
	    
	    if ( ylabel != "" )
		target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return theGroup.regHist( target );
	}
    
    template <class histClass>
    StatusCode registerHistIR ( 
	MonGroup & theGroup, histClass * & target, const std::string & name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	const std::string & xlabel = "", const std::string & ylabel = "" )
	{
	    target = new histClass( name.c_str(), title.c_str(), 
				    nbinsx, xlow, xhi );
	    
	    if ( xlabel != "" )
		target->GetXaxis()->SetTitle( xlabel.c_str() );

	    if ( ylabel != "" )
		target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return theGroup.regHist( target );
	}
    
    template <class histClass>
    StatusCode registerHistIR ( 
	MonGroup & theGroup, histClass * & target, const std::string& name, const std::string & title,
	int nbinsx, double xlow, double xhi,
	int nbinsy, double ylow, double yhi,
	const std::string & xlabel = "", const std::string & ylabel = "" )
	{
	    target = new histClass( name.c_str(), title.c_str(), 
				    nbinsx, xlow, xhi,
				    nbinsy, ylow, yhi );
	    
	    if ( xlabel != "" )
		target->GetXaxis()->SetTitle( xlabel.c_str() );
	    
	    if ( ylabel != "" )
		target->GetYaxis()->SetTitle( ylabel.c_str() );
	    
	    return theGroup.regHist( target );
	}
    


};

#endif
