/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetTrackSystematicsTools/JetTrackFilterTool.h"
#include "InDetTrackSystematicsTools/InDetTrackTruthOriginTool.h"
#include "InDetTrackSystematicsTools/InDetTrackTruthOriginDefs.h"
#include "xAODTracking/TrackParticleContainer.h"

#include "FourMomUtils/xAODP4Helpers.h"
#include "PathResolver/PathResolver.h"

#include <TH2.h>
#include <TRandom3.h>
#include <TFile.h>
#include <utility>

namespace InDet {

  static const CP::SystematicSet FilterSystematics = 
    {
      InDet::TrackSystematicMap.at(TRK_EFF_LOOSE_TIDE),
      InDet::TrackSystematicMap.at(TRK_FAKE_RATE_TIGHT_TIDE),
      InDet::TrackSystematicMap.at(TRK_FAKE_RATE_LOOSE_TIDE)
    };

  JetTrackFilterTool::JetTrackFilterTool(const std::string& name) :
    InDetTrackSystematicsTool(name),
    m_trackOriginTool("InDet::InDetTrackTruthOriginTool")
  {

#ifndef XAOD_STANDALONE
    declareInterface<IJetTrackFilterTool>(this);
#endif

    declareProperty("Seed", m_seed, "Seed used to initialize the RNG");
    declareProperty("DeltaR", m_deltaR, "Delta-R cut in which to apply jet-track efficiency rejection");
    declareProperty("trkEffSystScale", m_trkEffSystScale, "Option to scale the effect of the systematic (default 1)");
    declareProperty("FakeUncertainty",  m_fakeUncertTIDE, "Option to set the fake uncertainty");
    declareProperty("calibFileNomEff", m_calibFileNomEff = "InDetTrackSystematicsTools/CalibData_21.2_2018-v15/TrackingRecommendations_final_rel21.root");
    declareProperty("calibFileJetEff", m_calibFileJetEff = "InDetTrackSystematicsTools/CalibData_21.2_2018-v15/TIDErejectProbv3.root");
    declareProperty("trackOriginTool", m_trackOriginTool);
  }

  StatusCode JetTrackFilterTool::initialize()
  {

    m_rnd = std::make_unique<TRandom3>(m_seed);

    ATH_CHECK ( initTIDEEffSystHistogram( m_trkEffSystScale,
               m_effForJetPt,
               m_calibFileJetEff,
               "h1") );
    ATH_CHECK( initObject<TH2>( m_trkNomEff,
               m_calibFileNomEff, 
               "EfficiencyVSEtaPt_AfterRebinningNominal_Loose" ) );

    ATH_MSG_INFO( "Using for nominal track efficiency the calibration file " << PathResolverFindCalibFile(m_calibFileNomEff) );
    ATH_MSG_INFO( "Using for jet track efficiency the calibration file " << PathResolverFindCalibFile(m_calibFileJetEff) );

    ATH_CHECK ( m_trackOriginTool.retrieve() );

    ATH_CHECK ( InDetTrackSystematicsTool::initialize() );

    return StatusCode::SUCCESS;
  }


  JetTrackFilterTool::~JetTrackFilterTool()
  {
    delete m_effForJetPt; m_effForJetPt = nullptr;
    delete m_trkNomEff; m_trkNomEff = nullptr;
  }

  bool JetTrackFilterTool::accept(const xAOD::TrackParticle* track, const xAOD::Jet* jet) const
  {

    if ( track == nullptr) {
      ATH_MSG_DEBUG( "Pointer to track is null!" );
      return false;
    }
    if ( jet == nullptr ) {
      ATH_MSG_DEBUG( "Pointer to jet is null." );
      return true;
    }

    // if the track is outside of the range of this jet, then it is allowed to pass
    constexpr bool useRapidity = false; // use eta instead of rapidity - the default for this function is true
    if ( !xAOD::P4Helpers::isInDeltaR( *track, *jet, m_deltaR, useRapidity ) ) return true;

    if ( isActive( TRK_EFF_LOOSE_TIDE ) ) {
      // the probability to drop a track scales with the tracking efficiency,
      // on account of the method used to derive the uncertainties
      float probDrop = std::fabs(m_trkEffSystScale); // default is one; adjust this parameter to increase / decrease the effect
      probDrop *= getFractionDropped( 1.0, m_effForJetPt, jet );
      probDrop *= getNomTrkEff( track );
      if ( m_rnd->Uniform(0, 1) < probDrop ) return false;
    }

    int origin = m_trackOriginTool->getTrackOrigin(track);

    if( isActive( TRK_FAKE_RATE_LOOSE_TIDE ) ){
      if ( InDet::TrkOrigin::isFake(origin) ) {
        if(m_rnd->Uniform(0, 1) <  m_fakeUncertTIDE) return false;
      }
    }

    if( isActive( TRK_FAKE_RATE_TIGHT_TIDE ) ){
      if ( InDet::TrkOrigin::isFake(origin) ) {
        ATH_MSG_DEBUG("Track fakes in jets uncertainty (Tight) covered by inclusive (Tight) uncertainty - operating in pass-through mode...");
        return true;
      }
    }

    return true;
  }

  bool JetTrackFilterTool::accept( const xAOD::TrackParticle* track, const xAOD::JetContainer* jets ) const
  {
    if ( jets == nullptr ) {
      ATH_MSG_DEBUG( "Pointer to jet container is null." );
      return true;
    }
    // check that the track passes every jet
    for ( const auto* jet : *jets ) {
      if ( !accept( track, jet ) ) return false;
    }
    return true;
  }

  StatusCode JetTrackFilterTool::initTIDEEffSystHistogram(float scale, TH1 *&histogram, std::string rootFileName, std::string histogramName) const
  {
    ATH_CHECK( initObject<TH1>(histogram, rootFileName, histogramName) );

    for(int binx=1; binx<=histogram->GetNbinsX(); binx++) {

      float content = histogram->GetBinContent(binx); // get the bin content
      content *= scale; // scale systematic
      if(content > 1) content = 1; // protection: no larger than 100% uncertainty

      histogram->SetBinContent(binx, content);
    }

    return StatusCode::SUCCESS;
  }

  float JetTrackFilterTool::getFractionDropped(float fDefault, TH1 *histogram, const xAOD::Jet* jet) const
  {

    if( histogram == nullptr ) {
      return fDefault;
    }

    auto pt = jet->pt()*1.e-3;
    // get the appropriate ptBin, but also use the highest valid bin in case of overflow
    int ptBin = std::min(histogram->FindBin(pt), histogram->GetNbinsX());

    //    ATH_MSG_VERBOSE( "pt / bin no: " << pt << " " << ptBin );
    float frac = histogram->GetBinContent(ptBin);
    if( frac > 1. ) {
      ATH_MSG_WARNING( "Fraction from histogram " << histogram->GetName()
           << " is greater than 1. Setting fraction to 1." );
      frac = 1.;
    }
    if( frac < 0. ) {
      ATH_MSG_WARNING( "Fraction from histogram " << histogram->GetName()
           << " is less than 0. Setting fraction to 0." );
      frac = 0.;
    }
    //    ATH_MSG_DEBUG( "fraction dropped (from jet-pt histogram): " << frac );
    return frac;
  }

  float JetTrackFilterTool::getNomTrkEff(const xAOD::TrackParticle* track) const
  {
    if (m_trkNomEff == nullptr) {
      ATH_MSG_ERROR( "Nominal track efficiency histogram is not property initialized!" );
      return 0.;
    }
    // this histogram has pt on the x-axis and eta on the y-axis, unlike some other histograms used in this package
    // make sure to convert to GeV
    return m_trkNomEff->GetBinContent(std::as_const(m_trkNomEff)->FindBin(track->pt()*1e-3, track->eta()));
  }

  bool JetTrackFilterTool::isAffectedBySystematic( const CP::SystematicVariation& syst ) const
  {
    return InDetTrackSystematicsTool::isAffectedBySystematic( syst );
  }

  CP::SystematicSet JetTrackFilterTool::affectingSystematics() const
  {
    return FilterSystematics;
  }

  CP::SystematicSet JetTrackFilterTool::recommendedSystematics() const
  {
    return InDetTrackSystematicsTool::recommendedSystematics();
  }

  StatusCode JetTrackFilterTool::applySystematicVariation( const CP::SystematicSet& systs )
  {
    return InDetTrackSystematicsTool::applySystematicVariation(systs);
  }


} // namespace InDet
