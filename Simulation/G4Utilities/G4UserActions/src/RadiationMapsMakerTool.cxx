/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "G4UserActions/RadiationMapsMakerTool.h"

namespace G4UA{ 

  RadiationMapsMakerTool::RadiationMapsMakerTool(const std::string& type, 
						 const std::string& name,
						 const IInterface* parent)
    : UserActionToolBase<RadiationMapsMaker>(type, name, parent),
      m_radMapsFileName("RadMaps.root")  
  {
    /// Output Filename for the Radiation Maps
    declareProperty("RadMapsFileName", m_radMapsFileName);
    /// map granularities 
    /// number of bins in r and z for all 2D maps
    declareProperty("NBinsR"    , m_config.nBinsr);
    declareProperty("NBinsZ"    , m_config.nBinsz);
    /// number of bins in r, z and phi for all 3D maps
    declareProperty("NBinsR3D"  , m_config.nBinsr3d);
    declareProperty("NBinsZ3D"  , m_config.nBinsz3d);
    declareProperty("NBinsPhi3D", m_config.nBinsphi3d);
    /// map ranges
    /// for Zoomed area in 2D and 3D
    declareProperty("RMinZoom"  , m_config.rMinZoom);
    declareProperty("RMaxZoom"  , m_config.rMaxZoom);
    declareProperty("ZMinZoom"  , m_config.zMinZoom);
    declareProperty("ZMaxZoom"  , m_config.zMaxZoom);
    /// for Full detector in 2D
    declareProperty("RMinFull"  , m_config.rMinFull);
    declareProperty("RMaxFull"  , m_config.rMaxFull);
    declareProperty("ZMinFull"  , m_config.zMinFull);
    declareProperty("ZMaxFull"  , m_config.zMaxFull);
    /// for Zoomed area in 3D 
    declareProperty("PhiMinZoom", m_config.phiMinZoom);
    declareProperty("PhiMaxZoom", m_config.phiMaxZoom);
  }

  //---------------------------------------------------------------------------
  // Initialize Configurable Properties
  //---------------------------------------------------------------------------
  StatusCode RadiationMapsMakerTool::initialize()
  {
    ATH_MSG_INFO( "Initializing  " << name() << "\n" <<
                  "OutputFile:   " << m_radMapsFileName   << "\n"                << 
                  "2D Maps:      " << m_config.nBinsz     << " |z|-bins, "       << 
                                      m_config.nBinsr     << " r-bins"           << "\n"                << 
                  "Zoom:         " << m_config.zMinZoom   << " < |z|/cm < "      << m_config.zMaxZoom   << ", " << 
                                      m_config.rMinZoom   << " < r/cm < "        << m_config.rMaxZoom   << "\n" << 
                  "Full:         " << m_config.zMinFull   << " < |z|/cm < "      << m_config.zMaxFull   << ", " << 
                                      m_config.rMinFull   << " < r/cm < "        << m_config.rMaxFull   << "\n" << 
                  "3D Maps:      " << m_config.nBinsz3d   << " |z|-bins, "       << 
                                      m_config.nBinsr3d   << " r-bins, "         << 
                                      m_config.nBinsphi3d << " phi-bins"         << "\n"                << 
                  "Zoom:         " << m_config.zMinZoom   << " < |z|/cm < "      << m_config.zMaxZoom   << ", " << 
                                      m_config.rMinZoom   << " < r/cm < "        << m_config.rMaxZoom   << ", " <<
                                      m_config.phiMinZoom << " < phi/degrees < " << m_config.phiMaxZoom );
      
    return StatusCode::SUCCESS;
  }

  //---------------------------------------------------------------------------
  // Merge results from all threads
  //---------------------------------------------------------------------------
  StatusCode RadiationMapsMakerTool::finalize()
  {
    ATH_MSG_DEBUG( "Finalizing " << name() );

    // first make sure the vectors are empty

    RadiationMapsMaker::Report maps;

    maps.m_rz_tid .resize(0);
    maps.m_rz_eion.resize(0);
    maps.m_rz_niel.resize(0);
    maps.m_rz_h20 .resize(0);
    
    maps.m_full_rz_tid .resize(0);
    maps.m_full_rz_eion.resize(0);
    maps.m_full_rz_niel.resize(0);
    maps.m_full_rz_h20 .resize(0);
    
    maps.m_3d_tid .resize(0);
    maps.m_3d_eion.resize(0);
    maps.m_3d_niel.resize(0);
    maps.m_3d_h20 .resize(0);

    // then resize to proper size and initialize with 0's 

    maps.m_rz_tid .resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_rz_eion.resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_rz_niel.resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_rz_h20 .resize(m_config.nBinsz*m_config.nBinsr,0.0);
    
    maps.m_full_rz_tid .resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_full_rz_eion.resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_full_rz_niel.resize(m_config.nBinsz*m_config.nBinsr,0.0);
    maps.m_full_rz_h20 .resize(m_config.nBinsz*m_config.nBinsr,0.0);
    
    maps.m_3d_tid .resize(m_config.nBinsz3d*m_config.nBinsr3d*m_config.nBinsphi3d,0.0);
    maps.m_3d_eion.resize(m_config.nBinsz3d*m_config.nBinsr3d*m_config.nBinsphi3d,0.0);
    maps.m_3d_niel.resize(m_config.nBinsz3d*m_config.nBinsr3d*m_config.nBinsphi3d,0.0);
    maps.m_3d_h20 .resize(m_config.nBinsz3d*m_config.nBinsr3d*m_config.nBinsphi3d,0.0);

    // merge radiation map vectors from threads
    // Accumulate the results across threads
    m_actions.accumulate(maps, &RadiationMapsMaker::getReport,
                         &RadiationMapsMaker::Report::merge);

    TFile * f = new TFile(m_radMapsFileName.c_str(),"RECREATE");

    TH2D * h_rz_tid  = new TH2D("rz_tid" ,"rz_tid" ,m_config.nBinsz,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr,m_config.rMinZoom,m_config.rMaxZoom);
    TH2D * h_rz_eion = new TH2D("rz_eion","rz_eion",m_config.nBinsz,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr,m_config.rMinZoom,m_config.rMaxZoom);
    TH2D * h_rz_niel = new TH2D("rz_niel","rz_niel",m_config.nBinsz,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr,m_config.rMinZoom,m_config.rMaxZoom);
    TH2D * h_rz_h20  = new TH2D("rz_h20" ,"rz_h20" ,m_config.nBinsz,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr,m_config.rMinZoom,m_config.rMaxZoom);

    h_rz_tid  ->SetXTitle("|z| [cm]");
    h_rz_eion ->SetXTitle("|z| [cm]");
    h_rz_niel ->SetXTitle("|z| [cm]");
    h_rz_h20  ->SetXTitle("|z| [cm]");

    h_rz_tid  ->SetYTitle("r [cm]");
    h_rz_eion ->SetYTitle("r [cm]");
    h_rz_niel ->SetYTitle("r [cm]");
    h_rz_h20  ->SetYTitle("r [cm]");

    h_rz_tid  ->SetZTitle("TID [Gy]");
    h_rz_eion ->SetZTitle("E_{ion}/V [MeV/cm^{3}]");
    h_rz_niel ->SetZTitle("NIEL [n_{eq}/cm^{2}]");
    h_rz_h20  ->SetZTitle("SEE [h_{>20 MeV}/cm^{2}]");

    TH2D *h_full_rz_tid  = new TH2D("full_rz_tid" ,"full_rz_tid" ,m_config.nBinsz,m_config.zMinFull,m_config.zMaxFull,m_config.nBinsr,m_config.rMinFull,m_config.rMaxFull);
    TH2D *h_full_rz_eion = new TH2D("full_rz_eion","full_rz_eion",m_config.nBinsz,m_config.zMinFull,m_config.zMaxFull,m_config.nBinsr,m_config.rMinFull,m_config.rMaxFull);
    TH2D *h_full_rz_niel = new TH2D("full_rz_niel","full_rz_niel",m_config.nBinsz,m_config.zMinFull,m_config.zMaxFull,m_config.nBinsr,m_config.rMinFull,m_config.rMaxFull);
    TH2D *h_full_rz_h20  = new TH2D("full_rz_h20" ,"full_rz_h20" ,m_config.nBinsz,m_config.zMinFull,m_config.zMaxFull,m_config.nBinsr,m_config.rMinFull,m_config.rMaxFull);

    h_full_rz_tid  ->SetXTitle("|z| [cm]");
    h_full_rz_eion ->SetXTitle("|z| [cm]");
    h_full_rz_niel ->SetXTitle("|z| [cm]");
    h_full_rz_h20  ->SetXTitle("|z| [cm]");

    h_full_rz_tid  ->SetYTitle("r [cm]");
    h_full_rz_eion ->SetYTitle("r [cm]");
    h_full_rz_niel ->SetYTitle("r [cm]");
    h_full_rz_h20  ->SetYTitle("r [cm]");

    h_full_rz_tid  ->SetZTitle("TID [Gy]");
    h_full_rz_eion ->SetZTitle("E_{ion}/V [MeV/cm^{3}]");
    h_full_rz_niel ->SetZTitle("NIEL [n_{eq}/cm^{2}]");
    h_full_rz_h20  ->SetZTitle("SEE [h_{>20 MeV}/cm^{2}]");

    TH3D * h_3d_tid  = new TH3D("h3d_tid" ,"h3d_tid" ,m_config.nBinsz3d,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr3d,m_config.rMinZoom,m_config.rMaxZoom,m_config.nBinsphi3d,m_config.phiMinZoom,m_config.phiMaxZoom);
    TH3D * h_3d_eion = new TH3D("h3d_eion","h3d_eion",m_config.nBinsz3d,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr3d,m_config.rMinZoom,m_config.rMaxZoom,m_config.nBinsphi3d,m_config.phiMinZoom,m_config.phiMaxZoom);
    TH3D * h_3d_niel = new TH3D("h3d_niel","h3d_niel",m_config.nBinsz3d,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr3d,m_config.rMinZoom,m_config.rMaxZoom,m_config.nBinsphi3d,m_config.phiMinZoom,m_config.phiMaxZoom);
    TH3D * h_3d_h20  = new TH3D("h3d_h20" ,"h3d_h20" ,m_config.nBinsz3d,m_config.zMinZoom,m_config.zMaxZoom,m_config.nBinsr3d,m_config.rMinZoom,m_config.rMaxZoom,m_config.nBinsphi3d,m_config.phiMinZoom,m_config.phiMaxZoom);

    h_3d_tid  ->SetXTitle("|z| [cm]");
    h_3d_eion ->SetXTitle("|z| [cm]");
    h_3d_niel ->SetXTitle("|z| [cm]");
    h_3d_h20  ->SetXTitle("|z| [cm]");

    h_3d_tid  ->SetYTitle("r [cm]");
    h_3d_eion ->SetYTitle("r [cm]");
    h_3d_niel ->SetYTitle("r [cm]");
    h_3d_h20  ->SetYTitle("r [cm]");

    h_3d_tid  ->SetZTitle("#phi [#circ]");
    h_3d_eion ->SetZTitle("#phi [#circ]");
    h_3d_niel ->SetZTitle("#phi [#circ]");
    h_3d_h20  ->SetZTitle("#phi [#circ]");

    h_3d_tid  ->SetTitle("TID [Gy]");
    h_3d_eion ->SetTitle("E_{ion}/V [MeV/cm^{3}]");
    h_3d_niel ->SetTitle("NIEL [n_{eq}/cm^{2}]");
    h_3d_h20  ->SetTitle("SEE [h_{>20 MeV}/cm^{2}]");


    // normalize to volume element per bin
    for(int i=0;i<h_rz_tid->GetNbinsX();i++) { 
      for(int j=0;j<h_rz_tid->GetNbinsY();j++) { 
	int iBin = h_rz_tid->GetBin(i+1,j+1); 
	int vBin = m_config.nBinsr*i+j;
	double r0=h_rz_tid->GetYaxis()->GetBinLowEdge(j+1);
	double r1=h_rz_tid->GetYaxis()->GetBinUpEdge(j+1);
	double z0=h_rz_tid->GetXaxis()->GetBinLowEdge(i+1);
	double z1=h_rz_tid->GetXaxis()->GetBinUpEdge(i+1); 
	double vol=2*(z1-z0)*M_PI*(r1*r1-r0*r0); 
	double val;
	// TID
	val =maps.m_rz_tid[vBin];
	h_rz_tid->SetBinContent(iBin,val/vol);
	// EION
	val =maps.m_rz_eion[vBin];
	h_rz_eion->SetBinContent(iBin,val/vol);
	// NIEL
	val =maps.m_rz_niel[vBin];
	h_rz_niel->SetBinContent(iBin,val/vol);
	// SEE
	val =maps.m_rz_h20[vBin];
	h_rz_h20->SetBinContent(iBin,val/vol);
      }
    }
    h_rz_tid->Write();
    h_rz_eion->Write();
    h_rz_niel->Write();
    h_rz_h20->Write();

    // normalize to volume element per bin
    for(int i=0;i<h_full_rz_tid->GetNbinsX();i++) { 
      for(int j=0;j<h_full_rz_tid->GetNbinsY();j++) { 
	int iBin = h_full_rz_tid->GetBin(i+1,j+1); 
	int vBin = m_config.nBinsr*i+j;
	double r0=h_full_rz_tid->GetYaxis()->GetBinLowEdge(j+1);
	double r1=h_full_rz_tid->GetYaxis()->GetBinUpEdge(j+1);
	double z0=h_full_rz_tid->GetXaxis()->GetBinLowEdge(i+1);
	double z1=h_full_rz_tid->GetXaxis()->GetBinUpEdge(i+1); 
	double vol=2*(z1-z0)*M_PI*(r1*r1-r0*r0); 
	double val;
	// TID
	val =maps.m_full_rz_tid[vBin];
	h_full_rz_tid->SetBinContent(iBin,val/vol);
	// EION
	val =maps.m_full_rz_eion[vBin];
	h_full_rz_eion->SetBinContent(iBin,val/vol);
	// NIEL
	val =maps.m_full_rz_niel[vBin];
	h_full_rz_niel->SetBinContent(iBin,val/vol);
	// SEE
	val =maps.m_full_rz_h20[vBin];
	h_full_rz_h20->SetBinContent(iBin,val/vol);
      }
    }
    h_full_rz_tid->Write();
    h_full_rz_eion->Write();
    h_full_rz_niel->Write();
    h_full_rz_h20->Write();

    // normalize to volume element per bin
    for(int i=0;i<h_3d_tid->GetNbinsX();i++) { /* |z| */
      for(int j=0;j<h_3d_tid->GetNbinsY();j++) { /* r */
	for(int k=0;k<h_3d_tid->GetNbinsZ();k++) { /* phi */
	  int vBin = m_config.nBinsr3d*m_config.nBinsphi3d*i+m_config.nBinsphi3d*j+k;
	  int iBin = h_3d_tid->GetBin(i+1,j+1,k+1); 
	  double phi0=h_3d_tid->GetZaxis()->GetBinLowEdge(k+1);
	  double phi1=h_3d_tid->GetZaxis()->GetBinUpEdge(k+1);
	  double r0=h_3d_tid->GetYaxis()->GetBinLowEdge(j+1);
	  double r1=h_3d_tid->GetYaxis()->GetBinUpEdge(j+1);
	  double z0=h_3d_tid->GetXaxis()->GetBinLowEdge(i+1);
	  double z1=h_3d_tid->GetXaxis()->GetBinUpEdge(i+1); 
	  double vol=2*(z1-z0)*M_PI*(r1*r1-r0*r0)*(phi1-phi0)/360.; 
	  double val;
	  // TID
	  val =maps.m_3d_tid[vBin];
	  h_3d_tid->SetBinContent(iBin,val/vol);
	  // EION
	  val =maps.m_3d_eion[vBin];
	  h_3d_eion->SetBinContent(iBin,val/vol);
	  // NIEL
	  val =maps.m_3d_niel[vBin];
	  h_3d_niel->SetBinContent(iBin,val/vol);
	  // SEE
	  val =maps.m_3d_h20[vBin];
	  h_3d_h20->SetBinContent(iBin,val/vol);
	}
      }
    }
    h_3d_tid->Write();
    h_3d_eion->Write();
    h_3d_niel->Write();
    h_3d_h20->Write();

    f->Close();

    return StatusCode::SUCCESS;
  }

  std::unique_ptr<RadiationMapsMaker>  RadiationMapsMakerTool::makeAndFillAction(G4AtlasUserActions& actionList){
    ATH_MSG_INFO("Making a RadiationMapsMaker action");
    auto action = std::make_unique<RadiationMapsMaker>(m_config);
    actionList.runActions.push_back( action.get() );
    actionList.steppingActions.push_back( action.get() );
    return action;
  }

} // namespace G4UA 
