//
//   @file    comparitor.cxx         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2012 M.Sutton (sutt@cern.ch)    
//
//   $Id: comparitor.cxx, v0.0   Fri 12 Oct 2012 13:39:05 BST sutt $


// #include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h> 

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "TrigInDetAnalysis/Efficiency.h"

#include "ReadCards.h"

#include "Resplot.h"
#include "utils.h"
#include "label.h"
#include "DrawLabel.h" 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TColor.h"

#include "computils.h"

#include "AtlasStyle.h"
#include "AtlasLabels.h"


bool fulldbg = false;


int usage(const std::string& name, int status) { 
  std::ostream& s = std::cout;
  s << "Usage: " << name << " [OPTIONS]  test.root reference.root    chain1 chain2 chain2 ...\n\n"; 
  s << "  TIDA \'" << name << "\' plots comparison histograms\n\n"; 
  //  s << "Configuration: \n";
  //  s << "    -o filename   \tname of output grid (filename required)\n\n";
  s << "Options: \n";
  s << "    -c,  --config value       \t configure which histograms to plot from config file,\n";
  s << "    -t,  --tag value          \t appends tag 'value' to the end of output plot names, \n";
  s << "    -k,  --key value          \t prepends key 'value' to the from of output plots name, \n";
  s << "    -d,  --dir value          \t creates output files into directory, \"value\" \n";
  s << "    -e,  --efficiencies       \t make test efficiencies with respect to reference \n";
  s << "    -es, --effscale value     \t scale efficiencies to value\n";
  s << "    -er, --effscaleref  value \t scale efficiencies to value\n";
  s << "    -r,  --refit              \t refit all resplots\n";
  s << "    -l,  --labels             \t use specified labels for key\n";
  s << "         --taglabels          \t use specified additional labels \n";
  s << "    -nb  --nobayes            \t do not calculate Basyesian efficiency uncertaintiesr\n";
  s << "    -as, --atlasstyle         \t use ATLAS style\n";
  s << "    -al, --atlaslable value   \t set value for atlas label\n";
  s << "         --yrange  min max    \t use specified y axis range\n";  
  s << "    -ns, --nostats            \t do not show stats for mean and rms\n";
  s << "    -nm, --nomeans            \t do not show stats for the means\n";
  s << "    -nr, --noref              \t do not plot reference histograms\n";
  s << "         --normref            \t normalise the reference counting histograms to test histograms\n";
  s << "    -rc, --refchains values ..\t allow different reference chains for comparison\n";
  s << "    -s,  --swap pattern regex \t swap \"pattern\" in the reference chains name by \"regex\"\n";
  s << "    -np, --noplots            \t do not actually make any plot\n";
  s << "         --unscalepix         \t do not scale the number of pixels by 0.5 (scaled by default)\n";
  s << "         --chi2               \t show the chi2 with respect to the reference\n";
  s << "    -C,  --Cfiles             \t write C files also\n"; 
  s << "    -nw, --nowatermark        \t do not plot the release watermark\n"; 
  s << "         --nopng              \t do not print png files\n"; 
  s << "         --deleteref          \t delete unused reference histograms\n";
  s << "    -xo, --xoffset            \t relative x offset for the key\n"; 
  //  s << "         --fe,            \t relative x offset for the key\n"; 
  s << "    -h,  --help              \t this help\n";
  //  s << "\nSee " << PACKAGE_URL << " for more details\n"; 
  //  s << "\nReport bugs to <" << PACKAGE_BUGREPORT << ">";
  s << std::endl;
  return status;
}


void binwidth( TH1F* h ) { 
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    double w = h->GetBinLowEdge(i+1) - h->GetBinLowEdge(i);
    h->SetBinContent( i, h->GetBinContent(i)/w );
    h->SetBinError( i, h->GetBinError(i)/w );
  }
}

void ascale( TH1F* h, double s_ ) { 
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    h->SetBinContent( i, h->GetBinContent(i)*s_ );
    h->SetBinError( i, h->GetBinError(i)*s_ );
  }
}



// replace from a string
std::string fullreplace( std::string s, const std::string& s2, const std::string& s3) {
  if ( s2=="" || s2==s3 ) return s;  /// cowardly, don't replace string by itself 
  std::string::size_type pos;
  while ( (pos=s.find(s2)) != std::string::npos )  s.replace(pos, s2.size(), s3);
  return s;
} 



// std::string replace( std::string s, const std::string& pattern, const std::string& regex ) { 
//  if ( pattern!="" && regex!="" && s.find(pattern)!=std::string::npos ) s.replace( s.find(pattern), pattern.size(), regex );
//  return s;
// }


/// zero the contents of a 2d histogram 
void zero( TH2* h ) { 
  int Nx = h->GetNbinsX();
  int Ny = h->GetNbinsY();
  for ( int i=1 ; i<=Nx ; i++ ) { 
    for ( int j=1 ; j<=Ny ; j++ ) { 
      if ( h->GetBinContent( i, j )==0 ) { 
	//	h->SetBinContent( h->FindBin( i, j ), 1e-10 );
	//	h->SetBinError( h->FindBin( i, j ), 1e-10 );
	h->SetBinContent( i, j, 1e-10 );
	h->SetBinError( i, j, 1e-10 );
      }
    }
  }
}

double chi2( TH1* h0, TH1* h1 ) { 
  double c2 = 0;
  for ( int i=0 ; i<h0->GetNbinsX() ; i++ ) {
 
    double d0 = h0->GetBinContent(i+1);
    double d1 = h1->GetBinContent(i+1);
 
    double e0 = h0->GetBinError(i+1);
    double e1 = h1->GetBinError(i+1);

    double  e2 = e0*e0+e1*e1;
    
    if ( e2>0 ) c2 += (d0-d1)*(d0-d1)/e2; 

  }
  
  return c2;
}


int main(int argc, char** argv) { 

  if ( argc<4 ) return usage(argv[0], -1);

  /// control stuff
  
  std::string tag = "";
  std::string key = "";

  std::string dir = "";

  std::string ftestname = "";
  std::string frefname  = "";

  TFile* _ftest = 0;
  TFile* _fref  = 0;

  bool   effset = false;
  double effmin =  90;
  double effmax = 102;


  /// control flags and labels

  std::vector<std::string> usrlabels;
  bool uselabels    = false;
  bool addinglabels = false;

  std::vector<std::string> taglabels;
  bool addingtags  = false;

  bool make_ref_efficiencies = false;
  bool refit_resplots        = false;

  bool _bayes      = true;
  bool nopng       = false;
  bool nostats     = false;
  bool nomeans     = false;
  bool noref       = false;
  bool atlasstyle  = false;
  bool deleteref   = false;
  bool nowatermark = false;
  bool noplots     = false;
  bool Cfile       = false;
  bool notitle     = false;
  bool dochi2      = false;
  bool normref     = false;
  bool scalepix    = true;

  std::string atlaslabel = "for approval";

  double scale_eff     = -1;
  double scale_eff_ref = -1;
  
  std::string configfile = "";

  double xoffset = 0;

  
  double _ypos = 0;

  std::string pattern = "";
  std::string regex   = "";

  std::vector<std::string> chains;
  std::vector<std::string> refchains;

  bool addingrefchains = false;

  for(int i=1; i<argc; i++){
    std::string arg  = argv[i];

    if ( arg.find("-")!=0 && addinglabels ) { 
      usrlabels.push_back( arg );
      continue;
    }
    else addinglabels = false;

    if ( arg.find("-")!=0 && addingrefchains ) { 
      refchains.push_back( arg );
      continue;
    }
    else addingrefchains = false;

    if ( arg.find("-")!=0 && addingtags ) { 
      taglabels.push_back( arg );
      continue;
    }
    else addingtags = false;

    if ( arg=="-h" || arg=="--help" ) { 
       return usage(argv[0], 0);
    }
    else if ( arg=="-c" || arg=="--config" ) { 
      if ( ++i<argc ) configfile=argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg=="-t" || arg=="--tag" ) { 
      if ( ++i<argc ) tag=std::string("-")+argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg=="-l" || arg=="--labels" ) { 
      addinglabels = true;
    }
    else if ( arg=="-k" || arg=="--key" ) { 
      if ( ++i<argc ) key=argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg=="-d" || arg=="--dir" ) { 
      if ( ++i<argc ) dir=argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg=="--taglabels" ) { 
      addingtags = true;
    }
    else if ( arg=="--unscalepix" ) { 
      scalepix = false;
    }
    else if ( arg=="-yrange" ) { 
      effset = true;
      if ( ++i<argc ) effmin=std::atof(argv[i]);
      else return usage(argv[0], -1);
      if ( ++i<argc ) effmax=std::atof(argv[i]);
      else return usage(argv[0], -1);
    }
    else if ( arg=="-e" || arg=="--efficiencies" ) { 
      make_ref_efficiencies = true;
    }
    else if ( arg=="-r" || arg=="--refit" ) { 
      refit_resplots = true;
    }
    else if ( arg=="-nw" || arg=="--nowatermark" ) {
      nowatermark = true;
      Plots::setwatermark(!nowatermark);
    }
    else if ( arg=="--chi2" ) { 
      dochi2 = true;
    }
    else if ( arg=="-ns" || arg=="--nostats" ) { 
      nostats = true;
    }
    else if ( arg=="-nm" || arg=="--nomeans" ) { 
      nomeans = true;
    }
    else if ( arg=="-nt" || arg=="--notitle" ) { 
      notitle = true;
    }
    else if ( arg=="-nr" || arg=="--noref" ) { 
      Plotter::setplotref(false);
      noref = true;
    }
    else if ( arg=="--normref" ) { 
      normref = true;
    }
    else if ( arg=="-rc" || arg=="--refchains" ) { 
      addingrefchains = true;
    }
    else if ( arg=="-nb" || arg=="--nobayes" ) { 
      _bayes = false;
    }
    else if ( arg=="-es" || arg=="--effscale" ) { 
      if ( ++i<argc ) scale_eff=std::atof(argv[i]);
      else return usage(argv[0], -1);
    } 
    else if ( arg=="-er" || arg=="--effscaleref" ) { 
      if ( ++i<argc ) scale_eff_ref=std::atof(argv[i]);
      else return usage(argv[0], -1);
    } 
    else if ( arg=="-np" || arg=="--noplots" ) { 
      noplots = true;
    }
    else if ( arg=="-C" || arg=="--Cfiles" ) { 
      Cfile = true;
    }
    else if ( arg=="--deleteref" ) { 
      deleteref = true;
    }
    else if (              arg=="--nopng" ) { 
      nopng = true;
    }
    else if ( arg=="-as" || arg=="--atlasstyle" ) { 
      atlasstyle = true;
    }
    else if ( arg=="-al" || arg=="--atlaslabel" ) { 
      if ( ++i<argc ) atlaslabel=argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg=="-xo" || arg=="--xoffset" ) { 
      if ( ++i<argc ) xoffset=std::atof(argv[i]);
      else return usage(argv[0], -1);
    }
    else if ( arg=="-yp" || arg=="--ypos" ) { 
      if ( ++i<argc ) _ypos=std::atof(argv[i]);
      else return usage(argv[0], -1);
    }
    else if ( arg=="-s" || arg=="--swap" ) { 
      if ( ++i<argc ) pattern=argv[i];
      else return usage(argv[0], -1);
      if ( ++i<argc ) regex=argv[i];
      else return usage(argv[0], -1);
    }
    else if ( arg.find("-")==0 ) {
      std::cerr << "unknown option: " << arg << "\n" << std::endl;
      return usage(argv[0], -4);
    }
    else { 
      if ( _ftest==0 ) { 
	if ( exists(arg) ) { 
	  ftestname = arg;
	  // _ftest = new TFile( ftestname.c_str() );
	  _ftest = TFile::Open( ftestname.c_str() );

	}
	else { 
	  std::cerr << "main(): test file " << arg << " does not exist" << std::endl;
	  return -1;
	}	    
      }
      else if ( _fref==0 ) { 
	if ( exists(arg) ) {
	  frefname = arg;
	  // _fref = new TFile( frefname.c_str() );
	  _fref = TFile::Open( frefname.c_str() );
	}
	else { 
	  if ( _ftest ) delete _ftest;
	  std::cerr << "main(): ref file  " << arg << " does not exist" << std::endl;
	  return -1;
	}	    
      }
      else { 
	std::string chain = arg;
	replace ( chain, ":", "_" );
	replace ( chain, ";", "_" );
	chains.push_back(chain);
      }
    }
  }

  if ( scale_eff     == -1 ) scale_eff     = 100;
  if ( scale_eff_ref == -1 ) scale_eff_ref = scale_eff;


  if ( refchains.size()>0 && refchains.size()!=chains.size() ) return usage(argv[0], -1);
  
  if ( refchains.size()==0 ) refchains = chains;
  

  std::cout << argv[0] << " options:" << std::endl;
  std::cout << "\tATLAS style:                 " << ( atlasstyle ? "true" : "false" ) << std::endl; 
  std::cout << "\tBayesian uncertainties:      " << ( _bayes ? "true" : "false" ) << std::endl; 
  std::cout << "\trefit resplot uncertainties: " << ( refit_resplots ? "true" : "false" ) << std::endl; 
  std::cout << "\tsuppress mean and rms stats: " << ( nostats ? "true" : "false" ) << std::endl;  
  if ( !nostats ) std::cout << "\tsuppress meanstats:          " << ( nomeans ? "true" : "false" ) << std::endl;  
  std::cout << "\tsuppress png output:         " << ( nopng ? "true" : "false" ) << std::endl;  
  std::cout << "\tsuppress reference output:   " << ( noref ? "true" : "false" ) << std::endl;  
  if ( usrlabels.size()>0 )    std::cout << "\tlabels:       " << usrlabels << std::endl;  
  if ( taglabels.size()>0 )    std::cout << "\textra text:   " << taglabels << std::endl;  


  if ( atlasstyle ) { 
    SetAtlasStyle();
    
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.15);
    
  }
  else {
    gStyle->SetPadLeftMargin(0.105);
    gStyle->SetPadBottomMargin(0.105);
  }

  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.05);

  if ( chains.size()<4 ) { 
    std::cout << "Chains: " << chains << std::endl;
  }
  else {  
    std::cout << "Chains: " << std::endl;
    for ( unsigned ic=0 ; ic<chains.size() ; ic++ ) std::cout << "\t" << chains[ic] << std::endl;
  }    

  if ( usrlabels.size()>0 ) std::cout << "labels " << usrlabels << std::endl;

  if ( usrlabels.size()>0 && usrlabels.size()==chains.size() ) uselabels = true;


  /// get release data

  
  TTree*   dataTree = 0;
  TString* releaseData = new TString("");
  std::vector<std::string> release_data;
  
  if ( !nowatermark ) { 
    dataTree = (TTree*)_ftest->Get("dataTree");
  
    if ( dataTree ) { 
      dataTree->SetBranchAddress( "ReleaseMetaData", &releaseData);
      
      for (unsigned int i=0; i<dataTree->GetEntries() ; i++ ) {
	dataTree->GetEntry(i);      
	release_data.push_back( releaseData->Data() );
	std::cout << "main() release data: " << release_data.back() << " " << *releaseData << std::endl;
      }
    }
  
  
    if ( release_data.size()>0 ) { 
      if ( release_data.size()>1 ) std::cerr << "main() more than one release - using only the first" << std::endl;  
      
      //    std::cout << "release: " << chop(release_data[0], " " ) << std::endl;
      //    std::cout << "release: " << chop(release_data[0], " " ) << std::endl;
      
      release += "  (" + chop(release_data[0], " " );
      release += " " + chop(release_data[0], " " ) + ")";
    }
  }
  
  
  // Make output directory                                                                                                                           
  if (dir != "") {
    dir += "/";
    if ( mkdir( dir.c_str(), 0777 ) ) std::cerr << "main() couldn't create directory " << dir << std::endl;
    else                              std::cout << "main() output will be sent to directory " << dir << std::endl; 
  }

  TFile& ftest = *_ftest;
  TFile& fref  = *_fref;

  std::vector<std::string> chainnames;


  chainnames.resize(chains.size());
  chainnames.clear();

  int NeventTest=0;
  int NeventRef=0;


  std::vector<std::string> savedhistos;

  
  if ( !nowatermark ) { 
    TH1D* htestev = (TH1D*)ftest.Get("event") ;
    TH1D* hrefev  = (TH1D*)fref.Get("event") ;

    std::cout << "htestev " << htestev << " " << hrefev  << std::endl;
    
    if ( htestev ) NeventTest = htestev->GetEntries();
    if ( hrefev )  NeventRef  = hrefev->GetEntries();
  
    savedhistos.push_back("event");
  }
  else { 
    NeventTest = 1; /// get from Chains? 
    NeventRef  = 1; 
  }



  int    Nhistos = 48;
  const  int __Nhistos = 48;
  std::string __histos[__Nhistos][6] = { 

    /// distributions - 4
    //  { "pT",  "p_{T}",     "xaxis:lin:0.7:100",  "Offline p_{T} [GeV]",   "yaxis:log:auto",  ""  },
    { "pT",      "p_{T}",     "xaxis:lin:auto:1:100",     "Offline p_{T} [GeV]",   "yaxis:log:auto",  ""  },
    { "pT_rec",  "p_{T} rec", "xaxis:lin:auto:1:100",   "Trigger p_{T} [GeV]",   "yaxis:log:auto",  ""  },
    { "a0",      "a0",        "xaxis:lin:-2:2",     "Offline a_{0} [mm]",    "yaxis:log:auto",  ""  },
    { "a0_rec",  "a0 rec",    "xaxis:lin:-2:2",     "Trigger a_{0} [mm]",    "yaxis:log:auto",  ""  },
    { "z0",      "z0",        "xaxis:lin:-250:250", "z_{0} [mm]",            "yaxis:log:auto",  ""  },

    /// efficiencies - 10 
    //    { "pT_eff", "Efficiency p_{T}", "xaxis:log:0.7:100",     "Offline p_{T} [GeV]",          "yaxis:lin:90:102",       "Efficiency [%]" },       
    { "pT_eff",       "Efficiency p_{T}", "xaxis:log:auto:1:100",        "Offline track p_{T} [GeV]",    "yaxis:lin:auto:90:102",  "Efficiency [%]" },       
    { "eta_eff",      "Efficiency #eta",  "xaxis:lin",             "Offline track #eta",           "yaxis:lin:auto:90:102",  "Efficiency [%]" },       
    { "phi_eff",      "Efficiency #phi",  "xaxis:lin",             "Offline track #phi",           "yaxis:lin:auto:90:102",  "Efficiency [%]" },       
    { "d0_eff",       "Efficiency d0",    "xaxis:lin:autosym",     "Offline track d_{0} [mm]",     "yaxis:lin:auto:90:102",  "Efficiency [%]" },       
    //    { "a0_eff",  "Efficiency a0",   "xaxis:lin:-2:2",        "Offline track d_{0} [mm]",     "yaxis:lin:90:102",       "Efficiency [%]" },        
    { "a0_eff",       "Efficiency a0",    "xaxis:lin:autosym",     "Offline track d_{0} [mm]",     "yaxis:lin:auto:90:102",  "Efficiency [%]" },      
    { "z0_eff",       "Efficiency z0",    "xaxis:lin:-250:250",    "Offline track z_{0} [mm]",     "yaxis:lin:auto:90:102",  "Efficiency [%]" },        
 
    { "eff_vs_mu",    "Efficiency <#mu>",            "xaxis:lin:auto",       "<#mu>",              "yaxis:lin:90:102",       "Efficiency [%]" },       
    { "roi_dphi_eff", "Efficiency #Delta#phi(RoI)",  "xaxis:lin:-0.6:0.6",   "#Delta#phi (RoI)",   "yaxis:lin:90:102",       "Efficiency [%]" },
    { "roi_deta_eff", "Efficiency #Delta#eta(RoI)",  "xaxis:lin:-0.6:0.6",   "#Delta#eta (RoI)",   "yaxis:lin:90:102",       "Efficiency [%]" },       
    { "roi_dR_eff",   "Efficiency #DeltaR(RoI)",     "xaxis:lin:0:0.6",      "#Delta R (RoI)",     "yaxis:lin:90:102",       "Efficiency [%]" },       

    /// standard residuals - 5 
    { "ipT_res",    "Residual 1/p_{T}",  "xaxis:lin:-0.15:0.2",     "#Delta 1/p_{T} [GeV^{-1}]",    "yaxis:log:auto",    "Normalised entries" },       
    { "eta_res",    "Residual #eta",     "xaxis:lin:-0.05:0.05",    "#Delta#eta",                   "yaxis:log:auto",    "Normalised entries" },       
    { "phi_res",    "Residual #phi",     "xaxis:lin:-0.05:0.05",    "#Delta#phi",                   "yaxis:log:auto",    "Normalised entries" },
    //  { "z0_res", "Residual z0",       "xaxis:lin:-7:10",         "#Delta z_{0} [mm]",            "yaxis:lin:0:0.035", "Normalised entries" },
    { "z0_res",     "Residual z0",       "xaxis:lin:-10:10",        "#Delta z_{0} [mm]",            "yaxis:log:auto",    "Normalised entries" },
    { "a0_res",     "Residual a0",       "xaxis:lin:-1:1",          "#Delta d_{0} [mm]",            "yaxis:log:auto",    "Normalised entries" },       
 
    /// residuals vs track parameters - 17
    //    { "rd0_vs_pt/sigma",    "Residual d vs p_{T}",          "xaxis:lin:0:100",     "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "d_{0} resolution [mm]" },
    { "rd0_vs_pt/sigma",          "Residual d vs p_{T}",          "xaxis:lin:auto",      "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "d_{0} resolution [mm]" },
    { "rd0_vs_signed_pt/sigma",   "Residual d vs signed p_{T}",   "xaxis:lin:-100:100",  "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "d_{0} resolution [mm]" },    
    // { "rd0_vs_ABS_pt/sigma",   "Residual d vs absolute p_{T}", "xaxis:log:1:100",     "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "d_{0} resolution [mm]" },       
    { "rd0_vs_eta/sigma",         "Residual d vs #eta",           "xaxis:lin",           "Offline #eta",          "yaxis:lin:auto",  "d_{0} resolution [mm]" },                    
    // { "rd0_vs_eta/sigma",      "Residual d vs #eta",           "xaxis:lin",           "Offline #eta",          "yaxis:lin:0.015:0.05",  "d_{0} resolution [mm]" },                    
    { "rd0_vs_ipt/sigma",         "Residual d vs 1/p_{T}",        "xaxis:lin",           "1/p_{T} [GeV^{-1}]",    "yaxis:lin:auto",  "d_{0} resolution [mm]" },                       

    //    { "ript_vs_pt/sigma",   "Residual 1/p_{T} vs p_{T}",    "xaxis:lin:1:100",    "Offline p_{T} [GeV]",   "yaxis:log:auto",  "1/p_{T} resolution [GeV^{-1}]" },
    { "ript_vs_pt/sigma",         "Residual 1/p_{T} vs p_{T}",    "xaxis:lin:auto",     "Offline p_{T} [GeV]",   "yaxis:log:auto",  "1/p_{T} resolution [GeV^{-1}]" },
    { "ript_vs_eta/sigma",        "Residual 1/p_{T} vs #eta",     "xaxis:lin",          "Offline #eta",          "yaxis:lin:auto",  "1/p_{T} resolution [GeV^{-1}]" },
    { "ript_vs_ipt/sigma",        "Residual 1/p_{T} vs 1/p_{T}",  "xaxis:lin",          "1/p_{T} [GeV^{-1}]",    "yaxis:lin:auto",  "1/p_{T} resolution [GeV^{-1}]" },
    //   { "rpt_vs_ipt/sigma",    "Residual p_{T} vs 1/p_{T}",    "xaxis:lin",          "1/p_{T} [GeV^{-1}]",    "",                "p_{T} resolution [GeV]" },      

    //    { "reta_vs_pt/sigma",   "Residual #eta p_{T}",          "xaxis:log:1:100",    "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "#eta resolution" },            
    { "reta_vs_pt/sigma",         "Residual #eta p_{T}",          "xaxis:log:auto",     "Offline p_{T} [GeV]",   "yaxis:lin:auto",  "#eta resolution" },            
    { "reta_vs_eta/sigma",        "Residual #eta vs #eta",        "xaxis:lin",          "Offline #eta",          "yaxis:lin:auto",  "#eta resolution" },            
    //   { "reta_vs_eta/sigma",   "Residual #eta vs #eta",        "xaxis:lin",          "Offline #eta",          "yaxis:lin:1e-4:1.4e-3",  "#eta resolution" },            
    { "reta_vs_ipt/sigma",        "Residual #eta vs 1/p_{T}",     "xaxis:lin",          "1/p_{T} [GeV^{-1}]",    "yaxis:lin:auto",  "#eta resolution" },            

    //    { "rphi_vs_pt/sigma",   "Residual #phi vs p_{T}",       "xaxis:lin:auto",     "p_{T} [GeV]",           "yaxis:lin:auto", "#phi resolution" },
    { "rphi_vs_pt/sigma",         "Residual #phi vs p_{T}",       "xaxis:lin:1:100",    "p_{T} [GeV]",           "yaxis:lin:auto", "#phi resolution" },
    { "rphi_vs_ipt/sigma",        "Residual #phi vs 1/p_{T}",     "xaxis:lin",          "1/p_{T} [GeV^{-1}]",    "yaxis:lin:auto", "#phi resolution" },

    { "rzed_vs_eta/sigma",        "Residual z vs #eta",           "xaxis:lin",          "Offline #eta",          "yaxis:lin:auto", "z_{0} resolution [mm]" },
    //    { "rzed_vs_pt/sigma",   "Residual z vs p_{T}",          "xaxis:log:1:100",    "Offline p_{T} [GeV]",   "yaxis:lin:auto", "z_{0} resolution [mm]" },
    { "rzed_vs_pt/sigma",         "Residual z vs p_{T}",          "xaxis:log:auto",     "Offline p_{T} [GeV]",   "yaxis:lin:auto", "z_{0} resolution [mm]" },
    { "rzed_vs_signed_pt/sigma",  "Residual z vs signed p_{T}",   "xaxis:lin:-100:100", "Offline p_{T} [GeV]",   "yaxis:lin:auto", "z_{0} resolution [mm]" },
    //  { "rzed_vs_ABS_pt/sigma", "Residual z vs absolute p_{T}", "xaxis:lin:1:100",    "Offline p_{T} [GeV]",   "yaxis:lin:auto", "z_{0} resolution [mm]" },
    { "rzed_vs_zed/sigma",        "Residual z vs z",              "xaxis:lin:-250:250", "Offline z [mm]",        "yaxis:lin:auto", "z_{0} resolution [mm]" },
    { "rzed_vs_ipt/sigma",        "Residual z vs 1/p_{T}",        "xaxis:lin",          "1/p_{T} [GeV^{-1}]",    "yaxis:lin:auto", "z_{0} resolution [mm]" },

    /// track multiplicity - 1
    { "ntracks_rec",             "number of reconstructed tracks", "xaxis:lin:auto",   "N trigger tracks",     "yaxis:log:auto", "Entries"  },

    /// hit multiplicity - 6
    { "npix_eta/mean",           "mean number of pixel hits",  "xaxis:lin",   "Offline #eta",   "yaxis:lin:3:6",  "Pixel hits"    },
    { "nsct_eta/mean",           "mean number of SCT hits",    "xaxis:lin",   "Offline #eta",   "yaxis:lin:7:10", "SCT clusters"  },

    //    { "npix_pt/mean",           "mean number of pixel hits",  "xaxis:lin:0.7:250",   "Offline p_{T} [GeV]",   "yaxis:lin:3:6",  "Pixel hits"    },
    //    { "nsct_pt/mean",           "mean number of SCT hits",    "xaxis:lin:0.7:250",   "Offline p_{T} [GeV]",   "yaxis:lin:7:10", "SCT clusters"  },

    //    { "npixh_pt/mean",           "mean number of pixel holes",  "xaxis:lin:0.7:250",   "Offline p_{T} [GeV]",   "yaxis:lin:-1:6",  "Pixel holes"    },
    //    { "nscth_pt/mean",           "mean number of SCT holes",    "xaxis:lin:0.7:250",   "Offline p_{T} [GeV]",   "yaxis:lin:-1:10", "SCT holes"  },

    { "npix_pt/mean",           "mean number of pixel hits",  "xaxis:lin:auto",   "Offline p_{T} [GeV]",   "yaxis:lin:3:6",  "Pixel hits"    },
    { "nsct_pt/mean",           "mean number of SCT hits",    "xaxis:lin:auto",   "Offline p_{T} [GeV]",   "yaxis:lin:7:10", "SCT clusters"  },

    { "npixh_pt/mean",           "mean number of pixel holes",  "xaxis:lin:auto",   "Offline p_{T} [GeV]",   "yaxis:lin:-1:6",  "Pixel holes"    },
    { "nscth_pt/mean",           "mean number of SCT holes",    "xaxis:lin:auto",   "Offline p_{T} [GeV]",   "yaxis:lin:-1:10", "SCT holes"  },


    /// chi2 and chi2 probability - 3
    { "Chi2prob/1d",           "Chi2 probability",       "xaxis:lin",         "track #chi^{2} Probability",   "yaxis:lin:auto",  "Entries"  },
    { "Chi2dof/1d",            "Chi2 per dof",           "xaxis:lin",         "track #chi^{2} per dof",       "yaxis:log:auto",  "Entries"  },
    { "Chi2/1d",               "Chi2",                   "xaxis:lin",         "Offline track #chi^{2}",       "yaxis:lin",       "Entries"  },

    { "Chi2prob/mean",         "Chi2 probability vs pT", "xaxis:log:auto",    "Offline p_{T} [GeV]",          "yaxis:lin",  "mean track #Chi^{2} Probability"  }

  };
  

  std::vector<std::vector<std::string> > _histos;

  /// read config in from a file if requested ...

  if ( configfile!="" ) { 
    if ( exists(configfile) ) { 

      std::cout << argv[0] << ":\treading histogram configuration from file " << configfile << std::endl; 
    
      ReadCards rc(configfile);
      
      std::vector<std::string> histos = rc.GetStringVector("histos");
      
      //  std::string histos[Nhistos];
      
      for ( unsigned i=0 ; i<histos.size() ; ) { 
	std::vector<std::string> duff;
	/// 6 entries are needed for every histogram ...
	for ( int j=0 ; j<6 && i<histos.size() ; j++, i++  ) duff.push_back( histos[i] );
	_histos.push_back( duff );
      }
      
      for ( unsigned i=0 ; i<_histos.size() ; i++ ) std::cout << "histos: " << i << "\t" << _histos[i] << std::endl;
      
    }
    else { 
      std::cerr << argv[0] << ":\t config file not found: " << configfile << std::endl;
      return -1;
    }
  }
  else { 
    for ( int i=0 ; i<__Nhistos ; i++ ) { 
      std::vector<std::string> duff;
      for ( int j=0 ; j<6 ; j++ ) duff.push_back( __histos[i][j] );
      _histos.push_back( duff );
    }
  }

  Nhistos = _histos.size();

  if ( _histos.size()==0 ) return usage(argv[0], -1);

  //  const int Nhistos2D = 2;
  //  int Nhistos2D = 0;

  std::vector<std::string> histos2D; 
  std::vector<std::string> histonames2D; 

  //  std::string> histos2D[2] = {
  //    "eta_phi_rec",
  //    "phi_d0_rec",
  //  };

  //  std::string histonames2D[Nhistos2D] = {
  //    "#phi vs #eta",
  //    "d_{0} vs #phi",
  //  };


  std::set<std::string> h2dset;

  for ( int i=0 ; i<Nhistos ; i++ ) {
    std::string hist = _histos[i][0];
    if ( contains( hist, "2d" ) && !contains( hist, "2dof" ) ) { 
      histos2D.push_back( _histos[i][0] ); 
      histonames2D.push_back( _histos[i][1] ); 
      h2dset.insert( hist );
    }
  }

  if ( fulldbg ) std::cout << __LINE__ << std::endl;

  //  Nhistos2D = histos2D.size();

  std::cout << "size of chains " << chains.size() << std::endl;

  if ( chains.size()<4 ) { 
    std::cout << "main() processing chains : " << chains << std::endl; 
  }
  else { 
    std::cout << "main() processing chains : " << std::endl;
    for ( unsigned ic=0 ; ic<chains.size() ; ic++ ) std::cout << "\t" << chains[ic] << std::endl;
  }


  if ( chains.empty() ) return 0;


  std::cout << "N chains " << chains.size() << "\t2d size " << histos2D.size() << std::endl; 


  /// create a better 2d colour pallette
  
  if ( true ) { 
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);


    const Int_t Number = 3;
    Double_t Red[Number] = { 0.00, 0.00, 1.00};
    Double_t Green[Number] = { 0.00, 5.00, 1.00};
    Double_t Blue[Number] = { 0.00, 0.50, 0.00};
    Double_t Length[Number] = { 0.00, 0.50, 1.00 };
    Int_t nb=50;
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  }
  else gStyle->SetPalette(1);

  if ( fulldbg ) std::cout << __LINE__ << std::endl;

  /// so we can easily limit the number of histograms drawn
  int _Nhistos = Nhistos;

  std::vector<std::string> histos(_Nhistos);


  if ( fulldbg ) std::cout << __LINE__ << std::endl;

  for ( int i=0 ; i<_Nhistos ; i++ ) histos[i] = _histos[i][0];

  if ( fulldbg )  std::cout << __LINE__ << std::endl; 



  double rightmargin = gStyle->GetPadRightMargin(); 
  gStyle->SetPadRightMargin(0.1); 


  for ( unsigned int j=0; j<chains.size(); j++)  {

      if ( fulldbg )  std::cout << __LINE__ << std::endl; 

      std::string refchain = fullreplace( refchains[j], pattern, regex );

      if ( fulldbg )  std::cout << __LINE__ << std::endl; 



      //      if ( chains[j]=="FTK_TrackParticle" ){
      { 
	
	for ( unsigned i=0 ; i<histos.size() ; i++ ) {

	  
	  if ( fulldbg )  std::cout << __LINE__ << std::endl;
	  
	  std::string hist = _histos[i][0];
	  
	  if ( !contains( hist, "2d" ) || contains( hist, "2dof" ) ) continue; 
	  
	  if ( fulldbg ) std::cout << __LINE__ << std::endl;
	  
	  TH2F* htest2D = (TH2F*)ftest.Get((chains[j]+"/"+hist ).c_str()) ;
	  TH2F* href2D  = (TH2F*)fref.Get((refchain+"/"+hist ).c_str()) ;
	  
	  
	  if ( htest2D==0 || href2D==0 ) continue;

	  zero( htest2D );
	  zero( href2D );
	  
	  AxisInfo xinfo( _histos[i][2] ); 
	  AxisInfo yinfo( _histos[i][4] ); 

	  //	      std::cout << "yinfo " << _histos[i][0] << ": " << _histos[i][4] << " :: " << yinfo << " " << yinfo.rangeset() << std::endl;

	  htest2D->GetYaxis()->SetTitle( _histos[i][5].c_str() );
	  href2D->GetYaxis()->SetTitle(  _histos[i][5].c_str() );
	  
	  htest2D->GetXaxis()->SetTitle( _histos[i][3].c_str() );
	  href2D->GetXaxis()->SetTitle(  _histos[i][3].c_str() );
	
	  if ( yinfo.rangeset() ) { 		
	    htest2D->GetYaxis()->SetRangeUser( yinfo.lo(), yinfo.hi() );
	    href2D->GetYaxis()->SetRangeUser( yinfo.lo(), yinfo.hi() );
	  }

	  if ( xinfo.rangeset() ) { 		
	    htest2D->GetXaxis()->SetRangeUser( xinfo.lo(), xinfo.hi() );
	    href2D->GetXaxis()->SetRangeUser( xinfo.lo(), xinfo.hi() );
	  }
	   

	  savedhistos.push_back( chains[j]+"/"+hist );
	  
	  TCanvas* c1 = new TCanvas(label( "canvas2D-%d", i ).c_str(),"histogram",1200,500);
	  c1->Divide(2,1);
	  c1->cd(1);
	  
	  std::string title = hist;
	  title += " test";
	  htest2D->SetTitle(title.c_str());
	  //              htest2D->Scale(1./NeventTest);
	  htest2D->Draw("COLZ");

	  c1->cd(2);

	  std::string title2 = hist;
	  title2 += "reference";
	  href2D->SetTitle(title2.c_str());
	  // gPad->SetLogz();                                                                                         
	  href2D->Draw("COLZ");

	  //	  gPad->SetLogz(true);

	  std::string plotname = chains[j];
	  plotname += "_";
	  //	  plotname += "HLT_";
	  plotname += hist.c_str();
	  plotname += ".png";
	  replace( plotname, "/", "_");

	  plotname = fullreplace( plotname, "InDetTrigTrackingxAODCnv", "" );
                               
	  std::cout << "plot 2D " << dir+plotname << std::endl;
	  
       	  c1->Print((dir+plotname).c_str());
	  
	  delete c1;
	}
      }
  }

  gStyle->SetPadRightMargin(rightmargin); 

  
  //  std::exit(0);

  /// so we can easily limit the number of histograms drawn
  ///  int _Nhistos = Nhistos;

  ///  std::vector<std::string> histos(_Nhistos);

  for ( int i=0 ; i<Nhistos ; i++ ) {

    if ( fulldbg )  std::cout << __LINE__ << std::endl;

    //    if ( i!=17 ) continue;

    //    histos[i] = _histos[i][0];

    /// don't process 2d histograms juet yet
    if ( h2dset.find( histos[i] )!=h2dset.end() ) continue; 

    //    std::cout << "main() processing histos[" << i <<  "] " << (i<10 ? " " : "" ) << (chains[0]+"/"+histos[i]) << std::endl;
    std::cout << "main() processing histos[" << i <<  "] " << (i<10 ? " " : "" ) << histos[i] << "\t" << AxisInfo( _histos[i][2] ) << std::endl;

    //    std::vector<Plotter> plots;

    Plots plots_eff;
    plots_eff.clear();

    Plots plots;
    plots.clear();
  
    TCanvas* c1 = new TCanvas(label("canvas-%d",i).c_str(),"histogram",800,600);
    c1->cd();

    


    /// legends ....
    Legend legend;
    Legend legend_eff;


    double xpos  = 0.15;
    double ypos  = 0.93;

    if ( contains(histos[i],"eff") || contains(histos[i],"Eff_") ) ypos = 0.15;

    if ( atlasstyle ) { 
      xpos  = 0.18;
      if ( ypos>0.5 ) ypos = 0.85; 
      else            ypos = 0.18;
    }

    if ( _ypos!=0 ) ypos = _ypos;

    double xpos_original = xpos;

    xpos += xoffset;

    /// calculate all the postions for the items in the legend

    int Nlines = chains.size() + taglabels.size();

    std::vector<double> ypositions;

    double deltay = (chains.size()*0.06-0.005)/chains.size();

    double ylo = ypos;
    double yhi = ypos;
    
    if ( ypos>0.5 ) ylo -= Nlines*deltay;
    else            yhi += Nlines*deltay;
    
    for ( int ilines=0 ; ilines<Nlines ; ilines++ ) { 
      ypositions.push_back( yhi - deltay*(ilines+0.5) );
    }


    // efficiencies or residuals?
    //    if ( contains(histos[i],"eff") ) {
    legend     = Legend( xpos, xpos+0.1, ylo, ylo+chains.size()*0.06-0.005 );
    legend_eff = Legend( xpos, xpos+0.1, ylo, ylo+chains.size()*0.06-0.005 );
    //    }
    //    else {
    //      legend     = Legend( xpos, xpos+0.1, ypos-chains.size()*0.06-0.01, ypos );
    //      legend_eff = Legend( xpos, xpos+0.1, ypos-chains.size()*0.06-0.01, ypos  );
    //    }


    std::vector<std::string> Mean;
    std::vector<std::string> RMS;
  
    Mean.clear();
    RMS.clear();

    std::vector<std::string> Chi2;
    std::vector<std::string> MeanRef;
    std::vector<std::string> RMSRef;
  
    //    int colours[6] = { 1, 2, 4, 6, 7, 8 };

    Chi2.clear();
    MeanRef.clear();
    RMSRef.clear();

    std::string plotname = ""; 


    int mean_power = 0;
    int rms_power  = 0;
    bool power_set = false;;
    

    std::string xaxis = _histos[i][3];
    std::string yaxis = _histos[i][5];

    AxisInfo xinfo( _histos[i][2] ); 
    AxisInfo yinfo( _histos[i][4] ); 

    //    bool uselogx = xinfo.log();
    //    bool uselogy = yinfo.log();

    for ( unsigned int j=0; j<chains.size(); j++)  {

      std::cout << "       processing chain[" << j << "]   " << chains[j] << std::endl;
   
      TH1F* htest = 0;
      TH1F* href  = 0;

      TH1F* htestnum = 0;
      TH1F* htestden = 0;
      TH1F* hrefnum  = 0;

      TGraphAsymmErrors* tgtest = 0;


      std::string refchain = fullreplace( refchains[j], pattern, regex );

  
      /// refit the resplots - get the 2d histogram and refit
     
      if ( refit_resplots && contains(histos[i],"/sigma") ) { 

    	    std::cout << "       refitting:  " << histos[i] << std::endl;

	    Resplot::setoldrms95(false);
	    Resplot::setscalerms95(true);

	    std::string _tmp = histos[i];
	    std::string base = chop( _tmp, "/sigma" );
  
	    TH2D* _htest2d = (TH2D*)ftest.Get((chains[j]+"/"+base+"/2d").c_str()) ;
	    TH2D* _href2d = (TH2D*)ftest.Get((refchain+"/"+base+"/2d").c_str()) ;

	    savedhistos.push_back( chains[j]+"/"+base+"/2d" );


	    if ( _htest2d==0 || ( !noref && _href2d==0 ) ) continue;

	    /// get the test histogram

	    //	    std::cout << "test " << _htest2d << std::endl;
	    Resplot rtest("tmp", _htest2d);
	    rtest.Finalise(Resplot::FitNull95);
	    htest = (TH1F*)rtest.Sigma()->Clone("rtest_sigma"); htest->SetDirectory(0);

	    /// NB: DON'T Refit the reference, since only the central values
	    ///     are plotted
	    //	 std::cout << "ref " << _href2d << std::endl;
	    //   Resplot rref("tmp", _href2d);
	    //   rref.Finalise(Resplot::FitNull95);
	    //   href = (TH1F*)rref.Sigma()->Clone("rref_sigma"); href->SetDirectory(0);

	    /// still get the reference histo

	    //	    href  = (TH1F*)fref.Get((chains[j]+"/"+histos[i]).c_str()) ;
	    TH1F* hreft  = (TH1F*)fref.Get( (refchain+"/"+histos[i]).c_str() );

	    if ( htest==0 || hreft==0 ) { 
	      std::cerr << "missing histogram: " << (refchain+"/"+histos[i]) << " " << hreft << std::endl; 
	      continue;
	    }

	    href = (TH1F*)hreft->Clone();
	    href->SetDirectory(0);

	    std::cout << "\tget " << (refchain+"/"+histos[i]) << "\t" << href << std::endl;

	    savedhistos.push_back( refchain+"/"+histos[i] );

      }
      else { 

	if ( fulldbg ) std::cout << __LINE__ << std::endl;

	/// everything else 
       	
	std::string reghist = histos[i];

	std::cout << "chains: " << (chains[j]+"/"+reghist) << std::endl;

        htest = (TH1F*)ftest.Get((chains[j]+"/"+reghist).c_str()) ;

	savedhistos.push_back( chains[j]+"/"+reghist );

	TH1F* hreft  = (TH1F*)fref.Get((refchain+"/"+reghist).c_str()) ;


	if ( htest==0 || hreft==0 ) { 
	  if ( htest==0 ) std::cerr << "missing histogram: " << (chains[j]+"/"+reghist) << " " << htest<< std::endl; 
	  if ( hreft==0 ) std::cerr << "missing histogram: " << (refchain+"/"+reghist)  << " " << hreft << std::endl; 
	  continue;
	}
	
	if ( fulldbg ) std::cout << "htest: " << htest << std::endl; 
	if ( fulldbg ) std::cout << "hreft: " << hreft << std::endl; 

	href = (TH1F*)hreft->Clone();
	href->SetDirectory(0);


	std::cout << "\tget " << ( refchain+"/"+reghist) << "\thref  " << href << std::endl;
	std::cout << "\tget " << (chains[j]+"/"+reghist) << "\thtest " << htest << std::endl;

	if ( htest==0 || ( !noref && href==0 ) ) continue;


	if ( fulldbg ) std::cout << __LINE__ << std::endl;
	
	if ( histos[i].find("1d")!=std::string::npos ) { 
	  if (        htest->GetNbinsX()>500 ) htest->Rebin(10);
	  if ( href && href->GetNbinsX()>500 ) href->Rebin(10);
	}


	if ( histos[i].find("zed_eff")!=std::string::npos ) { 
	  if (        htest->GetNbinsX()>100 ) htest->Rebin(5);
	  if ( href && href->GetNbinsX()>100 ) href->Rebin(5);
	}


        if ( fulldbg ) std::cout << __LINE__ << std::endl;

	if ( scalepix && std::string(htest->GetName()).find("npix")!=std::string::npos ) htest->Scale(0.5);

	if ( fulldbg ) std::cout << __LINE__ << std::endl;

	if ( notitle ) { 
	  htest->SetTitle("");
	  if( href ) href->SetTitle("");
	} 

        if ( fulldbg ) std::cout << __LINE__ << std::endl;

      }

      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      if ( make_ref_efficiencies ) { 
	if ( htest && href ) { 

	  //	  std::cout << "contains _eff " << contains( std::string(htest->GetName()), "eff" ) << std::endl;

	  if ( contains( std::string(htest->GetName()), "_eff" ) ) {

	    htestnum = (TH1F*)ftest.Get((chains[j]+"/"+histos[i]+"_n").c_str()) ;

	    TH1F* hrefnumt  = (TH1F*)fref.Get((refchain+"/"+histos[i]+"_n").c_str()) ;

	    hrefnum = (TH1F*)hrefnumt->Clone();
	    hrefnum->SetDirectory(0);


	    //	    savedhistos.push_back( chains[j]+"/"+histos[i]+"_n" );
	    savedhistos.push_back( refchain+"/"+histos[i]+"_n" );
	    
	    //	    std::cout << "numerator histos " << htestnum << " " << hrefnum << std::endl;

	  }
	}
      }


      if ( _bayes ) { 

	if ( htest && contains( std::string(htest->GetName()), "eff" ) ) {

	  //	  delete htest;

	  htestnum = (TH1F*)ftest.Get((chains[j]+"/"+histos[i]+"_n").c_str()) ;
	  htestden = (TH1F*)ftest.Get((chains[j]+"/"+histos[i]+"_d").c_str()) ;

	  savedhistos.push_back( chains[j]+"/"+histos[i]+"_n" );
	  savedhistos.push_back( chains[j]+"/"+histos[i]+"_d" );

	  std::cout << "Bayesian error calculation " << htestnum << " " << htestden << "\tscale " << scale_eff << std::endl;

	  if ( htestnum && htestden ) { 

#if 0
	    if ( contains( htest->GetName(), "_vs_lb" )  ) { 
	      std::cout << "rebin " << histos[i] << std::endl;
	      htestnum->Rebin(3);
	      htestden->Rebin(3);
	    }
#endif	
    
	    Efficiency e( htestnum, htestden, "", scale_eff );
	    tgtest = e.Bayes(scale_eff);

	    htest = e.Hist();

	  }
	
	  /// now recalculate reference

	  std::cout << "recalculating reference efficiencies ..." << std::endl; 

	  if ( href ) { 

	    //	  delete htest;

	    std::cout << "doin ..." << std::endl; 

	    TH1F* hrefnum = (TH1F*)fref.Get((refchain+"/"+histos[i]+"_n").c_str()) ;
	    TH1F* hrefden = (TH1F*)fref.Get((refchain+"/"+histos[i]+"_d").c_str()) ;

	    std::cout << "Bayesian error calculation " << htestnum << " " << htestden << "\tscale " << scale_eff     << std::endl;
	    std::cout << "Bayesian error calculation " << hrefnum  << " " << hrefden  << "\tscale " << scale_eff_ref << std::endl;
	    
	    if ( hrefnum && hrefden ) { 
	      Efficiency e( hrefnum, hrefden, "", scale_eff_ref );
	      // tgref = e.Bayes(scale_eff);

	      href = e.Hist();
	      
	    }
	  }
	}
	
      }

      
      if ( htest==0 ) { 
	std::cout << "       no test histogram :   " << (chains[j]+"/"+histos[i]) << std::endl;
	continue;
      }

      if ( !noref && href==0 ) { 
	std::cout << "       no ref histogram :    " << (chains[j]+"/"+histos[i]) << std::endl;
	continue;
      }


      htest->GetYaxis()->SetTitleOffset(1.5); 
      htest->GetXaxis()->SetTitleOffset(1.5); 
      htest->GetXaxis()->SetTitle(xaxis.c_str());
      htest->GetYaxis()->SetTitle(yaxis.c_str());

      if ( !noref ) { 
	href->GetYaxis()->SetTitleOffset(1.5);
	href->GetXaxis()->SetTitleOffset(1.5);
	href->GetXaxis()->SetTitle(xaxis.c_str());
	if ( contains(yaxis,"Efficiency") && !contains(yaxis,"%") && scale_eff==100 ) href->GetYaxis()->SetTitle((yaxis+" [%]").c_str());
	else href->GetYaxis()->SetTitle(yaxis.c_str());
      }	

      if ( fulldbg ) std::cout << __LINE__ << std::endl;

#if 1
      if ( contains(histos[i],"ntracks") ) {
	htest->Rebin(2);
        htest->Sumw2();
        if ( !noref ) { 
	  href->Rebin(2);
	  href->Sumw2();
	}
	//    htest->Scale(1./NeventTest);
	//    href->Scale(1./NeventRef);
      }
#endif

      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      /// only set the plot name once )  
      if ( plotname == "" ) { 

	if ( key!="" ) { 
	  htest->SetTitle("");
	  if ( href ) href->SetTitle("");
	  plotname = key+"_";
	}
	else if(contains(chains[j],"FTK")){
	  htest->SetTitle(("FTK "+ _histos[i][1]).c_str());
	  if ( href ) href->SetTitle(("FTK "+ _histos[i][1]).c_str());
	  plotname = "FTK_";
	}
	else if(fcontains(chains[j],"HLT_")){
	  htest->SetTitle("");
	  if ( href ) href->SetTitle("");
	  plotname = "HLT_";
	}
	else if(fcontains(chains[j],"EF_")){
	  htest->SetTitle("");
	  if ( href ) href->SetTitle("");
	  plotname = "EF_";
	}
	else if(fcontains(chains[j],"L2_")){
	  htest->SetTitle("");
	  if ( href ) href->SetTitle("");
	  plotname = "L2_";
	}
	
	plotname += histos[i]; 

	/// replace the "/" in the filename so we don't try to 
	/// make plots in subdirectories by accident  
	replace(plotname, "/", "_"); 
	
      }


      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      bool residual = false;
      
      if ( contains(histos[i],"_res") ||  contains(histos[i],"residual_") || contains(histos[i],"1d") ) residual = true; 

      std::string c = chains[j];

      if ( c.find("_IDTrkNoCut")!=std::string::npos ) c.erase( c.find("_IDTrkNoCut"), 11 );
      if ( c.find("_idperf")!=std::string::npos )     c.erase( c.find("_idperf"), 7 );
      if ( c.find("_bperf")!=std::string::npos )      c.erase( c.find("_bperf"), 6 );
      if ( c.find("xAODCnv")!=std::string::npos )     c.erase( c.find("xAODCnv"), 7 );
      if ( c.find("Tracking")!=std::string::npos ) c.replace( c.find("Tracking"), 8, "Trk" );    

      //      if ( c.find("_EFID")!=std::string::npos )     c.erase( c.find("_EFID"), 5 );

      replace( c, "_Tr", " :  " );
      replace( c, "_In", " :  " );

      c = "  " + c;

      /// calculate and set axis limits

      //      std::cout << "adding plot " << histos[i] << " " << htest->GetName() << std::endl;

      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      if ( uselabels )  plots.push_back( Plotter( htest, href, usrlabels[j], tgtest ) );
      else              plots.push_back( Plotter( htest, href, c, tgtest ) );

      if ( fulldbg ) std::cout << __LINE__ << std::endl;


      if ( make_ref_efficiencies ) { 
	
	if ( htestnum && hrefnum ) { 
	  Efficiency e( htestnum, hrefnum, "", scale_eff );

	  TH1* h = e.Hist();

	  double range = h->GetMaximum()-h->GetMinimum();

	  if ( range<0.2*scale_eff ) {
 
	    double _max = int( (h->GetMaximum() + 20)*0.1 )*0.1*scale_eff;
	    double _min = int( (h->GetMinimum() - 10)*0.1 )*0.1*scale_eff;
	    
	    if ( _max>1*scale_eff ) _max = 1.02*scale_eff;
	    if ( _min<0 )           _min = 0;
	    
	    h->SetMinimum(_min);
	    h->SetMaximum(_max);
	   
	  }

	  plots_eff.push_back( Plotter( e.Hist(), 0, c ) );
	  
	}     
      }
      

      if ( href ) Chi2.push_back( label( "chi2 = %lf", chi2( htest, href ) ) );

      if( residual ) {
	
	/// resolutions 

	std::cout << "calculating resolutions : " << histos[i] << " " << htest->GetName() << std::endl;

	TF1* d95 = Resplot::FitNull95( (TH1D*)htest );
	
	double   mean_95 = d95->GetParameter(1);
	double  dmean_95 = d95->GetParError(1);
	double    rms_95 = d95->GetParameter(2);
	double   drms_95 = d95->GetParError(2);
	
	std::cout <<  "\t\t" << histos[i] 
		  << "\tmean:     " << mean_95 << " +- " << dmean_95 
		  << "\trms:     "  <<  rms_95 << " +- " << drms_95 << std::endl; 

	/// calculate power

	if ( !power_set ) { 
	  for ( int ip=-2 ; ip<9 ; ip++ ) { 
	    if ( std::fabs(mean_95) >= std::pow( 10., double(-ip) ) ) { 
	      mean_power = ip;
	      break;
	    }
	  }
	  
	  for ( int ip=-2 ; ip<9 ; ip++ ) { 
	    if ( std::fabs(rms_95)  >= std::pow( 10., double(-ip) ) ) { 
	      rms_power = ip;
	      break;
	    }
	  }
	}

	power_set = true;

	std::cout <<  "\t\t" << histos[i] 
		  << "\tmean:     " << mean_95 << " +- " << dmean_95 << " : pow " << mean_power
		  << "\trms:     "  <<  rms_95 << " +- " << drms_95  << " : pow " << rms_power << std::endl; 


	if ( mean_power == 0 ) { 
	  Mean.push_back(label("mean_{95}  = %4.2lf #pm %4.2lf", mean_95, dmean_95) );
	}
	else { 
	  Mean.push_back(label("mean_{95}  = ( %4.2lf #pm %4.2lf ) #times 10^{%d}", 
			       mean_95*std::pow(10.,double(mean_power)), dmean_95*std::pow(10,double(mean_power)), -mean_power ) );
	  
	}

	if ( rms_power == 0 ) { 
	  RMS.push_back(label( "rms_{95}   = %4.2lf #pm %4.2lf", rms_95,  drms_95 ) ); 
	}
	else { 
	  RMS.push_back(label( "rms_{95}   = ( %4.2lf #pm %4.2lf ) #times 10^{%d}", 
			       rms_95*std::pow(10.,double(rms_power)),  drms_95*std::pow(10,double(rms_power)), -rms_power ) );
	}

	if ( href ) { 
	  TF1* d95ref = Resplot::FitNull95( (TH1D*)href );
	  
	  double   mean_95ref = d95ref->GetParameter(1);
	  double  dmean_95ref = d95ref->GetParError(1);
	  double    rms_95ref = d95ref->GetParameter(2);
	  double   drms_95ref = d95ref->GetParError(2);
	  
	  std::cout <<  "\t\t" << histos[i]
		    << "\tmean ref: " << mean_95ref << " +- " << dmean_95ref << " : pow " << mean_power
		    << "\trms ref: "  <<  rms_95ref << " +- " << drms_95ref  << " : pow " << rms_power << std::endl;
	  
	  if ( mean_power == 0 ) { 
	    MeanRef.push_back(label("mean_{95} ref = %4.2lf #pm %4.2lf", mean_95ref, dmean_95ref) );
	  }
	  else { 
	    MeanRef.push_back(label("mean_{95} ref = ( %4.2lf #pm %4.2lf ) #times 10^{%d}", 
				    mean_95ref*std::pow(10,double(mean_power)), dmean_95ref*std::pow(10,double(mean_power)), -mean_power ) );
	    
	  }
	  
	  
	  if ( rms_power == 0 ) { 
	    RMSRef.push_back(label( "rms_{95}  ref = %4.2lf #pm %4.2lf", rms_95ref,  drms_95ref ) ); 
	  }
	  else { 
	    RMSRef.push_back(label( "rms_{95}  ref = ( %4.2lf #pm %4.2lf ) #times 10^{%d}", 
				    rms_95ref*std::pow(10,double(rms_power)),  drms_95ref*std::pow(10,double(rms_power)), -rms_power ) );
	  }
	}
	
        htest->Sumw2();
        if ( href ) href->Sumw2();
        htest->Scale(1./NeventTest);
        if ( href ) href->Scale(1./NeventRef);

      }
     
      if ( yinfo.normset() ) { 
	Norm( htest );
	if ( href ) Norm( href );
      }

      if ( yinfo.binwidth() ) { 
	binwidth( htest );
	if ( href ) binwidth( href );
      }
    


      if ( !noref && normref && !contains( histos[i], "mean") && !contains( histos[i], "sigma" ) && !contains( histos[i], "eff" ) ) { 
	double entries = Entries( htest );
	Norm( href, entries );
      }

    }

    if ( fulldbg ) std::cout << __LINE__ << std::endl;

    if ( !noplots ) { 

      /// try to localise all axis range setting, log, lin scales etc
      /// to this one place 

      if ( fulldbg ) std::cout << __LINE__ << std::endl;
    
      /// sort out the range settings for the xaxis ...
      plots.sortx( xinfo );

      if ( fulldbg ) std::cout << __LINE__ << std::endl;


      double  yminset = 0;
      double  ymaxset = 0;

      if ( yinfo.autoset() ) { 
	
	double rmin = 0;
	double rmax = 0;
	
	if ( xinfo.rangeset() ) { 
	  rmin = plots.realmin( plots.lo(), plots.hi() );
	  rmax = plots.realmax( plots.lo(), plots.hi() );
	}
	else {
	  rmin = plots.realmin();
	  rmax = plots.realmax();
	}
	
	if ( yinfo.log() && rmin!=0 && rmax!=0 ) { 

	  double delta = std::log10(rmax)-std::log10(rmin);

	  if ( atlasstyle ) ymaxset =  rmax*std::pow(10,delta*0.15*(chains.size()+taglabels.size()+1));
	  else              ymaxset =  rmax*std::pow(10,delta*0.15*(chains.size()+taglabels.size())); 

	  yminset =  rmin*std::pow(10,-delta*0.1);

	  if ( yminset!=yminset ) { 
	    std::cerr << " range error " << delta << " " << yminset << " " << ymaxset << "\t(" << rmin << " " << rmax << ")" << std::endl;
	    std::exit(-1);
	  }

	}
	else { 
	  double delta = rmax-rmin;
	  ymaxset = rmax+delta*0.1*chains.size();
	  yminset = rmin-delta*0.1;
	}
	
      }
      else {  
	if ( yinfo.rangeset() ) { 
	  yminset = yinfo.lo();
	  ymaxset = yinfo.hi();
	}
      }
      
      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      //      std::cout <<  "yauto: " << yinfo.autoset() << "\tyrange " << yinfo.rangeset() << std::endl;

      //      std::cout << "yminset " << yminset << "\tymaxset " << ymaxset << std::endl;  

      if ( yinfo.autoset() && yinfo.rangeset() ) {

	if ( fulldbg ) std::cout << __LINE__ << std::endl;

	if ( yminset>yinfo.lo() ) yminset = yinfo.lo();
	if ( ymaxset<yinfo.hi() ) ymaxset = yinfo.hi();
      }
      

      if ( contains(histos[i],"_eff") ) {

	if ( fulldbg ) std::cout << __LINE__ << std::endl;
      
	if ( effset ) { 
	  ymaxset = effmax;
	  yminset = effmin;
	}
      }
      
      if ( ymaxset!=0 || yminset!=0 ) {

	if ( fulldbg ) std::cout << __LINE__ << std::endl;

	plots.Max( ymaxset );
	plots.Min( yminset );
      }

      if ( fulldbg ) std::cout << __LINE__ << std::endl;

      if ( yminset>0 ) plots.SetLogy(yinfo.log());
      else             plots.SetLogy(false);
   
      //      plots.SetLogy(false);

      if ( fulldbg ) std::cout << __LINE__ << std::endl;
      
      ///    if ( contains(histos[i],"_res"))  plots.xrange(true);
      //      if ( contains(histos[i],"_res") ) plots.MaxScale( 100 ); 
      
      //    plots.limits();
      
      /// actually draw the plot here ...
      
      if ( fulldbg ) if ( fulldbg ) std::cout << __LINE__ << std::endl;

      plots.Draw( legend );

      if ( fulldbg ) if ( fulldbg ) std::cout << __LINE__ << std::endl;

      if ( atlasstyle ) ATLASLabel( xpos, ypositions[0]+deltay, atlaslabel.c_str() );

      if ( fulldbg ) if ( fulldbg ) std::cout << __LINE__ << std::endl;
          
      for ( unsigned it=0 ; it<taglabels.size() ; it++ ) { 
	//      std::cout << "\ttaglabel " << ypositions[it] << "\t(" << histos[i] << ")" << std::endl;
	DrawLabel( xpos, ypositions[it], taglabels[it], kBlack, 0.04 );
      }
    }
    
    if ( fulldbg ) if ( fulldbg ) std::cout << __LINE__ << std::endl;


    if ( ( !nostats || !nomeans ) && !noplots ) { 
      if ( dochi2 ) for ( unsigned  j=0 ; j<Chi2.size() ; j++ ) DrawLabel( 0.75, 0.85-j*0.035, Chi2[j], colours[j%6] );
      if ( (contains(histos[i],"_res") || 
	    contains(histos[i],"1d")   ||
	    histos[i]=="pT"            || 
	    contains(histos[i],"residual_") ||
	    contains(histos[i],"vs_pt") ) && !contains(histos[i],"sigma") ) { 

	if ( contains(histos[i],"_res") || contains(histos[i],"residual_") || contains(histos[i],"1d") ){
	  for ( unsigned j=0 ; j<chains.size() ; j++ ) { 
	    if ( !noref ) { 
	      if ( j<MeanRef.size() ) {
		if ( !nomeans ) DrawLabel( xpos_original-0.02, (0.57-j*0.035), MeanRef[j], colours[j%6] );
		DrawLabel( xpos_original-0.02, (0.57-0.035*chains.size()-j*0.035)-0.01, RMSRef[j],  colours[j%6] );
	      }
	    }
	    if ( j<Mean.size() ) {
	      if ( !nomeans ) DrawLabel( 0.62, (0.57-j*0.035), Mean[j],  colours[j%6] );
	      DrawLabel( 0.62, (0.57-0.035*chains.size()-j*0.035)-0.01, RMS[j],  colours[j%6] );
	      //	      std::cout << "\tdraw stats " << histos[i] << "\tMean " << Mean[j] << std::endl;
	    }
	  }
	}
      }
      
    }

#if 0
    if ( uselogx ) c1->SetLogx(true);
    else           c1->SetLogx(false);

    if ( uselogy ) c1->SetLogy(true);
    else           c1->SetLogy(false);
#endif

    if ( fulldbg ) std::cout << __LINE__ << std::endl;

    if ( !noplots ) { 

      if ( plotname!="" ) { 

	if ( Cfile ) plots.back().Print( dir+plotname+tag+".C" );

	//      plots.back().Print( dir+plotname+tag+".C" );
	plots.back().Print( dir+plotname+tag+".pdf" );

	if ( !nopng ) plots.back().Print( dir+plotname+tag+".png" );

      }    

      
      if ( make_ref_efficiencies ) { 
	
	plots_eff.SetXaxisTitle( plots.GetXaxisTitle() ); 
	plots_eff.SetYaxisTitle( plots.GetYaxisTitle() ); 
	
	plots_eff.Draw( legend_eff );
	
	if ( plotname!="" ) {

	  plots_eff.back().Print( dir+plotname+tag+"_refeff.pdf" );

	  if ( !nopng ) plots_eff.back().Print( dir+plotname+tag+"_refeff.png" );
	
	  gPad->SetLogx(true);

	}    
	
      }
      
      delete c1;
      
    }

    if ( fulldbg )  std::cout << __LINE__ << std::endl;
 
  }


  if ( fulldbg ) std::cout << __LINE__ << std::endl;
  
  /// if deleting all non-used reference histograms 

  if ( deleteref ) {
    
    if ( _fref ) { 
    
      /// clean up reference file
      
      std::cout << "main() cleaning up reference file" << std::endl;

      //      TFile* newout = new TFile(".newout.root","recreate"); 
      TFile* newout = new TFile(".newout.root","recreate"); 
      newout->cd();
      
      /// copy the release tree if need be 

      if ( dataTree ) dataTree->Write("dataTree");
      
      TDirectory* base = gDirectory;
      
      for ( unsigned i=0 ; i<savedhistos.size() ; i++ ) { 
	
	//	std::cout << i << " " << savedhistos[i] << std::endl;
	
	std::vector<std::string> dirs = AxisInfo::split( savedhistos[i], "/" );
	
	for ( unsigned j=0 ; j<dirs.size()-1 ; j++ ) { 
	  std::cout << "\t" << dirs[j] << std::endl;
	  TDirectory* renedir = gDirectory->GetDirectory( dirs[j].c_str() );
	  if ( renedir==0 ) gDirectory->mkdir( dirs[j].c_str() );
	  gDirectory->cd( dirs[j].c_str() );
	}
	
	TH1* href  = (TH1*)fref.Get( savedhistos[i].c_str() );
	if ( !noref && href ) {
	  std::cout << i << " " << savedhistos[i] << " 0x" << href << std::endl;
	  href->Write( dirs.back().c_str() );
	}
      
          
	base->cd();
      }
    
      newout->Close();
    
    }
  }


  if ( fulldbg ) std::cout << __LINE__ << std::endl;
      
  /// close files

  if ( _ftest ) _ftest->Close();
  if ( _fref )  _fref->Close();
      

  /// now actually overwrite the old reference file

  if ( deleteref ) { 
    std::cout << "ref " << frefname << "\ttest " << ftestname << std::endl; 
    if ( frefname !=  ftestname ) { 
      std::string cmd = std::string("mv ") + frefname + " " + frefname + ".bak";
      std::system( cmd.c_str() );
      
      cmd = std::string("mv .newout.root ") + std::string(frefname);
      std::system( cmd.c_str() );
    }
    else {
      std::cerr << "reference file and test file are the same - not replacing" << std::endl;
    }
  }  
 
  if ( _ftest ) delete _ftest;
  if ( _fref )  delete _fref;

  return 0;
}
