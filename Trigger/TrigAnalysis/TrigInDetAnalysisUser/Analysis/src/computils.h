/* emacs: this is -*- c++ -*- */
/**
 **     @file    computils.h
 **
 **     @author  mark sutton
 **     @date    Sat Aug 30 2014 14:38:03 CEST  
 **
 **     Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
 **/

#ifndef COMPUTILS_H
#define COMPUTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <exception>

#include "label.h"
#include "utils.h"
#include "DrawLabel.h" 


#include "TStyle.h"
#include "TPad.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

#include "TLegend.h"
#include "TrigInDetAnalysis/Efficiency.h"


extern bool LINEF;
extern bool LINES;


void ATLASFORAPP_LABEL( double x, double y, int color, double size=0.06 ); 

void myText( Double_t x, Double_t y, Color_t color, const std::string& text, Double_t tsize);


/// return the current data and time
std::string stime();
static std::string release;

double integral( TH1* h );

void Norm( TH1* h, double scale=1 );

double Entries( TH1* h );

/// does a string contain the substring
bool contains( const std::string& s, const std::string& p);
bool contains( const std::string& s, char p) noexcept;

/// does a string contain the substring at the beginning of the string
bool fcontains( const std::string& s, const std::string& p);

/// does a file exist
bool exists( const std::string& filename );

/// tail of a string 
std::string tail( std::string s, const std::string& pattern );

/// head of a string 
std::string head( std::string s, const std::string& pattern );

/// match a file name
std::string globbed( const std::string& s );

void contents( std::vector<std::string>& keys, 
	       TDirectory* td, 
	       const std::string& directory="", 
	       const std::string& pattern="", 
	       const std::string& path="" );

void contents( std::vector<std::string>& keys, 
	       TDirectory* td, 
	       const std::string& directory="", 
	       const std::vector<std::string>& patterns=std::vector<std::string>(), 
	       const std::string& path="" );


double realmax( TH1* h, bool include_error=true, double lo=0, double hi=0 );
double realmin( TH1* h, bool include_error=true, double lo=0, double hi=0 );

std::string findcell( std::string name, const std::string& regex, const std::string& splitex="/" );

std::string findrun( TFile *f );

double plotable( TH1* h ); // , double xlo=-999, double xhi=-999 );


class data_mismatch : public std::exception { 
public:
  data_mismatch(const std::string& s) { std::cerr << "exception:" << what() << " " << s << std::endl; }; 
  virtual const char* what() const throw() { return "data don't match"; }
};



template<typename T>
std::ostream& operator<<( std::ostream& s, std::vector<T>& v) { 
  for ( unsigned i=0 ; i<v.size() ; i++ ) s << "\t" << v[i];
  return s;
}


/// automatically set the xrange on a histogram
std::vector<int>    findxrange(TH1* h, bool symmetric=false );
std::vector<double> findxrangeuser(TH1* h, bool symmetric=false );



void trim_tgraph( TH1* h, TGraphAsymmErrors* t );  


void xrange(TH1* h, bool symmetric=true );

/// copy the TTree of release info from one directory to another
void copyReleaseInfo( TFile* finput, TFile* foutdir );


class true_mean { 

public: 
  
  true_mean( TH1F* h );

  double mean() const  { return m_mean; } 
  double error() const { return m_error; } 

private:
  
  double m_mean;
  double m_error;
  
};


/// class to store information about axes, limits, whether it is 
/// log or linear scale etc

class AxisInfo { 

public:

  AxisInfo( const std::string& s ) : 
    m_info(s), 
    m_log(false),
    m_autoset(false),
    m_symmetric(false),
    m_rangeset(false),
    m_lo(0),
    m_hi(0),
    m_norm(false),
    m_refnorm(false),
    m_binwidth(false),
    m_offset(0),
    m_trim(false)
  { 
    //    std::cout << "AxisInfo::info" << m_info << std::endl;

    std::vector<std::string> keys = split( s, ":" );
    
    //    std::cout << "\n\n" << s << "\nnkeys " << keys.size() << std::endl; 

    if ( keys.size()>0 ) m_tag = keys[0];
    
    bool minset = false;
    //  bool maxset = false;
    
    for ( size_t i=1 ; i<keys.size() ; i++ ) { 
      
      if       ( keys[i]=="lin" )   m_log       = false;
      else if  ( keys[i]=="log" )   m_log       = true;
      else if  ( keys[i]=="sym" )   m_symmetric = true; 
      else if  ( keys[i]=="norm" )  m_norm      = true;
      else if  ( keys[i]=="refn" )  m_refnorm   = true;
      else if  ( keys[i]=="width" ) m_binwidth  = true;
      else if  ( keys[i]=="auto" )  m_autoset   = true;
      else if  ( keys[i]=="trim" )  m_trim      = true;
      else if  ( keys[i].find("offset")==0  )  { 

	std::cout << "offset:" << std::endl;
	std::cout << "\tkey: " << keys[i] << std::endl;
	std::cout << "\tpos: " << keys[i].find("offset") << std::endl;

	std::string offset = keys[i];
	m_offset=std::atof(offset.substr(6,offset.size()-6).c_str());

	std::cout << "m_offset: " << m_offset << std::endl;
      }
      else if  ( keys[i]=="auton" )  {
	m_autoset = true;
	m_norm    = true;
      }
      else if  ( keys[i]=="autow" )  {
	m_autoset  = true;
	m_binwidth = true;
      }
      else if  ( keys[i]=="autown" || keys[i]=="autonw" ) {
	m_autoset  = true;
	m_norm     = true;
	m_binwidth = true;
      }
      else if  ( keys[i]=="autosym" ) { 
	m_autoset = true; 
	m_symmetric = true; 
      }
      else if  ( keys[i]=="normw" ||  keys[i]=="widthn" )  {
	m_norm     = true;
	m_binwidth = true;
      }
      else if  ( !minset )  { 
	m_lo = std::atof(keys[i].c_str());
	i++;
	if ( i<keys.size() ) m_hi = std::atof(keys[i].c_str());
	else {
	  std::cerr << "not enough values for the axis range: " << s << std::endl;
	  std::exit(-1);
	}
	minset = true;
	// maxset = true;
	m_rangeset = true;
      }
    }
        
#if 0
    std::cout << "AxisInfo::info" << m_info << "\n";
    std::cout << "\tlog   " << m_log       << "\n";
    std::cout << "\tauto  " << m_autoset   << "\n";
    std::cout << "\tsym   " << m_symmetric << "\n";
    std::cout << "\trange " << m_rangeset << " : " << m_lo << " - " << m_hi << std::endl;
#endif

  }
  
  /// accessors 

  std::string tag() const { return m_tag; }

  bool   log()  const { return m_log; }

  bool   autoset() const { return m_autoset; }

  bool   normset()    const { return m_norm; }

  bool   refnormset() const { return m_refnorm; }
  
  bool   symmetric() const { return m_symmetric; }

  bool   rangeset() const { return m_rangeset; }

  bool   trim() const { return m_trim; }


  double offset() const { return m_offset; }


  double lo() const { return m_lo; } 
  double hi() const { return m_hi; } 
  
  double binwidth() const { return m_binwidth; }

  std::string c_str() const { return m_info; }


public:

  static std::vector<std::string> split( const std::string& s, const std::string& t=":"  ) {
    
    std::string sc = s;
    size_t pos = sc.find(t);
    
    std::vector<std::string> tags;
    
    while ( pos!=std::string::npos ) { 
      tags.push_back( chop(sc,t) );
      pos = sc.find(t);
    }
    
    tags.push_back(sc);
    
    return tags;
  } 
  

public:

  std::string m_info;
  
  std::string m_tag;

  bool   m_log;
  bool   m_autoset;
  bool   m_symmetric;

  bool   m_rangeset;
  double m_lo; 
  double m_hi;

  bool   m_norm;
  bool   m_refnorm;
  
  bool   m_binwidth;
  
  double m_offset;

  bool   m_trim;

};


inline std::ostream& operator<<( std::ostream& s, const AxisInfo& a ) { 
  s << "[ " << a.tag() << ( a.log() ? " : log" : "" ) << " ";
  if      (  a.autoset() ) s << " : auto";
  else if ( a.rangeset() ) s << " : range " << a.lo() << " - " << a.hi();
  s << " ]";
  return s; 
}







/// slightly more convenient legend class
class Legend { 

public:

  Legend() : m_leg(nullptr){ 
    m_x[0]=0.0;
    m_y[0]=0.0;
    m_x[1]=0.0;
    m_y[1]=0.0;
	  } 

  Legend(double x1, double x2, double y1, double y2): m_leg(nullptr) { 
    m_x[0]=x1;
    m_y[0]=y1;
    m_x[1]=x2;
    m_y[1]=y2;
  }


  // Legend( const Legend& leg ) : mleg((TLegend*)leg.mleg->Clone()) { } 

 Legend(const Legend& legend) : m_leg(legend.m_leg) { } 

  ~Legend() { } 

  TLegend* legend() { return m_leg; } 

  size_t size() const { 
    return m_entries.size();
  }  

  double TextSize() const { return m_leg->GetTextSize(); }

  int TextFont() const { return m_leg->GetTextFont(); }

  double height() const { return m_y[1]-m_y[0]; }

  double width() const { return m_x[1]-m_x[0]; }

  void AddEntry( TObject* tobj, const std::string& s, const std::string& type="p" ) { 
    m_obj.push_back( tobj );
    m_entries.push_back( s );
    m_type.push_back( type );
  }

  void Draw() { 
    
    /// ha ! don't actually create the legend until we want to draw it, 
    /// then we can determine the size etc automatically

    double y0 = 0;

    if ( m_y[0]>0.5 ) { 
      y0 = m_y[1] - m_entries.size()*0.05;
    }
    else { 
      y0 = m_y[0];
      m_y[1] = y0 + m_entries.size()*0.05;
    }
      
    m_leg = new TLegend( m_x[0], y0, m_x[1], m_y[1] );
     
    m_leg->SetBorderSize(0);
    m_leg->SetTextFont(42);
    m_leg->SetTextSize(0.04);
    m_leg->SetFillStyle(3000);
    m_leg->SetFillColor(0);
    m_leg->SetLineColor(0);
    
    for ( size_t i=0 ; i<m_entries.size() ; i++ ) { 
      m_leg->AddEntry( m_obj[i], m_entries[i].c_str(), m_type[i].c_str() );
    }

    m_leg->Draw();
  }


private:

  TLegend* m_leg;

  double m_x[2];
  double m_y[2];

  std::vector<TObject*>    m_obj;
  std::vector<std::string> m_entries;
  std::vector<std::string> m_type;

}; 


extern int   colours[6]; // = {  1,  2, kBlue-4,  6, kCyan-2,  kMagenta+2 };
extern int   markers[6]; // = { 20, 24,      25, 26,      25,          22 };
extern double msizes[6]; // = {  1,  1,       1,  1,       1,            1 };


template<typename T>
void setParameters( T* h, TGraphAsymmErrors* tg ) { 
  tg->SetLineColor( h->GetLineColor() );
  tg->SetMarkerColor( h->GetMarkerColor() );
  tg->SetMarkerStyle( h->GetMarkerStyle() );
  tg->SetMarkerSize( h->GetMarkerSize() );
  tg->SetLineWidth( h->GetLineWidth() );
  tg->SetLineStyle( h->GetLineStyle() );
}


template<typename T>
void zeroErrors( T* h ) {
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) h->SetBinError( i, 1e-100 ); 
} 



/// generic plotter class - better to have one of these - make 
/// sure it can be configured however you like, line styles, 
/// marker types, legends etc 
/// now a template so can be used for TH1D and TH2D etc
template<typename T>
class tPlotter { 

public: 
  
  tPlotter(T* _htest=0, T* _href=0, const std::string& s="", TGraphAsymmErrors* _tgtest=0, TGraphAsymmErrors* _tgref=0 ) : 
    m_htest(_htest), m_href(_href),
    m_tgtest(_tgtest), m_tgref(_tgref),
    m_plotfilename(s),
    m_max_entries(4),
    m_entries(0), 
    m_trim_errors(false)
  {
    // plotref = true;
  }

  
  tPlotter(const tPlotter& p) : 
    m_htest(p.m_htest),   m_href(p.m_href), 
    m_tgtest(p.m_tgtest), m_tgref(p.m_tgref), 
    m_plotfilename(p.m_plotfilename),
    m_max_entries(p.m_max_entries),
    m_entries(0),
    m_trim_errors(p.m_trim_errors) {     
  }

  

  /// sadly, root manages all the histograms (does it really? 
  /// who can tell) so we mustn't delete anything just in case
  // ~tPlotter() { if ( m_tg ) delete m_tg; }
  /// NO, NO, NO, don't delete the objects, root need them because 
  /// it is TOO STUPID to allow objects to be used as actual objects
  ~tPlotter() { } 

  std::string plotfilename() const { return m_plotfilename; }

  void trim_errors(bool b) { m_trim_errors=b; } 
  
  bool  trim_errors() const { return m_trim_errors; } 
  
  void Draw( int i, Legend* lleg, bool mean=false, bool first=true, bool drawlegend=false ) { 
    
    if ( htest() ) {
      gStyle->SetOptStat(0);
      if ( href() ) { 
	href()->SetLineColor(colours[i%6]);
	href()->SetLineStyle(2);
	href()->SetMarkerStyle(0);
      }
      
      if ( LINEF ) htest()->SetLineColor(colours[i%6]);

      htest()->SetLineStyle(1);

      if ( LINEF ) htest()->SetMarkerColor(htest()->GetLineColor());

      if ( LINEF ) htest()->SetMarkerStyle(markers[i%6]);
      // if ( LINEF && htest()->GetMarkerStyle() == 20 ) 

      if ( LINEF ) htest()->SetMarkerSize( msizes[i%6]*htest()->GetMarkerSize() );

      if ( htest() ) std::cout << "\tentries: " << plotable( htest() );

      std::cout << std::endl;

      if ( first )  {

	if ( tgtest() ) { 
	    zeroErrors(htest());
	    htest()->GetXaxis()->SetMoreLogLabels(true);
	    if ( trim_errors() ) trim_tgraph( htest(), tgtest() );


	    htest()->Draw("ep");
	    if ( LINES ) htest()->Draw("lhistsame");
	    setParameters( htest(), tgtest() );
	    tgtest()->Draw("esame");
	}
	else { 
	    htest()->GetXaxis()->SetMoreLogLabels(true);
	    htest()->Draw("ep");
	    if ( LINES ) htest()->Draw("lhistsame");
	}
       
      }

      if ( plotref && href() ) { 
	if ( contains(href()->GetName(),"_vs_")  || 
	     contains(href()->GetName(),"sigma") || 
	     contains(href()->GetName(),"mean") || 
	     contains(href()->GetName(),"_eff")  ||
	     contains(href()->GetName(),"Res_") || 
	     contains(href()->GetName(),"Eff_") ) href()->Draw("hist same][");
	else                                      href()->Draw("hist same");
      }

      if ( tgtest() ) { 
	zeroErrors(htest());

	if ( trim_errors() ) trim_tgraph( htest(), tgtest() );
	setParameters( htest(), tgtest() );
	tgtest()->Draw("e1same");
	if ( LINES ) tgtest()->Draw("lsame");

      }
      
#if 0
      /// solid white centres for the marker types
      if ( htest()->GetMarkerStyle()>23 ) { 
	TH1D* hnull = (TH1D*)htest()->Clone("duff"); hnull->SetDirectory(0);
	zeroErrors( hnull );
	hnull->SetLineColor(kWhite);
	hnull->SetMarkerStyle( htest()->GetMarkerStyle()-4 );
	//	hnull->SetMarkerStyle( 0 ); 
	hnull->SetMarkerColor(kWhite);
	hnull->SetMarkerSize( htest()->GetMarkerSize()*0.75 );
	// hnull->SetMarkerSize( 0 );
	hnull->DrawCopy("l same");
	delete hnull;
      }
#endif

      htest()->Draw("ep same");
      if ( LINES ) htest()->Draw("lhist same");

      // href()->Draw("lhistsame");
      // htest()->Draw("lhistsame");

      std::string key = m_plotfilename;

      static TH1D* hnull = new TH1D("hnull", "", 1, 0, 1);
      hnull->SetMarkerColor(kWhite);
      hnull->SetLineColor(kWhite);
      hnull->SetMarkerStyle(0);
      hnull->SetLineStyle(0);
      hnull->SetLineWidth(0);
      hnull->SetMarkerSize(0);

      
      if ( lleg ) { 
      
	Legend& leg = *lleg;
	
	if ( mean ) { 
	
	  char meanrefc[64];
	  bool displayref = false;
	  if ( meanplotref && href() ) { 
	    displayref = true;
	    true_mean muref( href() );
	    std::sprintf( meanrefc, " <t> = %3.2f #pm %3.2f ms (ref)", muref.mean(), muref.error() );
	  }
	  else { 
	    std::sprintf( meanrefc, "%s", "" );
	  }
	  
	  true_mean mutest( htest() );
	  char meanc[64];
	  std::sprintf( meanc, " <t> = %3.2f #pm %3.2f ms", mutest.mean(), mutest.error() );
	  
	  std::string dkey = key;
	  
	  std::string remove[7] = { "TIME_", "Time_", "All_", "Algorithm_", "Class_", "HLT_", "Chain_HLT_" };
	  
	  if ( dkey.find("Chain")!=std::string::npos ) {
	    if ( dkey.find("__")!=std::string::npos ) dkey.erase( 0, dkey.find("__")+2 );
	  } 
	  
	  
	  for ( int ir=0 ; ir<7 ; ir++ ) { 
	    if ( dkey.find( remove[ir] )!=std::string::npos ) dkey.erase( dkey.find( remove[ir]), remove[ir].size() );
	  } 
	  
	  std::string rkey = dkey;
	  
	  
	  if ( LINEF || leg.size() < m_max_entries ) { 
	    dkey += std::string(" : ");
	    
	    if ( (dkey+meanc).size()>58 ) { 
	      leg.AddEntry( htest(), dkey, "p" );
	      leg.AddEntry( hnull,   meanc,        "p" ); 
	    }
	    else { 
	      leg.AddEntry( htest(), (dkey+meanc).c_str(), "p" );
	    }	  
	    
	    if ( displayref ) { 
	      rkey += std::string(" : ");
	      /// not quite yet ...
	      //  leg.AddEntry( hnull, "", "l" );
	      
	      if ( (rkey+meanrefc).size()>58 ) { 
		leg.AddEntry( href(), rkey, "l" );
		leg.AddEntry( hnull,  meanrefc,     "l" ); 
	      }
	      else {
		leg.AddEntry( href(), (rkey+meanrefc).c_str(), "l" );
	      }
	    }
	  }
	  
	}
	else { 
	  if ( LINEF || leg.size()<m_max_entries ) leg.AddEntry( htest(), key, "p" );
	}
	
	m_entries++;
	
	if ( drawlegend ) leg.Draw();
      }
      
    }
  }



  void DrawLegend( int i, Legend& leg, bool mean=false, bool first=true, bool drawlegend=false ) { 
    
    if ( htest() ) {
      gStyle->SetOptStat(0);
      if ( href() ) { 
	href()->SetLineColor(colours[i%6]);
	href()->SetLineStyle(2);
	href()->SetMarkerStyle(0);
      }

      if ( LINEF ) htest()->SetLineColor(colours[i%6]);
      htest()->SetLineStyle(1);
      if ( LINEF ) htest()->SetMarkerColor(htest()->GetLineColor());
      if ( LINEF ) htest()->SetMarkerStyle(markers[i%6]);

      if ( htest() ) std::cout << "\tentries: " << plotable( htest() );
      std::cout << std::endl;

      if ( first )  {

	if ( tgtest() ) { 
	    zeroErrors(htest());
	    htest()->GetXaxis()->SetMoreLogLabels(true);
	    if ( trim_errors() ) trim_tgraph( htest(), tgtest() );


	    // htest()->Draw("ep");
	    if ( LINES ) htest()->Draw("lhistsame");
	    setParameters( htest(), tgtest() );
	    //   tgtest()->Draw("esame");
	}
	else { 
	    htest()->GetXaxis()->SetMoreLogLabels(true);
	    //    htest()->Draw("ep");
	    //  if ( LINES ) htest()->Draw("lhistsame");
	}
       
      }


#if 0
      if ( plotref && href() ) { 
	if ( contains(href()->GetName(),"_vs_")  || 
	     contains(href()->GetName(),"sigma") || 
	     contains(href()->GetName(),"mean") || 
	     contains(href()->GetName(),"_eff")  ||
	     contains(href()->GetName(),"Res_") || 
	     contains(href()->GetName(),"Eff_") ) href()->Draw("hist same][");
	else                                      href()->Draw("hist same");
      }

      if ( tgtest() ) { 
	zeroErrors(htest());

	if ( trim_errors() ) trim_tgraph( htest(), tgtest() );
	setParameters( htest(), tgtest() );
	tgtest()->Draw("e1same");
	if ( LINES ) tgtest()->Draw("lsame");

      }
#endif

      
#if 0
      /// solid white centres for the marker types
      if ( htest()->GetMarkerStyle()>23 ) { 
	TH1D* hnull = (TH1D*)htest()->Clone("duff"); hnull->SetDirectory(0);
	zeroErrors( hnull );
	hnull->SetLineColor(kWhite);
	hnull->SetMarkerStyle( htest()->GetMarkerStyle()-4 );
	//	hnull->SetMarkerStyle( 0 ); 
	hnull->SetMarkerColor(kWhite);
	hnull->SetMarkerSize( htest()->GetMarkerSize()*0.75 );
	// hnull->SetMarkerSize( 0 );
	hnull->DrawCopy("l same");
	delete hnull;
      }
#endif

      //      htest()->Draw("ep same");
      //      if ( LINES ) htest()->Draw("lhist same");

      // href()->Draw("lhistsame");
      // htest()->Draw("lhistsame");

      std::string key = m_plotfilename;

      static TH1D* hnull = new TH1D("hnull", "", 1, 0, 1);
      hnull->SetMarkerColor(kWhite);
      hnull->SetLineColor(kWhite);
      hnull->SetMarkerStyle(0);
      hnull->SetLineStyle(0);
      hnull->SetLineWidth(0);
      hnull->SetMarkerSize(0);

      
      if ( mean ) { 

	char meanrefc[64];
	bool displayref = false;
	if ( meanplotref && href() ) { 
	  displayref = true;
	  true_mean muref( href() );
	  std::sprintf( meanrefc, " <t> = %3.2f #pm %3.2f ms (ref)", muref.mean(), muref.error() );
	}
	else { 
	  std::sprintf( meanrefc, "%s", "" );
	}


	true_mean mutest( htest() );
	char meanc[64];
	std::sprintf( meanc, " <t> = %3.2f #pm %3.2f ms", mutest.mean(), mutest.error() );
	
	std::string dkey = key;

	std::string remove[7] = { "TIME_", "Time_", "All_", "Algorithm_", "Class_", "HLT_", "Chain_HLT_" };

	if ( dkey.find("Chain")!=std::string::npos ) {
	  if ( dkey.find("__")!=std::string::npos ) dkey.erase( 0, dkey.find("__")+2 );
	} 
	  

	for ( int ir=0 ; ir<7 ; ir++ ) { 
	  if ( dkey.find( remove[ir] )!=std::string::npos ) dkey.erase( dkey.find( remove[ir]), remove[ir].size() );
	} 

	std::string rkey = dkey;


	if ( LINEF || leg.size() < m_max_entries ) { 
	  dkey += std::string(" : ");

	  if ( (dkey+meanc).size()>58 ) { 
	    leg.AddEntry( htest(), dkey, "p" );
	    leg.AddEntry( hnull,   meanc,        "p" ); 
	  }
	  else { 
	    leg.AddEntry( htest(), (dkey+meanc).c_str(), "p" );
	  }	  

	  if ( displayref ) { 
	    rkey += std::string(" : ");
	    //  leg.AddEntry( hnull, "", "l" );

	    if ( (rkey+meanrefc).size()>58 ) { 
	      leg.AddEntry( href(), rkey, "l" );
	      leg.AddEntry( hnull,  meanrefc,     "l" ); 
	    }
	    else {
	      leg.AddEntry( href(), (rkey+meanrefc).c_str(), "l" );
	    }
	  }
	}

      }
      else { 
	if ( LINEF || leg.size()<m_max_entries ) leg.AddEntry( htest(), key, "p" );
      }

      m_entries++;

      if ( drawlegend ) leg.Draw();

    }
  }

  

  /// print the output 
  void Print(const std::string& s="") {
    if ( s!="" ) gPad->Print(s.c_str());
    else         gPad->Print(m_plotfilename.c_str());
  } 

  T* htest() { return m_htest; }
  T* href()  { return m_href; }


  TGraphAsymmErrors* tgtest() { return m_tgtest; }
  TGraphAsymmErrors* tgref()  { return m_tgref; }


  void max_entries( int i ) { m_max_entries = i; } 

public:

  static void setplotref( bool b )     { plotref=meanplotref=b; }
  static void setmeanplotref( bool b ) { meanplotref=b; }

private:

  /// actual histograms 
  T* m_htest;
  T* m_href;

  TGraphAsymmErrors* m_tgtest;
  TGraphAsymmErrors* m_tgref;

  std::string m_plotfilename;

  static bool plotref;
  static bool meanplotref;

  size_t  m_max_entries;
  size_t  m_entries;

  bool m_trim_errors;

};

template<typename T>
bool tPlotter<T>::plotref = true;


template<typename T>
bool tPlotter<T>::meanplotref = true;




typedef tPlotter<TH1F> Plotter;

bool empty( TH1* h );



inline void hminus(TH1* h) { 
  std::cout << __FUNCTION__ << std::endl; 
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    double duff = h->GetBinContent(i);
    if (  duff<0 ) { 
      std::cout<< "\t\t" << __FUNCTION__ << " " << h->GetName() << "  " << i << " " << h->GetBinContent(i) << " " << (duff*1e6) << std::endl;  
    }
  }
  h->DrawCopy();
  gPad->Print( (std::string("duff-")+h->GetName()+".pdf").c_str() );
}


/// set of generic plots
class Plots : public std::vector<Plotter> { 
  
public:
  
  Plots(const std::string& s="", bool errors=false ) : 
    m_name(s), 
    m_logx(false), m_logy(false), 
    m_maxset(false), m_max(0),
    m_minset(false), m_min(0),
    m_rangeset(false), 
    m_lo(0.0), m_hi(0.0),
    m_trim_errors(errors)
  { }

  double realmin( double lo=0, double hi=0 ) {
    bool first = true;
    double _min = 1000;
    for ( unsigned i=0 ; i<size() ; i++ ) {
      double rmtest = ::realmin( at(i).htest(), false, lo, hi );
      if ( rmtest!=0 && ( first || _min>rmtest ) ) _min = rmtest;
      if ( rmtest!=0 ) first = false;
    }
    return _min;
  }

  double realmax(double lo=0, double hi=0) {
    bool first = true;
    double _max = 0;
    for ( unsigned i=0 ; i<size() ; i++ ) {
      //      double rmref  = realmin( at(i).href(), false );
      double rmtest = ::realmax( at(i).htest(), false, lo, hi );
      if ( rmtest!=0 && ( first || _max<rmtest ) ) _max = rmtest;
      if ( rmtest!=0 ) first = false;
    }
    return _max;
  }


  
  void MaxScale(double scale=1.1, double lo=0, double hi=0) {  
    
    if ( size()<1 ) return;

    double tmax = realmax(lo,hi);
    double tmin = realmin(lo,hi);
   
    m_max = scale*tmin;
    
    if ( m_logy ) m_min = tmin;

    for ( unsigned i=0 ; i<size() ; i++ ) {
      if ( at(i).href() ) at(i).href()->SetMaximum(scale*tmax);
      at(i).htest()->SetMaximum(scale*tmax);
    } 

  }


  void MinScale( double scale, double lo=0, double hi=0 ) { 

    if ( size()<1 ) return;

    if ( scale==0 ) { 
      for ( unsigned i=0 ; i<size() ; i++ ) {
	if ( at(i).href() ) at(i).href()->SetMinimum(0);
	at(i).htest()->SetMinimum(0);
      }  
      return;
    }

    
    double tmin = 0;

    if ( lo!=hi ) tmin = realmin( lo, hi );
    else          tmin = realmin();

    m_min = scale*tmin;

    for ( unsigned i=0 ; i<size() ; i++ ) {
      if ( at(i).href() ) at(i).href()->SetMinimum(scale*tmin);
      at(i).htest()->SetMinimum(scale*tmin);
    }
  }
 

  void Min( double scale ) { 
    m_minset = true;
    for ( unsigned i=0 ; i<size() ; i++ ) {
      if ( at(i).href() ) at(i).href()->SetMinimum(scale);
      at(i).htest()->SetMinimum(scale);
      if ( m_logy ) { 
	if ( at(i).href() ) if ( at(i).href()->GetMinimum()<=0 )  at(i).href()->GetMinimum(1e-4);
	if ( at(i).htest()->GetMinimum()<=0 ) at(i).htest()->GetMinimum(1e-4);
      }
    }
  }
  

  void Max( double scale ) { 
    m_maxset = true;
    for ( unsigned i=0 ; i<size() ; i++ ) {
      if ( at(i).href() ) at(i).href()->SetMaximum(scale);
      at(i).htest()->SetMaximum(scale);
    }
  }

  std::vector<double> findxrange( bool symmetric=false ) { 

    /// don't use empty histograms to find the bin limits
    /// unless they are *all* empty 
    
    std::vector<double> v(2,0);
    
    TH1F* hf = at(0).htest();

    double vlo  =  1e21;
    double vhi  = -1e21;
    
    if ( hf->GetBinLowEdge(1)<vlo )                 vlo = hf->GetBinLowEdge(1);
    if ( hf->GetBinLowEdge(hf->GetNbinsX()+1)>vhi ) vhi = hf->GetBinLowEdge( hf->GetNbinsX()+1 );
      
    if ( size()>0 ) v = ::findxrangeuser( hf, symmetric );

    bool first = false;

    for ( unsigned i=1 ; i<size() ; i++ ) { 
  
      hf = at(i).htest();

      if ( ::empty( hf ) ) continue;

      if ( hf->GetBinLowEdge(1)<vlo )                 vlo = hf->GetBinLowEdge(1);
      if ( hf->GetBinLowEdge(hf->GetNbinsX()+1)>vhi ) vhi = hf->GetBinLowEdge( hf->GetNbinsX()+1 );

      std::vector<double> limits = ::findxrangeuser( hf, symmetric );

      double lo = limits[0];
      double hi = limits[1];

      if ( first ) { 
	v[0] = lo;
	v[1] = hi;
      }
      else { 
	if ( v[0]>lo ) v[0] = lo;
	if ( v[1]<hi ) v[1] = hi;
      }

      first = false;
    }

    double upper = ( v[1]-v[0] )*1.1 + v[0];
    double lower = v[0] - ( v[1]-v[0] )*0.1; 

    if ( m_logx ) {
      double dx = std::log10(v[1])-std::log10(v[0]);
      upper = std::pow(10,dx*1.1 + std::log10(v[0]));
      lower = std::pow(10,std::log10(v[0]) - dx*0.1); 
    }

    if ( lower<vlo ) lower = vlo;
    if ( upper>vhi ) upper = vhi;

    v[0] = lower;
    v[1] = upper;
    
    return v;
  }
  


  void sortx( const AxisInfo& xinfo ) {
    
    if ( xinfo.rangeset() ) { 
      m_lo = xinfo.lo();
      m_hi = xinfo.hi();
    }
    
    if ( xinfo.autoset() && size() > 0 ) {
      std::vector<double> limits = findxrange( xinfo.symmetric() );
      if ( xinfo.rangeset() ) { 
	if ( limits[0]<m_lo ) m_lo = limits[0];
	if ( limits[1]>m_hi ) m_hi = limits[1];
      }
      else { 
	m_lo = limits[0];
	m_hi = limits[1];
      }
    }
    else if (size() == 0)
    {
        std::cout<<"Warning in computils.h::sortx() size=0.  Setting m_lo/m_hi to 0/1.  You will probably have empty figures."<<std::endl;
        m_lo = 0;
        m_hi = 1;
    }

    if ( xinfo.rangeset() || xinfo.autoset() ) { 
      SetRangeUser( m_lo, m_hi );
      if ( xinfo.log() && m_lo>0 ) SetLogx(true);
      else                         SetLogx(false);
    }

  }


  double lo() const { return m_lo; }
  double hi() const { return m_hi; }


  void xrange(bool symmetric=false) { 
    m_rangeset = false;
    for ( unsigned i=0 ; i<size() ; i++ ) { 
      if ( at(i).href() ) ::xrange( at(i).href(), symmetric );
      ::xrange( at(i).htest(), symmetric );
    }
  }

  void SetRangeUser( double lo, double hi ){ 
    m_rangeset = true;
    m_lo = lo;
    m_hi = hi;
    for ( unsigned i=0 ; i<size() ; i++ ) { 
      if ( at(i).href() ) at(i).href()->GetXaxis()->SetRangeUser( m_lo, m_hi );
      at(i).htest()->GetXaxis()->SetRangeUser( m_lo, m_hi );
    }
  }
  

  void limits() { 
    double rmax = realmax();
    double rmin = realmin();    
    if ( rmin<0 ) { 
      std::cout << "\tlimits \t" << m_name << "\tmin " << rmin << "\tmax " << rmax << std::endl; 
    }
  }


   
  void Draw( Legend& leg, bool means=false ) {  
    Draw_i( leg, means );
    if ( LINES ) { 
      LINES=false;
      Draw_i( leg, means );
      LINES=true;
    }
  }

  void Draw_i( Legend& leg, bool means=false ) {  

    bool first = true;
    
    if ( m_logy ) {      /// increase the number of log labels if only a few decades
      for ( unsigned i=0 ; i<size() ; i++, first=false ) { 
	double ymax = at(i).htest()->GetMaximum();
	double ymin = at(i).htest()->GetMinimum();
	at(i).htest()->GetYaxis()->SetMoreLogLabels(true);
	if ( ymax/ymin>1e6 ) at(i).htest()->GetYaxis()->SetMoreLogLabels(false);
	break;
      }
    }

    for ( unsigned i=0 ; i<size() ; i++ ) at(i).trim_errors( m_trim_errors );

    /// still can't get this working correctly - for the time being leave these alternative
    /// loop orders in place but commented out until we can find a better solution ...
  
    //  for ( unsigned i=0 ; i<size() ; i++,  first=false ) at(i).DrawLegend( i, leg, means, first, (i==size()-1) );
    //  for ( unsigned i=0 ; i<size() ; i++,  first=false ) at(i).Draw( i, &leg, means, first, (i==size()-1) );
    for ( unsigned i=size() ; i-- ;  first=false ) at(i).Draw( i, &leg, means, first, i==0 );

    if ( watermark ) DrawLabel(0.1, 0.02, "built "+stime()+release, kBlack, 0.03 );

    gPad->SetLogy(m_logy);
    gPad->SetLogx(m_logx);
  }

  void SetLogx( bool b=true ) { m_logx=b; } 
  void SetLogy( bool b=true ) { m_logy=b; } 

  std::string GetXaxisTitle() { 
    if ( size()>0 ) return at(0).htest()->GetXaxis()->GetTitle();
    return "";
  }

  std::string GetYaxisTitle() { 
    if ( size()>0 ) return at(0).htest()->GetYaxis()->GetTitle();
    return "";
  }

  void SetXaxisTitle(const std::string& s) { 
    if ( size()>0 ) at(0).htest()->GetXaxis()->SetTitle(s.c_str());
  }

  void SetYaxisTitle(const std::string& s) { 
    if ( size()>0 ) at(0).htest()->GetYaxis()->SetTitle(s.c_str());
  }


  /// this is so that we can update the stats as we go along, 
  /// but no updating isn't done at the moment 
  void push_back(const Plotter& val) {
    std::vector<Plotter>::push_back( val );
  }

public:

  static void setwatermark(bool b) { watermark = b; }

private:
  
  std::string m_name;

  /// canvas log setting
  bool m_logx;
  bool m_logy;
  
  /// yaxis range setting
  bool   m_maxset;
  double m_max;
  
  bool   m_minset;
  double m_min;

  /// xaxis range setting
  bool   m_rangeset;
  double m_lo;
  double m_hi;

  bool m_trim_errors;

private:

  static bool watermark;

};




/// details of the histogram axes etc

class HistDetails { 

public:

  HistDetails( const std::vector<std::string>& v ) : m_extra(""), m_xinfo(v[2]), m_yinfo(v[4]) { 
    if ( v.size() < 6 ) throw std::exception();
    m_details.reserve(6);
    for ( size_t i=0 ; i<6 ; i++ ) m_details.push_back(v[i]); 
    getextra();
  }

  HistDetails( const std::string* vp ) : m_extra(""), m_xinfo(vp[2]), m_yinfo(vp[4]) { 
    m_details.reserve(6);
    for ( size_t i=0 ; i<6 ; i++ ) m_details.push_back(vp[i]); 
    getextra();
  }

  std::string  name() const { return m_details[0]; } 

  std::string  detail() const { return m_extra; }

  std::string  info() const { return m_details[1]; } 

  std::string xtitle() const { return m_details[3]; }
  std::string ytitle() const { return m_details[5]; }

  const AxisInfo& xaxis() const { return m_xinfo; }
  const AxisInfo& yaxis() const { return m_yinfo; }

private:

  void getextra() { 
    if ( contains( m_details[0], "-" ) ) { 
      m_extra = m_details[0].substr( m_details[0].find('-'), m_details[0].size() );
      m_details[0] = m_details[0].substr( 0, m_details[0].find('-') );
    }
    if ( contains( m_details[0], "+" ) ) { 
      m_extra = m_details[0].substr( m_details[0].find('+'), m_details[0].size() );
      m_details[0] = m_details[0].substr( 0, m_details[0].find('+') );
    }
  }

private:

  std::vector<std::string> m_details;

  std::string m_extra;

  AxisInfo m_xinfo;
  AxisInfo m_yinfo;

};



inline std::ostream& operator<<( std::ostream& s, const HistDetails& h ) { 
  return s << "[ " << h.name() << "  \tx: \"" << h.xtitle() << "\"   " << h.xaxis() << "\t\ty: \"" << h.ytitle() << "\"   " << h.yaxis() << " ]"; 
}



// plot panel inforamtion

class Panel { 

public:

  /// don't know how many rows or total hists yet,
  /// but do know number of columns 
  Panel( const std::string& s, int nc ) : 
    m_name(s), m_nhist(-1), m_nrows(-1), m_ncols(nc) { 
    m_hist.reserve( nc );
  } 

  /// know number of rows and columns
  Panel( const std::string& s, int nr, int nc ) : 
    m_name(s), m_nhist(nr*nc), m_nrows(nr), m_ncols(nc) { 
    m_hist.reserve( m_nhist );
  } 

  void push_back( const HistDetails& h ) { 
    m_hist.push_back( h );
    m_nhist = m_hist.size();
    m_nrows = m_nhist/m_ncols + (m_nhist%m_ncols ? 1 : 0 );
  }

  std::string name() const { return m_name; }

  size_t size() const { return m_hist.size(); }

  const HistDetails& operator[](int i) const { return m_hist.at(i); }
  HistDetails&       operator[](int i)       { return m_hist.at(i); }

  const HistDetails& back() const { return m_hist.back(); }
   HistDetails&      back()       { return m_hist.back(); }

  int nrows() const { return m_nrows; }
  int ncols() const { return m_ncols; }

private:

  std::string               m_name;

  int                       m_nhist;

  int                       m_nrows;
  int                       m_ncols;

  std::vector<HistDetails> m_hist;

};



inline std::ostream& operator<<( std::ostream& s, const Panel& p ) { 
  s << "Panel: " << p.name();
  if ( p.size() == 1 ) s << "\t" << p[0];
  else for ( size_t i=0 ; i<p.size() ; i++ ) s << "\n\t" << p[i];
  return s;
}





#endif // COMPUTILS_H
