// emacs: this is -*- c++ -*-
//
//   @file    lumiParser.h        
//
//                   
//  
//   Copyright (C) 2016 M.Sutton (sutt@cern.ch)    
//
//   $Id: lumiparser.h, v0.0   Tue 12 Jul 2016 13:47:00 CEST sutt $


#ifndef  LUMIPARSER_H
#define  LUMIPARSER_H

#include <iostream>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <cmath>

#include "lumiList.h"


class lumiParser : public lumiList {

public:
  
  lumiParser() { }

  lumiParser( const std::string& file ) { read( file ); }

  void read( const std::string& file ) {

    /// open the file
    std::ifstream input( file );

    int run;
    std::vector<int> lbstart;
    std::vector<int> lbend;
    
    bool first = true;
    
    /// go through all the lines ...
    for( std::string line; getline( input, line ); ) { 
      
      //    std::string s;
      //    std::cout << line << std::endl;
      
      std::string inputline = tolower(line);
      
      /// got a run 
      if ( inputline.find("run:")!=std::string::npos || 
	   inputline.find("<run>")!=std::string::npos ) {
	
	/// parse the run number

	/// really need to check for upper or lower 
	/// case in addition
	std::string label = line;
	chop(  label, "Run:" );
	chop(  label, "<Run>" );
	chomp( label, "</Run>" );
	clean( label, " " );
	
	run = std::atoi(label.c_str());
	lbstart.clear();
	lbend.clear();
	
	//      std::cout << "\n\n" << line << "\trun: " << label << std::endl;

	/// parse the lumiblock ranges
	for( std::string fline ; getline( input, fline ) ; ) { 
	  if ( tolower(fline).find("from:")!=std::string::npos ) { 
	    chop(  fline, "From:" );
	    clean( fline, " ");
	    std::string start = fline; 
	    chomp( start, "to");
	    chop(  fline, "to" );
	    lbstart.push_back( std::atoi(start.c_str()) );
	    lbend.push_back( std::atoi(fline.c_str()) );
	  }
	  else if ( tolower(fline).find("lbrange")!=std::string::npos ) {
	    clean( fline, "\"" );
	    clean( fline, " " );
	    chop(  fline, "LBRange");
	    chomp( fline, "/>" );
	    chop(  fline, "Start=" );
	    std::string start = fline; 
	    chomp( start, "End");
	    chop(  fline, "End=");
	    lbstart.push_back( std::atoi(start.c_str()) );
	    lbend.push_back( std::atoi(fline.c_str()) );
	  }
	  else break;
	}
		
	for ( unsigned  il=0 ; il<lbstart.size() ; il++ ) this->addRange( run, lbstart[il], lbend[il] );
	
      }
    }    
  }

  virtual ~lumiParser() { } 

private:

  /// all these string manipulation routines should all be available
  /// from elsewhere
   
  void chomp( std::string& s, const std::string& regex ) { 
    size_t pos = s.find( regex );
    if ( pos!=std::string::npos ) s.erase( pos, s.size() );
  }
    
  void chop( std::string& s, const std::string& regex ) { 
    size_t pos = s.find( regex );
    if ( pos!=std::string::npos ) s.erase( 0, pos+regex.size() );  
  }
  
  
  void clean( std::string& s, const std::string& regex ) { 
    size_t pos = s.find( regex );
    while( pos!=std::string::npos ) { 
      s.erase( pos, regex.size() );
      pos = s.find( regex );
    }
  }
  
  
  std::string toupper( const std::string& s ) { 
    std::string t = s;
    for ( unsigned i=0 ; i<t.size() ; i++ ) if ( t[i]>='a' && t[i]<='z' ) t[i] += 'A'-'a';  
    return t;
  }
  
  
  std::string tolower( const std::string& s ) { 
    std::string t = s;
    for ( unsigned i=0 ; i<t.size() ; i++ ) if ( t[i]>='A' && t[i]<='Z' ) t[i] -= 'A'-'a';  
    return t;
  }
  
  
};



inline std::ostream& operator<<( std::ostream& s, const lumiParser& _l ) { 
  return s;
}


#endif  // LUMIPARSER_H 










