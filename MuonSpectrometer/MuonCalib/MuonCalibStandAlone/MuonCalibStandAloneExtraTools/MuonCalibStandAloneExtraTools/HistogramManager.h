/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 13.08.2008, AUTHOR: MAURO IODICE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#ifndef SRC_HISTOGRAMMANAGER_H
#define SRC_HISTOGRAMMANAGER_H

#include "TMath.h"
#include <string.h>
#include <stdio.h>

#include <sstream>
#include <iostream>
#include <fstream> 

#include <TFile.h>
#include <TString.h>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"

#include "MuonCalibStandAloneExtraTools/StringUtil.h"
#include "MuonCalibStandAloneExtraTools/MDTName.h"

#define UPDATETIME 2000

using namespace std;
using namespace MuonCalib;

class MdtIdHelper ;

class MdtChamber  {
public:
  std::string m_region;
  std::string m_side;
  int m_sector;
  int m_absEta;
  int m_eta_id;
  int m_phi_id;
  std::string m_chamberType;
  std::string m_chamberName;
  std::string m_athenaChamberName;

  void Fill(std::string region, std::string side, int sector, std::string chamberType, int absEta_id) {
       ToString _ts ;
       string sectorString=_ts(sector);
       if(sector<10) sectorString="0"+sectorString;
       m_region = region ;
       m_side = side ;
       m_sector = sector ;
       m_chamberType = chamberType ;
       m_absEta = absEta_id ;
       int eta_id = absEta_id ;
       if (side=="C") eta_id = -absEta_id ;
       m_eta_id = eta_id ;
       int phi_id = (int) (sector+1)/2 ;
       m_phi_id = phi_id ;
       string chName = chamberType+_ts(absEta_id)+side+sectorString ;
       m_chamberName = chName ;
       string chAthenaName = chamberType+"_"+_ts(phi_id)+"_"+_ts(eta_id) ;
       m_athenaChamberName = chAthenaName ;
       return ;
  }
};

struct sortMdtChambersByName
{
     bool operator()(MDTName a ,MDTName b)
     {
       string t1=a.getOnlineName();
       string t2=b.getOnlineName();
       return (t1 < t2);
     }
};

class HistogramManager {
public:
  HistogramManager();
  HistogramManager(const MdtIdHelper * mdtIdHelper);
  ~HistogramManager();
  void buildGlobalHistos();
  void buildTrackHistos(); 
  void buildDebugHistos();
  void buildTopLevel(string region, string side, int sectorMin, int sectorMax);
  void buildSector(string region,string side, int sector);
  void buildChamberHistos(MDTName) ;

  void setChamberCutOut(string chamber, TH1F * href ) ;
  void setChamberDisconnectedTubes(string chamber, TH1F * href ) ;
  int  GetTubeOffsetML1(string chamberName) ;
  int  GetTubeOffsetAtEndML1(string chamberName) ;

  std::vector<MDTName> GetChamberList(string region, string side, int sector) ;

  bool openOutputFile(string filename);
  bool openReadOnlyFile(string filename);
  bool openUpdateFile(string filename);
  void WriteAndCloseFile();  

  void SetDoTracks(bool);

  TObject * GetMdtHisto(string histo_name);
  TObject * GetMdtHisto(string histo_name, string region, string side);
  TObject * GetMdtHisto(string histo_name, string region, string side, int sector);
  TObject * GetMdtHisto(string histo_name, MDTName);

  string GetMdtDirectoryName(); 
  string GetMdtDirectoryName(string region, string side);
  string GetMdtDirectoryName(string region, string side, int sector);
  string GetMdtDirectoryName(string region, string side, int sector, string chamberType, int eta);
  string GetMdtDirectoryName(MDTName); 

  TObject * GetTDaqHisto(string histo_name, string region) ;
  TObject * GetTDaqHisto(string histo_name, string region, string side);
  TObject * GetTDaqHisto(string histo_name, string region, string side, int sector);
  TObject * GetTDaqHisto(string histo_name, string region, string side, int sector, string chamberType, int eta);

  string GetTDaqDirectoryName(string region); 
  string GetTDaqDirectoryName(string region, string side);
  string GetTDaqDirectoryName(string region, string side, int sector);
  string GetTDaqDirectoryName(string region, string side, int sector, string chamberType, int eta);

  TObject * GetHisto(string main_dir, string histo_name);

  void ReadChamberMapFile (string chamberName, int * chamberGeoParams, int numParams);
  int GetChamberNumOfML (string chamberName);
  int GetChamberNumOfMezzPerML (string chamberName);
  int GetChamberTubesPerMezz (string chamberName);

  // std::vector<string> GetChamberTypeList(string region, string side, int sector) ;
  int GetEtaMax(string region, string side, int sector, string chamberType) ;

  // void buildAll(int sectorMin, int sectorMax);
  // void test();
 
  TFile * rootFile(){return m_rootfile;} ;
  TFile * m_rootfile ;

  private:

  TObjArray m_hList;
  const MdtIdHelper * m_MdtIdHelper ;

  bool m_doTracks;

  // ToString _ts;
  // ToChar _tc;
  
};

////////////////////////////////////////////////////////////////////
//
// UTILITIES of General usage
//

/*********************
//
class ToString{
  public:
  template< class T >
    std::string operator()( const T& i )
  {
    std::ostringstream os;
    os << i;
    return os.str();
  }
};
*********************************/


/*********************************
class ToChar {
 public:
    char* operator()( const string i )
  {
    return i.c_str();
  }
    char* operator()( const int i )
  {
    ToString ts;
    return ts(i).c_str();
  }
    char* operator()( const float i )
  {
    ToString ts;
    return ts(i).c_str();
  }
    char* operator()( const char * i )
  {
    return i;
  }

};


*********************************/

#endif //SRC_HISTOMANAGER_H

