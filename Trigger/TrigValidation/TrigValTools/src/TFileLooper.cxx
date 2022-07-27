/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file   TFileLooper.cxx
 * @brief  TFileLooper implementation
 * @author Frank Winklmeier
 *
 */

#include "TrigValTools/TFileLooper.h"

#include <iostream>
#include <sstream>
#include "TError.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TKey.h"
#include "TPRegexp.h"

using namespace std;

Int_t TFileLooper::run(const char* filename, const char* rootDir)
{
  m_errorCode = 0;

  if (filename==0) {
    cout << "Invalid file name (0)" << endl;
    m_errorCode = 1;
    return m_errorCode;
  }
    
  beginJob();  
  processFile(filename,rootDir);
  endJob();
  return m_errorCode;
}


// Convert one file into html
void TFileLooper::processFile(const char* filename, const char* rootDir)
{    
  m_file = new TFile(filename);
  if (m_file->IsZombie()) {
    cout << "Cannot open "<<filename << endl;
    m_errorCode = 1;
    delete m_file;    
    return;
  }

  if (rootDir) {
    if (!m_file->cd(rootDir)) {
      cout << "Cannot change to directory " << rootDir << endl;
      m_errorCode = 1;
      return;
    }
    m_rootDir = rootDir;
  }
  else m_file->cd();
  
  beforeFile();
  processDir(*gDirectory);
  afterFile();
  
  delete m_file;   // calls Close()
}


// Save all histograms in dir into png files
void TFileLooper::processDir(TDirectory& dir)
{
  TString s(getPathFromDir(dir));

  if (skipDir(dir)) {
    cout << "Skipping " << s << endl;
    return;
  }

  if (m_verbose) cout << "Reading directory "<< dir.GetPath() << endl;
        
  // Sort directory content
  TList* dirList = dir.GetListOfKeys();
  dirList->Sort();

  TIter next(dirList);
  TKey* key;
  while ((key = (TKey*)next())) {
    //    if (skipObject(getKeyPath(dir,*key))) continue;

    TString className(key->GetClassName());
    if (className=="TDirectoryFile" || className=="TDirectory") {
      dir.cd(key->GetName());
      beforeDir();      
      processDir(*gDirectory);
      afterDir();
    }
    else {
      if (skipObject(getKeyPath(dir,*key))) {
        m_skippedObjects.push_back(getKeyPath(dir,*key).Data());
        continue;
      }
      processKey(dir, *key);
    }
  }
}

void TFileLooper::processKey(TDirectory& dir, TKey& key)
{
  cout << "--- processKey: " << key.GetName()
       << " in directory " << dir.GetPath() << endl;
}


// Check if we have to skip this directory
Bool_t TFileLooper::skipDir(const TDirectory& dir)
{
  if (m_skipDirs.Contains(dir.GetName())) return kTRUE;
  else return kFALSE;
}


// Check if we have to skip this object
// Logic: 1) check if name matches any of the m_failRE
//        2) check if name matches any of the m_passRE
//        3) Return true if 1)=true and 2)=false
Bool_t TFileLooper::skipObject(const char* name)
{
  Bool_t failMatch(kFALSE);
  Bool_t passMatch(kFALSE);
  
  for (auto& re : m_failRE) {
    if (re.Match(name)>0) {
      failMatch = kTRUE;
      break;
    }
  }

  // give object another chance to match any of the m_passRE
  for (auto& re : m_passRE) {
    if (re.Match(name)>0) {
      passMatch = kTRUE;
      break;
    }
  }

  bool result;
  if (m_passBeforeFail) result = (!passMatch || failMatch);
  else result = (failMatch && !passMatch);

  if (m_verbose && result) cout << "Skipping " << name << endl;
  return result;
    
}


// Skip keys that match this regexp...
void TFileLooper::addFailRegexp(const char* regexp)
{
  if (regexp) {
    m_failRE.emplace_back(regexp);
  }
}


// ...and do not match this regexp.
void TFileLooper::addPassRegexp(const char* regexp)
{
  if (regexp) {
    m_passRE.emplace_back(regexp);
  }
}

// Create file system path from TDirectory path
// Example: "../Tau-monitoring.root:/TrigSteering" -> "Tau-monitoring.root/TrigSteering"
TString TFileLooper::getPathFromDir(const TDirectory& dir)
{
  TString s(dir.GetPath());
  if (s.Index("../")==0) s.Remove(0,3);
  s.ReplaceAll(":","");
  return s;
}

// return key name with path
// Example: myHist in directory myDir -> "myDir/myHist"
TString TFileLooper::getKeyPath(const TDirectory& dir, const TKey& key)
{
  TString s(dir.GetPath()+TString("/")+key.GetName());
  // Remove all chars from beginning including ":/"
  s.Replace(0,s.First(":")+2,0);

  return s;
}
