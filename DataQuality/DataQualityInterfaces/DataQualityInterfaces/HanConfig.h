/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef dqiHanConfig_h
#define dqiHanConfig_h

#include <string>
#include <set>
#include <map>
#include <TList.h>
#include <TMap.h>

#include <TObject.h>

#include "DataQualityInterfaces/MiniConfigTreeNode.h"
#include "DataQualityInterfaces/HanConfigAssessor.h"
#include "DataQualityInterfaces/HanConfigCompAlg.h"
#include "DataQualityInterfaces/DatabaseConfig.h"

#ifndef __CINT__
#include <boost/shared_ptr.hpp>
#endif

class TDirectory;
class TFile;
class TKey;

namespace dqm_core {
  class Input;
  class Output;
  class Region;
}


namespace dqi {

class MiniConfig;
class HanConfigAssessor;
class HanConfigGroup;
class HanOutput;
class HanInputRootFile;

class HanConfig : public TObject {
public:

  HanConfig();
  virtual ~HanConfig();


  virtual void AssembleAndSave( std::string infileName, std::string outfileName,
                                std::string connectionString="sqlite://;schema=/afs/cern.ch/user/a/atlasdqm/dqmdisk1/cherrypy-devel/RefDB.db;dbname=REFDB",
                                long runNumber=2147483646, bool bulk=false);

  virtual void BuildMonitors( std::string configName, HanInputRootFile& input, HanOutput& output );
#ifndef __CINT__
  virtual boost::shared_ptr<dqm_core::Region> BuildMonitorsNewRoot( std::string configName, HanInputRootFile& input, dqm_core::Output& output );
#endif
  virtual void BuildConfigOutput( std::string configName, TFile* inputFile, std::string path,
                                  std::map<std::string,TSeqCollection*>* outputMap, TSeqCollection *outputList );

  virtual TObject* GetReference( std::string& groupName, std::string& name );
  virtual const HanConfigAssessor* GetAssessor( std::string& groupName, std::string& name ) const;

  virtual void GetRegexList( std::set<std::string>& regexlist );

protected:

  typedef std::map<std::string,TDirectory*> DirMap_t;


  class RefVisitor : public MiniConfigTreeNode::Visitor {
  public:
    RefVisitor( TFile* outfile_, HanConfig::DirMap_t& directories_, TMap* refsourcedata );
    virtual void Visit( const MiniConfigTreeNode* node ) const;
  protected:
    TFile* m_outfile;
    HanConfig::DirMap_t& m_directories;
    TMap* m_refsourcedata;
  };

  class RefWriter : public MiniConfigTreeNode::Writer {
  public:
    RefWriter( DatabaseConfig& databaseConfig_, const bool bulk);
    virtual void Write( MiniConfigTreeNode* node);
  protected:
    DatabaseConfig& m_databaseConfig;
    const bool m_bulk;
  };

  class AssessmentVisitorBase : public MiniConfigTreeNode::Visitor {
  public:
    AssessmentVisitorBase( HanConfigGroup* root_, const MiniConfig& algConfig_,
                           const MiniConfig& thrConfig_, const MiniConfig& refConfig_,
                           TFile* outfile_, HanConfig::DirMap_t& directories_,
			   TMap* refsourcedata_ );

  protected:

    void GetAlgorithmConfiguration( HanConfigAssessor* dqpar, const std::string& algID,
                                    const std::string& assessorName = "" ) const;

    HanConfigGroup* m_root;
    const MiniConfig& m_algConfig;
    const MiniConfig& m_thrConfig;
    const MiniConfig& m_refConfig;
    TFile* m_outfile;
    HanConfig::DirMap_t& m_directories;
    TMap* m_refsourcedata;
    // File cache
    mutable std::map<std::string, std::shared_ptr<TFile> > m_filecache;
    // following is so we can skip repeated attempts to open nonexistent files
    mutable std::unordered_set<std::string> m_badPaths;
    // following is a cache of the set of keys in each TFile
    // vector since we are going to iterate through them
    mutable std::map<std::string, std::vector<std::string>> m_keycache;
    std::shared_ptr<TFile> GetROOTFile(std::string& fname) const;
    void PopulateKeyCache(std::string& fname, std::shared_ptr<TFile> file) const;
    void EnsureKeyCache(std::string& fname) const;
  };


  class RegionVisitor : public AssessmentVisitorBase {
  public:
    RegionVisitor( HanConfigGroup* root_, const MiniConfig& algConfig_,
                   const MiniConfig& thrConfig_, const MiniConfig& refConfig_,
                   HanConfig::DirMap_t& directories_ );
    virtual void Visit( const MiniConfigTreeNode* node ) const;
  };


  class AssessmentVisitor : public AssessmentVisitorBase {
  public:
    AssessmentVisitor( HanConfigGroup* root_, const MiniConfig& algConfig_,
                       const MiniConfig& thrConfig_, const MiniConfig& refConfig_,
                       TFile* outfile_, HanConfig::DirMap_t& directories_,
		       TMap* refsourcedata_ );
    virtual void Visit( const MiniConfigTreeNode* node ) const;
  };

#ifndef __CINT__
  class RegexVisitor : public HanConfigAssessor::Visitor {
  public:
    RegexVisitor( std::set<std::string>& regexes_ );
    virtual boost::shared_ptr<dqm_core::Node>
      Visit( const HanConfigAssessor* node, boost::shared_ptr<dqm_core::Region> ) const;
  protected:
    std::set<std::string>& m_regexes;
  };

  class ConfigVisitor : public HanConfigAssessor::Visitor {
  public:
    ConfigVisitor( TFile* file_, dqm_core::Output* output_ );
    virtual boost::shared_ptr<dqm_core::Node>
      Visit( const HanConfigAssessor* node, boost::shared_ptr<dqm_core::Region> dqParent ) const;
  protected:
    TFile* m_file;
    dqm_core::Output* m_output;
  };
#endif

  class CompAlgVisitor : public MiniConfigTreeNode::Visitor {
  public:
    CompAlgVisitor( TFile* outfile_ , const MiniConfig& compAlgConfig_);
    virtual void Visit( const MiniConfigTreeNode* node ) const;
  protected:
    TFile* m_outfile;
    const MiniConfig& m_compAlgConfig;
  };

  class MetadataVisitor : public MiniConfigTreeNode::Visitor {
  public:
    MetadataVisitor( TFile* outfile_ , const MiniConfig& metadataConfig_);
    virtual void Visit( const MiniConfigTreeNode* node ) const;
  protected:
    TFile* m_outfile;
    const MiniConfig& m_metadataConfig;
  };


  bool Initialize( const std::string& configName );


  TFile*             m_config;
#ifndef __CINT__
  boost::shared_ptr<dqm_core::Region>  m_dqRoot;
#endif
  HanConfigGroup*    m_top_level;
  TSeqCollection*    m_metadata;

//Get rid of Root macros that confuse Doxygen
///\cond CLASSDEF
  ClassDef( HanConfig, 0 ) // Creates a Han configuration from a MiniConfig
///\endcond

private:

  static TKey* GetObjKey( TDirectory* dir, std::string path );
  static TDirectory* ChangeInputDir( TDirectory* dir, const std::string& path );
  static TDirectory* ChangeOutputDir( TFile* file, const std::string& path, DirMap_t& directories );

};

} // namespace dqi

#endif
