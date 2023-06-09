/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//$Id: FCdeletePFN.cpp 509054 2012-07-05 13:33:16Z mnowak $
/** FCdeletePFN.cpp -- FileCatalog command line tool to delete the selected PFN 
    @author: Zhen Xie
    @date: 02/03/2005 Z.X.
    set default logging to Warning if no POOL_OUTMSG_LEVEL is set; 
    separate logging stream to std::cerr, output stream to std::cout.
    @date: 07/04/2005 Z.X.
    adopt to split interface
*/
#include "FileCatalog/CommandLine.h"
#include "FileCatalog/IFileCatalog.h"
#include "FileCatalog/FCException.h"
#include "POOLCore/Exception.h"
#include "FileCatalog/IFCAction.h"
#include "FileCatalog/FCEntry.h"
#include "FileCatalog/IFCContainer.h"
#include "CoralBase/MessageStream.h"
#include "CoralBase/MessageStream.h"
#include <memory>
using namespace pool;

void printUsage(){
  std::cout<<"usage: FCdeletePFN [-q query -p pfname -u contactstring -h]" <<std::endl;
}

static const char* opts[] = {"p","q","u","h",0};


int main(int argc, char** argv)
{
  std::string  myuri;
  std::string  mypfn;
  std::string  myquery;
  try{
    CommandLine commands(argc,argv);
    commands.CheckOptions(opts);

    if( commands.Exists("u") ){
      myuri=commands.GetByName("u");
    }    
    if( commands.Exists("q") ){
      myquery=commands.GetByName("q");
    }
    if( commands.Exists("p") ){
      mypfn=commands.GetByName("p");
    }
    if( commands.Exists("h") ){
      printUsage();
      exit(0);
    }
  }catch(std::string& strError){
    std::cerr<< "Error: command parsing error "<<strError<<std::endl;
    exit(0);
  }
  
  if( mypfn.empty() && myquery.empty() ){
    printUsage();
    std::cerr<<"Error: must specify pfname using -p, query using -q"<<std::endl;
    exit(0);
  }
  try{
    std::auto_ptr<IFileCatalog> mycatalog(new IFileCatalog);
    mycatalog->setWriteCatalog(myuri);
    FCAdmin a;
    mycatalog->setAction(a);
    FClookup l;
    mycatalog->setAction(l);
    mycatalog->connect();
    mycatalog->start();
    if( !myquery.empty() ){
      PFNContainer pfns(mycatalog.get());
      FileCatalog::FileID fid;
      l.lookupPFNByQuery(myquery,pfns);
      while(pfns.hasNext()){
        PFNEntry pentry=pfns.Next();
        a.deletePFN(pentry.pfname());
      }
    }else if( !mypfn.empty() ) {
      a.deletePFN(mypfn);
    }
    mycatalog->commit();  
    mycatalog->disconnect();
  }catch (const pool::Exception& er){
    //er.printOut(std::cerr);
    //std::cerr << std::endl;
    std::cerr<<er.what()<<std::endl;
    exit(1);
  }catch (const std::exception& er){
    std::cerr<<er.what()<<std::endl;
    exit(1);
  }
}


