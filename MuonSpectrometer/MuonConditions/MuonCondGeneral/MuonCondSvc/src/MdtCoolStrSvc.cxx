/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "AthenaPoolUtilities/AthenaAttributeList.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeListSpecification.h"
// temporary includes to access CLOBs
//#include "CoolKernel/ExtendedAttributeListSpecification.h"
//#include "CoolKernel/PredefinedStorageHints.h"

#include "MuonCondSvc/MdtCoolStrSvc.h"
#include "MuonCondSvc/MdtStringUtils.h"
// SEAL
//#include "SealBase/Time.h"
// STL Streams
#include <fstream>
#include <sstream>

// technology types for storage
#define INLINE_STRING 0
#define INLINE_CLOB   1
namespace MuonCalib {

MdtCoolStrSvc::MdtCoolStrSvc(const std::string& name, ISvcLocator* svc) :
  AthService(name,svc),
  p_detstore(0)
{
  // declare properties
}

MdtCoolStrSvc::~MdtCoolStrSvc() {}

const InterfaceID& MdtCoolStrSvc::type() const
{
  return MdtICoolStrSvc::interfaceID();
}

StatusCode MdtCoolStrSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if (MdtICoolStrSvc::interfaceID().versionMatch(riid)) {
    *ppvInterface=(MdtICoolStrSvc*)this;
  } else {
    return AthService::queryInterface(riid,ppvInterface);
  }
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::initialize()
{
  // service initialisation 

  MsgStream log(msgSvc(),name());

  log << MSG::DEBUG << "in initialize()" << endreq;

  // get detector store

  
  if (StatusCode::SUCCESS!=service("DetectorStore",p_detstore)) {
    log << MSG::FATAL << "Detector store not found" << endreq; 
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::finalize()
{
  MsgStream log(msgSvc(),name());
  log << MSG::DEBUG << "in finalize()" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::putFileT0(const std::string& folder, 
   const std::string& filename, const int chan, const int tech) const {
  MsgStream log(msgSvc(),name());
  log << MSG::INFO << "PutFile for file " << filename << " folder " <<
    folder << " chan " << chan << " technology " << tech << endreq;
  
  std::ifstream f(filename.c_str());
  std::string sdata="";
  std::string type;
  std::string deli_data = ",";
  std::string delimiter = ".";
  std::vector<std::string> tokens;
  MuonCalib::MdtStringUtils::tokenize(filename,tokens,delimiter);
  sdata+=tokens[0]+deli_data;

  if (f != NULL) {
    std::string blobline;  
    std::string blob_header = "";
    std::string version, region, ntubes;
  
    while(getline(f,blobline)) {
     
      std::string deli_data = ",";
      std::string delimiter = " ";
      std::vector<std::string> tokens;
      MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
      type = tokens[0];
      
      if (type.find("v0.0")==0) {
	std::string delimiter = " ";
	std::vector<std::string> tokens;
	MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
	blob_header = blobline;
	version= tokens[0];
	region= tokens[1];
	ntubes= tokens[5];
	sdata+=version + deli_data + region + deli_data +ntubes  + deli_data + " ";
      } else{
	std::string delimiter = " ";
	std::vector<std::string> tokens;
	MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
	//	for (unsigned int i=0; i<tokens.size(); i++) {
	  std::string t0= tokens[7];
	  std::string  adc= tokens[8];
	  std::string  stat= tokens[9];
	  sdata+=t0 +  deli_data + stat + deli_data + adc + deli_data;
	  //}
      }
      
    }
    sdata+="  end";
    f.close();
    
    std::string sdata_t0;
    sdata_t0+= sdata;
    int size_fin = sdata_t0.size();
    std::cout << "size of fin " << size_fin << std::endl;
    
    putData(folder,filename,chan,tech,sdata_t0 ); 
  } else {
    log << MSG::INFO << "Cannot open file " << filename << endreq;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}	

StatusCode MdtCoolStrSvc::putFileRT(const std::string& folder, 
	          const std::string& filename, const int chan, const int tech) const {
  MsgStream log(msgSvc(),name());
  log << MSG::INFO << "PutFile for RT file " << filename << " folder " <<
    folder << " chan " << chan << " technology " << tech << endreq;
  // transform to a string
  FILE* f = fopen (filename.c_str(),"rb");
  int size;
  char string[255];
  if (f) {
      fseek (f, 0L, SEEK_END);
      size = ftell (f);
      fseek (f, 0L, SEEK_SET);
      
      fgets (string , 255 ,f );
  }

  //puts (string);
  char header[16000]={""};
 
  char * pch;
  int i=0;
  int toSkip=0;
 
  pch = strtok (string," ");
  // check data format
  if (!strcmp(pch,"v0.0")) {
    std::cout << "data format 0 " << pch <<std::endl;
    toSkip=3;
  } else if ( !strcmp(pch,"v1.0") ) {
    std::cout << "data format 1 " << pch <<std::endl;
    toSkip=10;
  } else {
    std::cout << "Unknown data format " << pch <<std::endl;
  }

  while (pch != NULL)
  {
    if (i==1) {
      int numberOfRts=atoi(pch);
      if(numberOfRts!=1) { 
         std::cout<<"ERROR ERROR ERROR ERROR ERROR ERROR ERROR"<<std::endl; 
         std::cout<<" more than one Rt "<<numberOfRts
           <<" in the same file, format currently not supported!!!"<<std::endl;
         std::cout<<"ERROR ERROR ERROR ERROR ERROR ERROR ERROR"<<std::endl; 
	  fclose (f);
         return StatusCode::FAILURE;
      }
    }
    if (i==toSkip) {
      // region number
      if (sizeof(header)>=sizeof(pch)) strcpy(header, pch);
      
      if (sizeof(header)>=sizeof(pch)+1) strcat(header, ",");
    }
    if (i==toSkip+1) {
      // number of rt points
      if (sizeof(header)>=sizeof(pch))  strcat(header, pch);
      if (sizeof(header)>=sizeof(pch)+1) strcat(header, " ");
    }
    pch = strtok (NULL, " ");
    i++;
    // in the header version + # region
  }
  log << "header finale "<< endreq;
  puts(header);

  if (f != NULL)   {
    log << MSG::INFO << "Input file size is " << size << endreq;
    char pack[1000];
    while(!feof(f))
      {
	float rad; float sigma; float time;
	
	int ret = fscanf(f,"%f %f %f", &rad, &time, &sigma);
	//printf("\n %8.3f %8.3f %8.3f \n",rad,time,sigma);
	if (ret!=0){
	  char * xmlt0;
	  asprintf (&xmlt0, "%f,%f,%f,", rad, time, sigma);
	  if (sizeof(pack)>=sizeof(xmlt0)) strcpy(pack,xmlt0);
	  if (sizeof(header)>=sizeof(pack)) strcat(header,pack);
	}
      }

    std::string sdata = strcat(header," end ");
    int size_rt = sdata.size();
    std::cout << size_rt << std::endl;
    fclose (f);
    std::string sdata_rt;

      sdata_rt+= sdata;
      int size_fin2 = sdata_rt.size();
      std::cout << "size of fin " << size_fin2 << std::endl;

    putData(folder,filename,chan,tech,sdata_rt );
  } else {
    log << MSG::INFO << "Cannot open file " << filename << endreq;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}// fine RT

// Start Align
StatusCode MdtCoolStrSvc::putFileAlignCorr(const std::string& folder, 
   const std::string& filename, const int chan, const int tech) const {
  MsgStream log(msgSvc(),name());
  log << MSG::INFO << "PutFile for file " << filename << " folder " <<
    folder << " chan " << chan << " technology " << tech << endreq;

   // Reads the ascii file and write it in a string
  std::ifstream f(filename.c_str());
  std::string sdata="";

  if (f != NULL) {

    std::string blobline;
    
    // Define the Header
    std::string blob_header = "";
    // Check the first word to see if it is a correction
    std::string type;
    //Get Timestamp from Header: define variables to store it
    std::string since_str;
    std::string till_str;

    while(getline(f,blobline)) {

      sdata += blobline;
      sdata += "\n";
      std::string delimiter = ":";
      std::vector<std::string> tokens;
      MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
      type = tokens[0];
      //Parse line
      if (type.find("#")==0) {
	//skip it
	continue;
      }
      //      std::cout << type ;
      if (type.find("Header")==0) {
	std::string delimiter = "|";
	std::vector<std::string> tokens;
	MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
	blob_header = blobline;
	since_str = tokens[1];
	till_str = tokens[2];
      }

      if (type.find("Corr")==0) {
	std::string delimiter = " ";
	std::vector<std::string> tokens;
	MuonCalib::MdtStringUtils::tokenize(blobline,tokens,delimiter);
	//	for (unsigned int i=0; i<tokens.size(); i++)
	//	  std::cout << " " << tokens[i];
	//	std::cout << std::endl;
      }
    
      
      std::cout << "End of line" << std::endl;
      
    } 

   //Retrieve TimeStamp information stored in the header: index 0 in tokens
   const char* since = since_str.c_str();
   const char* till = till_str.c_str();
   int year, month, day, hour, minute, second, ns;
   sscanf(since,"%d/%d/%d %d:%d:%d.%d",&year,&month,&day,&hour,&minute,&second,&ns);
   log << MSG::INFO << "since : " << year << "/" << month << "/" << day << " " 
	<< hour << ":" << minute << ":" << second << endreq;
   //  seal::Time since_T(year,month,day,hour,minute,second,0,true);

   sscanf(till,"%d/%d/%d %d:%d:%d.%d",&year,&month,&day,&hour,&minute,&second,&ns);
   //seal::Time till_T(year,month,day,hour,minute,second,0,true);

   //log << MSG::INFO << "SINCE Seal " << since_T.format(true,"%F %T") << endreq;
   //log << MSG::INFO << "TILL Seal " << till_T.format(true,"%F %T") << endreq;

   log << MSG::INFO << "End of blob parsing" << endreq;
   log << MSG::INFO << "Dumping BLOB : " << sdata << endreq;


   f.close();
   //return StatusCode::SUCCESS;
   putData(folder,filename,chan,tech,sdata );
  } else {
    log << MSG::INFO << "Cannot open file " << filename << endreq;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}// End Align



StatusCode MdtCoolStrSvc::putAligFromFile(const std::string& folder, 
    const std::string& filename, const int chan, const int tech) const {
  MsgStream log(msgSvc(),name());
  log << MSG::INFO << "PutFile for file " << filename << " folder " <<
    folder << " chan " << chan << " technology " << tech << endreq;
  FILE* f = fopen (filename.c_str(),"rb");
  fseek (f, 0L, SEEK_END);
  int size = ftell (f);
  fseek (f, 0L, SEEK_SET);
  
  std::string sdata;
  if (size!=0)   {
    log << MSG::INFO << "Input file size is " << size << endreq;
    char pack[1000];
    while(!feof(f)){
      
      char line[300];
      fgets (line , 300 ,f ); // header of the file
      //puts(line);
      
      int i=0;
      char * pch;
      //printf ("Splitting string \"%s\" in tokens:\n",line);
      pch = strtok (line," \n\t");
      std::string clob;  
      while (pch != NULL)
	{
	  if (i<11) {
	    //printf ("%s\n",pch);
	    //name=pch;
	    if (sizeof(pack)>=sizeof(pch)) strcpy(pack,pch);
	    if (sizeof(pack)>=sizeof(pch)+1)strcat(pack, ",");
	    if (pack != NULL) {
	      sdata += pack;
	    }
	  }
	  i++;
	  pch = strtok (NULL, " \n\t");
	}
    }
    
    //    int size_t = sdata.size();
    sdata += " end";
    std::cout<< sdata<< std::endl;
    fclose (f);
    
    putData(folder,filename,chan,tech,sdata );
  } else {
    log << MSG::INFO << "Cannot open file " << filename << endreq;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
  
}

//dead tube
StatusCode MdtCoolStrSvc::putFileTube(const std::string& folder, 
   const std::string& filename, const int chan, const int tech) const {
  MsgStream log(msgSvc(),name());
  log << MSG::INFO << "PutFile for file " << filename << " folder " <<
    folder << " chan " << chan << " technology " << tech << endreq;
   // transform to a string
  
  FILE* f = fopen (filename.c_str(),"rb");
  fseek (f, 0L, SEEK_END);
  int size = ftell (f);
  fseek (f, 0L, SEEK_SET);
  char string [100];
  fgets (string , 100 ,f ); // header of the file
  puts(string);
  //char * fname=filename.c_str();
  
  std::string sdata;
  // log << MSG::INFO << "entro nel loop"  << endreq;

  if (size!=0)   {
    log << MSG::INFO << "Input file size is " << size << endreq;
    
    while(!feof(f)){
      
      char line[200];
      fgets (line , 200 ,f ); // header of the file
      puts(line);
     
      int i=0;
      char * pch;
      printf ("Splitting string \"%s\" in tokens:\n",line);
      pch = strtok (line,",");
      char pack[1000];
      while (pch != NULL)
	{
	  if (i==1) {
	    printf ("%s\n",pch);
	    //name=pch;
	    if (sizeof(pack)>=sizeof(pch)) strcpy(pack,pch);
	    if (sizeof(pack)>=sizeof(pch)+1) strcat(pack, ",");
	  }
	   if (i==4) {
	    printf ("%s\n",pch);
	    //mlayer=pch;
	    if (sizeof(pack)>=sizeof(pch)) strcat(pack,pch);
	    if (sizeof(pack)>=sizeof(pch)+1) strcat(pack, ",");
	  }
	   if (i==5) {
	    printf ("%s\n",pch);
	    //layer=pch;
	    if (sizeof(pack)>=sizeof(pch)) strcat(pack,pch);
	    if (sizeof(pack)>=sizeof(pch)+1) strcat(pack, ",");
	  }
	    if (i==6) {
	    printf ("%s\n",pch);
	    //tube=pch;
	    if (sizeof(pack)>=sizeof(pch)) strcat(pack,pch);
	    if (sizeof(pack)>=sizeof(pch)+1) strcat(pack, ",");
	  }
	    i++;
	  pch = strtok (NULL, " ,");
	}
      
      if (sizeof(pack)!=0) sdata += pack;
      std::cout<< sdata<< std::endl;

    }
      
    sdata += " end";
    fclose (f);
    putData(folder,filename,chan,tech,sdata );
  } else {
    log << MSG::INFO << "Cannot open file " << filename << endreq;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::getString(const std::string& folder, 
		      const int chan, std::string& data) const {
  // dummy datafile name - this is not returned to the client
  std::string rfile;
  StatusCode sc=getData(folder,chan,rfile,data);
  return sc;
}

StatusCode MdtCoolStrSvc::getFile(const std::string& folder, const int chan,
				   const std::string& file) const {
  MsgStream log(msgSvc(),name());
  std::string rfile,data;
  if (StatusCode::SUCCESS!=getData(folder,chan,rfile,data)) {
    return StatusCode::FAILURE;
  }
  // default to filename stored with data
  if (file!="") rfile=file;

  FILE* f = fopen (rfile.c_str(),"wb");
  if (f!=NULL) {
    int size=data.size();
    fwrite(data.c_str(),size,1,f);
    log << MSG::INFO << "getFile: written data of length " << size <<
    " into file " << rfile << endreq;
  } else {
    log << MSG::ERROR << "Failed to open file " << rfile << " for write" <<
      endreq;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::putData(const std::string& folder, 
        const std::string& filename,const int chan,const int tech, 
				   const std::string& data) const  {
  MsgStream log(msgSvc(),name());

  // check technology type
  if (tech!=INLINE_STRING && tech!=INLINE_CLOB) {
    log << MSG::ERROR << "Bad technology type (must be 0 or 1)" << endreq;
    return StatusCode::FAILURE;
  }
  // check data length for strings
  if (tech==INLINE_STRING && data.size()>=4000) {
    log << MSG::ERROR << "String technology selected (" << INLINE_STRING << 
      ") but data size is >4000 characters " << endreq;
    return StatusCode::FAILURE;
  }
  // check channel zero is not being requested
  if (chan==0) {
    log << MSG::ERROR << "Channel 0 cannot be used" << endreq;
    return StatusCode::FAILURE;
  }

  // check if collection already exists, if not have to create in TDS
  log << MSG::INFO << "putData to folder " << folder << " channel " << chan <<
    " from file " << filename << " using technology ";
  if (tech==INLINE_STRING) log << INLINE_STRING << " (String)" << endreq;
  if (tech==INLINE_CLOB) log << INLINE_CLOB << " (CLOB<16M)" << endreq;
  CondAttrListCollection* atrc=0;
  if (!p_detstore->contains<CondAttrListCollection>(folder)) {
    log << MSG::DEBUG << "Creating new CondAttrListCollection for folder "
	<< folder << endreq;
    CondAttrListCollection* atrc=new CondAttrListCollection(true);
    if (StatusCode::SUCCESS!=p_detstore->record(atrc,folder)) {
    log << MSG::ERROR << "Could not create CondAttrListCollection " <<
     folder << endreq;
    return StatusCode::FAILURE;
    }
  }
  // do const cast here so we can add to already exisiting collections
  const CondAttrListCollection* catrc=0;
  log << MSG::DEBUG << "Attempting to retrieve collection (const)" << endreq;
  if (StatusCode::SUCCESS!=p_detstore->retrieve(catrc,folder)) {
    log << MSG::ERROR << "Could not retrieve CondAttrListCollection " <<
       folder << endreq;
    return StatusCode::FAILURE;
  }
  atrc=const_cast<CondAttrListCollection*>(catrc);
  if (atrc==0) {
    log << MSG::ERROR << "Could not retrieve non-const pointer to atrc" <<
      endreq;
    return StatusCode::FAILURE;
  }

  log << MSG::DEBUG << "About to create AttributeListSpecification" << endreq;
  coral::AttributeListSpecification* aspec=0;
  // different things for string and CLOB technlogies
  aspec=new coral::AttributeListSpecification();
  aspec->extend("tech","int");
  aspec->extend("file","string");
  aspec->extend("data","string");
  //  if (tech==INLINE_STRING) {
  //aspec=new coral::AttributeListSpecification();
  //aspec->extend("tech","int");
  //aspec->extend("file","string");
  //aspec->extend("data","string");
  //} else if (tech==INLINE_CLOB) {
  //boost::shared_ptr<cool::ExtendedAttributeListSpecification> 
  //  espec( new cool::ExtendedAttributeListSpecification());
  //espec->push_back("tech","int");
  //espec->push_back("file","string");
  //espec->push_back("data","string",
  //	     "PredefinedStorageHints::STRING_MAXSIZE_16M");
  //aspec=new coral::AttributeListSpecification(
  //		     espec->attributeListSpecification());
  //}
  // define data contents
  AthenaAttributeList alist(*aspec);
  alist["tech"].setValue(tech);
  alist["file"].setValue(filename);
  //data=*(static_cast<const std::string*>((atr["data"]).addressOfData()));
  alist["data"].setValue(data);
  CondAttrListCollection::ChanNum channum=chan;

  std::ostringstream atrstring;
  alist.print(atrstring);
  log << MSG::DEBUG << "About to add channel to: " << atrc << endreq;
  atrc->add(channum,alist);
  return StatusCode::SUCCESS;
}

StatusCode MdtCoolStrSvc::getData(const std::string& folder, 
	    const int chan, std::string& file, std::string& data) const {
  MsgStream log(msgSvc(),name());
  const CondAttrListCollection* atrc;
  data="";
  if (StatusCode::SUCCESS!=p_detstore->retrieve(atrc,folder)) {
    log << MSG::ERROR << "getData failed for folder " << folder << " channel "
	<< chan << endreq;
    return StatusCode::FAILURE;
  }
  CondAttrListCollection::ChanNum channum=chan;
  CondAttrListCollection::const_iterator itr=atrc->chanAttrListPair(channum);
  if (itr!=atrc->end()) {
    const coral::AttributeList& atr=itr->second;
    //std::string data;
    //data=*(static_cast<const std::string*>((atr["data"]).addressOfData()));
    //atr["data"].getValue(data);
    //atr["file"].getValue(file);
    //std::ostringstream atrstring;
    //atr.print(atrstring);
    data=atr["data"].data<std::string>();
    file=atr["file"].data<std::string>();
    if (log.level() < MSG::INFO) {
      std::ostringstream atrstring;
      atr.toOutputStream(atrstring);
      log << MSG::DEBUG << "read Channum " << channum << " atrlist: " << 
	atrstring.str() << endreq;
    } else {
      log << MSG::ERROR << "Invalid channel number" << endreq;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

}//namespace
