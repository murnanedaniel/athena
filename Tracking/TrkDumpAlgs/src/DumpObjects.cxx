#include "TrkDumpAlgs/DumpObjects.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "PixelReadoutGeometry/PixelModuleDesign.h"
#include "SCT_ReadoutGeometry/SCT_ModuleSideDesign.h"
#include "InDetPrepRawData/PixelClusterContainer.h"
#include "InDetPrepRawData/SCT_ClusterContainer.h"
#include "TrkSpacePoint/SpacePointOverlapCollection.h"
#include "InDetSimData/InDetSimDataCollection.h"
#include "ReadoutGeometryBase/SiLocalPosition.h"
#include "TrkSpacePoint/SpacePointContainer.h"
#include "GeneratorObjects/xAODTruthParticleLink.h"
#include "GeneratorObjects/McEventCollection.h"
#include "xAODTruth/TruthVertex.h"
#include "InDetPrepRawData/SiCluster.h"

#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include "HepPDT/ParticleDataTable.hh"

#include "TrkTrack/TrackCollection.h"
#include "TrkTrack/TrackInfo.h"
#include "TrkEventPrimitives/FitQuality.h" 
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
#include "InDetRIO_OnTrack/SiClusterOnTrack.h"
#include "InDetRIO_OnTrack/PixelClusterOnTrack.h"
#include "InDetRIO_OnTrack/SCT_ClusterOnTrack.h"

#include "TrkTruthData/TrackTruthCollection.h"
#include "TrkTruthData/DetailedTrackTruthCollection.h"

#include "GaudiKernel/ITHistSvc.h"
#include "TTree.h"

//-------------------------------------------------------------------------
DumpObjects::DumpObjects(const std::string& name, ISvcLocator *pSvcLocator)
//-------------------------------------------------------------------------
  : AthAlgorithm(name, pSvcLocator), 
  m_pixelID(nullptr),
  m_SCT_ID(nullptr),
  m_pixelManager(nullptr),
  m_SCT_Manager(nullptr),
  m_event(0),
  m_selected(0),
  m_outfile(nullptr),
  m_particlePropSvc("PartPropSvc",name),
  m_particleDataTable(0),
  m_offset(0),
  m_dumpTruthParticles(true),
  m_dumpClusters(true),
  m_dumpSpacePoints(true)
  {
    declareProperty("Offset", m_offset);
    declareProperty("FileName", m_name = "");
    declareProperty("DumpTruthParticles", m_dumpTruthParticles = true);
    declareProperty("DumpClusters", m_dumpClusters = true);
    declareProperty("DumpSpacePoints", m_dumpSpacePoints = true);
    //
    declareProperty("NtupleFileName", m_ntupleFileName);
    declareProperty("NtupleDirectoryName", m_ntupleDirName);
    declareProperty("NtupleTreeName", m_ntupleTreeName);
    declareProperty("maxCL", m_maxCL = 1500000);
    declareProperty("csvFile", m_csvFile);
    declareProperty("rootFile", m_rootFile);
  }


//----------------------------------
StatusCode DumpObjects::initialize() {
//----------------------------------
  m_event = m_offset;

  // Grab PixelID helper
  if (detStore()->retrieve(m_pixelID, "PixelID").isFailure()) {
    return StatusCode::FAILURE;
  }

  if (detStore()->retrieve(m_pixelManager, "Pixel").isFailure()) {
    return StatusCode::FAILURE;
  }

  // Grab SCT_ID helper
  if (detStore()->retrieve(m_SCT_ID, "SCT_ID").isFailure()) {
    return StatusCode::FAILURE;
  }

  if (detStore()->retrieve(m_SCT_Manager, "SCT").isFailure()) {
    return StatusCode::FAILURE;
  }

  // particle property service
  if (m_particlePropSvc.retrieve().isFailure())
  {
    ATH_MSG_ERROR("Can not retrieve " << m_particlePropSvc << " . Aborting ... " );
    return StatusCode::FAILURE;
  }

   // and the particle data table 
  m_particleDataTable = m_particlePropSvc->PDT();
  if (m_particleDataTable==0)
  {
    ATH_MSG_ERROR( "Could not get ParticleDataTable! Cannot associate pdg code with charge. Aborting. " );
    return StatusCode::FAILURE;
  }

  // Define the TTree
  //
  ITHistSvc *tHistSvc;
  StatusCode sc =  service("THistSvc", tHistSvc);
  if ( sc.isFailure() ) {
    ATH_MSG_ERROR( "Unable to retrieve pointer to THistSvc" );
    return sc;
  }
  m_nt = new TTree(TString(m_ntupleTreeName), "Athena Dump for GNN4ITk");
  // NB: we must not delete the tree, this is done by THistSvc
  std::string fullNtupleName =  m_ntupleFileName + m_ntupleDirName + m_ntupleTreeName;
  sc = tHistSvc->regTree(fullNtupleName, m_nt);
  if (sc.isFailure()) {
    ATH_MSG_ERROR( "Unable to register TTree: " << fullNtupleName );
    return sc;
  }

  if (m_csvFile) std::cout << "saving files in csv format" << std::endl;
  if (m_rootFile) std::cout << "saving files in root format" << std::endl;

  if (m_rootFile) {
    m_CLID = new int[m_maxCL];

    m_CLindex = new int[m_maxCL];
    m_CLhardware = new std::vector<std::string>;
    m_CLx = new double[m_maxCL];
    m_CLy = new double[m_maxCL];
    m_CLz = new double[m_maxCL];
    m_CLbarrel_endcap = new int[m_maxCL];
    m_CLlayer_disk = new int[m_maxCL];
    m_CLeta_module = new int[m_maxCL];
    m_CLphi_module = new int[m_maxCL];
    m_CLside = new int[m_maxCL];
    m_CLmoduleID = new uint64_t[m_maxCL];
    m_CLparticleLink_eventIndex = new std::vector<std::vector<int>>;
    m_CLparticleLink_barcode = new std::vector<std::vector<int>>;
    m_CLbarcodesLinked = new std::vector<std::vector<bool>>;
    m_CLphis = new std::vector<std::vector<int>>;
    m_CLetas = new std::vector<std::vector<int>>;
    m_CLtots = new std::vector<std::vector<int>>;
    m_CLloc_direction1 = new double[m_maxCL];
    m_CLloc_direction2 = new double[m_maxCL];
    m_CLloc_direction3 = new double[m_maxCL];
    m_CLJan_loc_direction1 = new double[m_maxCL];
    m_CLJan_loc_direction2 = new double[m_maxCL];
    m_CLJan_loc_direction3 = new double[m_maxCL];
    m_CLpixel_count = new int[m_maxCL];
    m_CLcharge_count = new float[m_maxCL];
    m_CLloc_eta = new float[m_maxCL];
    m_CLloc_phi = new float[m_maxCL];
    m_CLglob_eta = new float[m_maxCL];
    m_CLglob_phi = new float[m_maxCL];
    m_CLeta_angle = new double[m_maxCL];
    m_CLphi_angle = new double[m_maxCL];
    m_CLnorm_x = new float[m_maxCL];
    m_CLnorm_y = new float[m_maxCL];
    m_CLnorm_z = new float[m_maxCL];
    m_CLlocal_cov = new std::vector<std::vector<double>>;

    m_CLevent_number = new int[m_maxCL];
    m_CLbarcode = new int[m_maxCL];
    m_CLpx = new float[m_maxCL];
    m_CLpy = new float[m_maxCL];
    m_CLpz = new float[m_maxCL];
    m_CLpt = new float[m_maxCL];
    m_CLeta = new float[m_maxCL];
    m_CLvx = new float[m_maxCL];
    m_CLvy = new float[m_maxCL];
    m_CLvz = new float[m_maxCL];
    m_CLradius = new float[m_maxCL];
    m_CLstatus = new float[m_maxCL];
    m_CLcharge = new float[m_maxCL];
    m_CLpdg_id = new int[m_maxCL];
    m_CLpassed = new int[m_maxCL];
    m_CLvProdNin = new int[m_maxCL];
    m_CLvProdNout = new int[m_maxCL];
    m_CLvProdStatus = new int[m_maxCL];
    m_CLvProdBarcode = new int[m_maxCL];
    m_CLvParentID = new std::vector<std::vector<int>>;
    m_CLvParentBarcode = new std::vector<std::vector<int>>;

    m_SPindex = new int[m_maxCL];
    m_SPx = new double[m_maxCL];
    m_SPy = new double[m_maxCL];
    m_SPz = new double[m_maxCL];
    m_SPCL1_index = new int[m_maxCL];
    m_SPCL2_index = new int[m_maxCL];

    m_TRKindex = new int[m_maxCL];
    m_TRKtrack_fitter = new int[m_maxCL];
    m_TRKparticle_hypothesis = new int[m_maxCL];
    m_TRKproperties = new std::vector<std::vector<int>>;
    m_TRKpattern = new std::vector<std::vector<int>>;
    m_TRKndof = new int[m_maxCL];
    m_TRKmot = new int[m_maxCL];
    m_TRKoot = new int[m_maxCL];
    m_TRKchiSq = new float[m_maxCL];
    m_TRKmeasurementsOnTrack_pixcl_sctcl_index = new std::vector<std::vector<int>>;
    m_TRKoutliersOnTrack_pixcl_sctcl_index = new std::vector<std:: vector<int>>;
    m_TRKcharge = new int [m_maxCL];
    m_TRKperigee_position = new std::vector<std::vector<double>>;
    m_TRKperigee_momentum = new std::vector<std::vector<double>>;
    m_TTCindex = new int[m_maxCL];
    m_TTCevent_index = new int[m_maxCL];
    m_TTCparticle_link = new int[m_maxCL];
    m_TTCprobability = new float[m_maxCL];

    m_DTTindex = new int[m_maxCL];
    m_DTTsize = new int[m_maxCL];
    m_DTTtrajectory_eventindex = new std::vector<std::vector<int>>;
    m_DTTtrajectory_barcode = new std::vector<std::vector<int>>;
    m_DTTstTruth_subDetType = new std::vector<std::vector<int>>;
    m_DTTstTrack_subDetType = new std::vector<std::vector<int>>;
    m_DTTstCommon_subDetType = new std::vector<std::vector<int>>;

    m_nt->Branch("nCL_ID",&m_nCL_ID,"nCL_ID/I");
    m_nt->Branch("CLID", m_CLID,"CLID[nCL_ID]/I");

    m_nt->Branch("nCL",&m_nCL,"nCL/I");
    m_nt->Branch("CLindex",m_CLindex,"CLindex[nCL]/I");
    m_nt->Branch("CLhardware", &m_CLhardware);
    m_nt->Branch("CLx",m_CLx,"CLx[nCL]/D");
    m_nt->Branch("CLy",m_CLy,"CLy[nCL]/D");
    m_nt->Branch("CLz",m_CLz,"CLz[nCL]/D");
    m_nt->Branch("CLbarrel_endcap", m_CLbarrel_endcap,"CLbarrel_endcap[nCL]/I");
    m_nt->Branch("CLlayer_disk", m_CLlayer_disk,"CLlayer_disk[nCL]/I");
    m_nt->Branch("CLeta_module", m_CLeta_module,"CLeta_module[nCL]/I");
    m_nt->Branch("CLphi_module", m_CLphi_module,"CLphi_module[nCL]/I");
    m_nt->Branch("CLside", m_CLside,"CLside[nCL]/I");
    m_nt->Branch("CLmoduleID", m_CLmoduleID);//,"CLmoduleID[nCL]/g");
    m_nt->Branch("CLparticleLink_eventIndex", &m_CLparticleLink_eventIndex);
    m_nt->Branch("CLparticleLink_barcode", &m_CLparticleLink_barcode);
    m_nt->Branch("CLbarcodesLinked", &m_CLbarcodesLinked);
    m_nt->Branch("CLphis", &m_CLphis);
    m_nt->Branch("CLetas", &m_CLetas);
    m_nt->Branch("CLtots", &m_CLtots);
    m_nt->Branch("CLparticleLink_barcode", &m_CLparticleLink_barcode);
    m_nt->Branch("CLbarcodesLinked", &m_CLbarcodesLinked);
    m_nt->Branch("CLphis", &m_CLphis);
    m_nt->Branch("CLetas", &m_CLetas);
    m_nt->Branch("CLtots", &m_CLtots);
    m_nt->Branch("CLloc_direction1", m_CLloc_direction1,"CLloc_direction1[nCL]/D");
    m_nt->Branch("CLloc_direction2", m_CLloc_direction2,"CLloc_direction2[nCL]/D");
    m_nt->Branch("CLloc_direction3", m_CLloc_direction3,"CLloc_direction3[nCL]/D");
    m_nt->Branch("CLJan_loc_direction1", m_CLJan_loc_direction1,"CLJan_loc_direction1[nCL]/D");
    m_nt->Branch("CLJan_loc_direction2", m_CLJan_loc_direction2,"CLJan_loc_direction2[nCL]/D");
    m_nt->Branch("CLJan_loc_direction3", m_CLJan_loc_direction3,"CLJan_loc_direction3[nCL]/D");
    m_nt->Branch("CLpixel_count", m_CLpixel_count,"CLpixel_count[nCL]/I");
    m_nt->Branch("CLcharge_count", m_CLcharge_count,"CLcharge_count[nCL]/F");
    m_nt->Branch("CLloc_eta", m_CLloc_eta,"CLloc_eta[nCL]/F");
    m_nt->Branch("CLloc_phi", m_CLloc_phi,"CLloc_phi[nCL]/F");
    m_nt->Branch("CLglob_eta", m_CLglob_eta,"CLglob_eta[nCL]/F");
    m_nt->Branch("CLglob_phi", m_CLglob_phi,"CLglob_phi[nCL]/F");
    m_nt->Branch("CLeta_angle", m_CLeta_angle,"CLeta_angle[nCL]/D");
    m_nt->Branch("CLphi_angle", m_CLphi_angle,"CLphi_angle[nCL]/D");
    m_nt->Branch("CLnorm_x", m_CLnorm_x,"CLnorm_x[nCL]/F");
    m_nt->Branch("CLnorm_y", m_CLnorm_y,"CLnorm_y[nCL]/F");
    m_nt->Branch("CLnorm_z", m_CLnorm_z,"CLnorm_z[nCL]/F");
    m_nt->Branch("CLlocal_cov", &m_CLlocal_cov);

    m_nt->Branch("nPartEVT",&m_nPartEVT,"nPartEVT/I");
    m_nt->Branch("CLevent_number",m_CLevent_number,"CLevent_number[nPartEVT]/I");
    m_nt->Branch("CLbarcode",m_CLbarcode,"CLbarcode[nPartEVT]/I");
    m_nt->Branch("CLpx",m_CLpx,"CLpx[nPartEVT]/F");
    m_nt->Branch("CLpy",m_CLpy,"CLpy[nPartEVT]/F");
    m_nt->Branch("CLpz",m_CLpz,"CLpz[nPartEVT]/F");
    m_nt->Branch("CLpt",m_CLpt,"CLpt[nPartEVT]/F");
    m_nt->Branch("CLeta",m_CLeta,"CLeta[nPartEVT]/F");
    m_nt->Branch("CLvx",m_CLvx,"CLvx[nPartEVT]/F");
    m_nt->Branch("CLvy",m_CLvy,"CLvy[nPartEVT]/F");
    m_nt->Branch("CLvz",m_CLvz,"CLvz[nPartEVT]/F");
    m_nt->Branch("CLradius",m_CLradius,"CLradius[nPartEVT]/F");
    m_nt->Branch("CLstatus",m_CLstatus,"CLstatus[nPartEVT]/F");
    m_nt->Branch("CLcharge",m_CLcharge,"CLcharge[nPartEVT]/F");
    m_nt->Branch("CLpdg_id",m_CLpdg_id,"CLpdg_id[nPartEVT]/I");
    m_nt->Branch("CLpassed",m_CLpassed,"CLpassed[nPartEVT]/I");
    m_nt->Branch("CLvProdNin",m_CLvProdNin,"CLvProdNin[nPartEVT]/I");
    m_nt->Branch("CLvProdNout",m_CLvProdNout,"CLvProdNout[nPartEVT]/I");
    m_nt->Branch("CLvProdStatus",m_CLvProdStatus,"CLvProdStatus[nPartEVT]/I");
    m_nt->Branch("CLvProdBarcode", m_CLvProdBarcode,"CLvProdBarcode[nPartEVT]/I");
    m_nt->Branch("CLvParentID", &m_CLvParentID);
    m_nt->Branch("CLvParentBarcode", &m_CLvParentBarcode);

    m_nt->Branch("nSP",&m_nSP,"nSP/I");
    m_nt->Branch("SPindex",m_SPindex,"SPindex[nSP]/I");
    m_nt->Branch("SPx",m_SPx,"SPx[nSP]/D");
    m_nt->Branch("SPy",m_SPy,"SPy[nSP]/D");
    m_nt->Branch("SPz",m_SPz,"SP_z[nSP]/D");
    m_nt->Branch("SPCL1_index",m_SPCL1_index,"SPCL1_index[nSP]/I");
    m_nt->Branch("SPCL2_index",m_SPCL2_index,"NPCL2_index[nSP]/I");

    m_nt->Branch("nTRK",&m_nTRK,"nTRK/I");
    m_nt->Branch("TRKindex",m_TRKindex,"TRKindex[nTRK]/I");
    m_nt->Branch("TRKtrack_fitter",m_TRKtrack_fitter,"TRKtrack_fitter[nTRK]/I");
    m_nt->Branch("TRKparticle_hypothesis",m_TRKparticle_hypothesis,"TRKparticle_hypothesis[nTRK]/I");
    m_nt->Branch("TRKproperties",&m_TRKproperties);
    m_nt->Branch("TRKpattern",&m_TRKpattern);
    m_nt->Branch("TRKndof",m_TRKndof,"TRKndof[nTRK]/I");
    m_nt->Branch("TRKmot",m_TRKmot,"TRKmot[nTRK]/I");
    m_nt->Branch("TRKoot",m_TRKoot,"TRKoot[nTRK]/I");
    m_nt->Branch("TRKchiSq",m_TRKchiSq,"TRKchiSq[nTRK]/F");
    m_nt->Branch("TRKmeasurementsOnTrack_pixcl_sctcl_index", &m_TRKmeasurementsOnTrack_pixcl_sctcl_index);
    m_nt->Branch("TRKoutliersOnTrack_pixcl_sctcl_index", &m_TRKoutliersOnTrack_pixcl_sctcl_index);
    m_nt->Branch("TRKcharge",m_TRKcharge,"TRKcharge[nTRK]/I");
    m_nt->Branch("TRKperigee_position", &m_TRKperigee_position);
    m_nt->Branch("TRKperigee_momentum", &m_TRKperigee_momentum);
    m_nt->Branch("TTCindex",m_TTCindex,"TTCindex[nTRK]/I");
    m_nt->Branch("TTCevent_index",m_TTCevent_index,"TTCevent_index[nTRK]/I");
    m_nt->Branch("TTCparticle_link",m_TTCparticle_link,"TTCparticle_link[nTRK]/I");
    m_nt->Branch("TTCprobability",m_TTCprobability,"TTCprobability[nTRK]/F");

    m_nt->Branch("nDTT",&m_nDTT,"nDTT/I");
    m_nt->Branch("DTTindex",m_DTTindex,"DTTindex[nDTT]/I");
    m_nt->Branch("DTTsize",m_DTTsize,"DTTsize[nDTT]/I");
    m_nt->Branch("DTTtrajectory_eventindex", &m_DTTtrajectory_eventindex);
    m_nt->Branch("DTTtrajectory_barcode", &m_DTTtrajectory_barcode);
    m_nt->Branch("DTTstTruth_subDetType", &m_DTTstTruth_subDetType);
    m_nt->Branch("DTTstTrack_subDetType", &m_DTTstTrack_subDetType);
    m_nt->Branch("DTTstCommon_subDetType", &m_DTTstCommon_subDetType);
  }

  return StatusCode::SUCCESS;
}


//-------------------------------
StatusCode DumpObjects::execute() {
//-------------------------------
  m_selected = 0;
  m_event++;

  // create a container with HepMcParticleLink and list of clusters
  // particle barcode --> is accepted and number of clusters
  std::map < std::pair < int, int >, std::pair < bool, int > > allTruthParticles;

  const McEventCollection* mcCollptr = 0;
  ATH_CHECK (evtStore()->retrieve(mcCollptr, "TruthEvent"));

  int accepted = 0;
  int all_truth = 0;

  // dump out file for truth events

  bool doSaveFile = false;
  std::string filename = (m_name == "") ? "subevents_evt" : (m_name+"_subevents_evt");

  std::map < int, int > allSubEvents;

  m_nCL_ID = 0;

  bool duplicateSubeventID = false;
  for (unsigned int cntr = 0; cntr < mcCollptr->size(); ++cntr) {
    if (not doSaveFile) {
      if (m_csvFile)
        m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
      doSaveFile = true;
    }
    int ID = mcCollptr->at(cntr)->event_number();
    if (doSaveFile) {
      if (m_csvFile)
        *m_outfile<<ID<<std::endl;
      if (m_rootFile)
        m_CLID[m_nCL_ID++]=ID;

  	  if (m_nCL_ID==m_maxCL) {
        std::cout<<"===================================="<<std::endl;
        std::cout<<"DUMP : hit max number of subevent ID"<<std::endl;
        std::cout<<"===================================="<<std::endl;
    	  break;
      }
    }

    std::map < int, int >::iterator it = allSubEvents.find(ID);
    if (it==allSubEvents.end()) allSubEvents.insert(std::make_pair(ID,1));
    else {
      it->second++;
      duplicateSubeventID=true;
    }
  }

  if (duplicateSubeventID) {
    std::cout<<"Duplicate subevent ID in event "<<m_event<<"."<<std::endl;
    if (doSaveFile && m_csvFile)
      *m_outfile << "Duplicate subevent ID in event " << m_event << "." << std::endl;
  }

  if (doSaveFile && m_csvFile) m_outfile->close();

  // dump out file for truth particles
  doSaveFile = false;
  if (m_csvFile)
    filename = (m_name == "") ? "particles_evt" : (m_name+"_particles_evt");

  m_nPartEVT=0;

  if (m_rootFile) {
    (*m_CLvParentID).clear();
    (*m_CLvParentBarcode).clear();
  }

  for (unsigned int cntr = 0; cntr < mcCollptr->size(); ++cntr) {
    const HepMC::GenEvent* genEvt = (mcCollptr->at(cntr));
    for ( HepMC::GenEvent::particle_const_iterator p = genEvt->particles_begin(); p != genEvt->particles_end(); ++p ) {
      if (not doSaveFile) {
        if (m_csvFile)
          m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
        doSaveFile = true;
      }

      //*p is a GenParticle 
      float px, py, pz, pt, eta, vx, vy, vz, radius, status, charge = 0.;
      std::vector<int> vParentID;
      std::vector<int> vParentBarcode;

      int vProdNin, vProdNout, vProdStatus, vProdBarcode;
      bool passed = isPassed(*p, px, py, pz, pt, eta, vx, vy, vz, radius, status, charge, vParentID, vParentBarcode, vProdNin, vProdNout,  vProdStatus, vProdBarcode);
      allTruthParticles.insert(std::make_pair(std::make_pair(genEvt->event_number(), (*p)->barcode()), std::make_pair(passed, 0)));
      all_truth++;
      if (passed) accepted++;   
      if (doSaveFile) {
        // subevent, barcode, px, py, pz, pt, eta, vx, vy, vz, radius, status, charge
        if (m_csvFile) {
          *m_outfile<<genEvt->event_number()<<","<<(*p)->barcode()<<","<<px<<","<<py<<","<<pz<<","<<pt<<","<<eta<<","<<vx<<","<<vy<<","<<vz<<","<<radius<<","<<status<<","<<charge<<","<<(*p)->pdg_id()<<","<<(passed ? "YES" : "NO")<<",";
        	*m_outfile<<vProdNin<<","<<vProdNout<<","<<vProdStatus<<","<<vProdBarcode<<",#,";
          for(std::size_t ipar=0; ipar<vParentID.size(); ipar++)
	          *m_outfile<<"("<<vParentID[ipar]<<","<<vParentBarcode[ipar]<<"),";
      	  *m_outfile<<"#"<<std::endl;
        }
        if (m_rootFile) {
  	      m_CLevent_number[m_nPartEVT]=genEvt->event_number();
	        m_CLbarcode[m_nPartEVT]=(*p)->barcode();
	        m_CLpx[m_nPartEVT]=px;
	        m_CLpy[m_nPartEVT]=py;
	        m_CLpz[m_nPartEVT]=pz;
          m_CLpt[m_nPartEVT]=pt;
          m_CLeta[m_nPartEVT]=eta;
          m_CLvx[m_nPartEVT] = vx;
          m_CLvy[m_nPartEVT] = vy;
          m_CLvz[m_nPartEVT] = vz;
          m_CLradius[m_nPartEVT] = radius;
          m_CLstatus[m_nPartEVT] = status;
          m_CLcharge[m_nPartEVT] = charge;
          m_CLpdg_id[m_nPartEVT] = (*p)->pdg_id();
          m_CLpassed[m_nPartEVT] = (passed ? true : false);
          m_CLvProdNin[m_nPartEVT] = vProdNin;
          m_CLvProdNout[m_nPartEVT] = vProdNout;
          m_CLvProdStatus[m_nPartEVT] = vProdStatus;
          m_CLvProdBarcode[m_nPartEVT] = vProdBarcode;
          (*m_CLvParentID).push_back(vParentID);
          (*m_CLvParentBarcode).push_back(vParentBarcode);
        }

    	  m_nPartEVT++;
    	  if (m_nPartEVT==m_maxCL) {
          std::cout<<"======================================="<<std::endl;
         std::cout<<"DUMP : hit max number of particle events"<<std::endl;
          std::cout<<"======================================="<<std::endl;
    	    break;
        }
      }
    }
  }

  if (doSaveFile && m_csvFile) m_outfile->close();

  //   std::cout << "Having " << all_truth << " truth particles and accepting " << accepted << std::endl;

  const InDet::PixelClusterContainer* PixelClusterContainer = 0;     
  if( evtStore()->retrieve(PixelClusterContainer,"PixelClusters").isFailure()) {
    ATH_MSG_ERROR("Cannot retrieve Pixel PrepDataContainer PixelClusters");
    return StatusCode::FAILURE;
  }

  const InDet::SCT_ClusterContainer* SCT_ClusterContainer = 0;     
  if (evtStore()->retrieve(SCT_ClusterContainer,"SCT_Clusters").isFailure()) {
    ATH_MSG_ERROR("Cannot retrieve SCT PrepDataContainer SCT_Clusters");
    return StatusCode::FAILURE;
  }

  doSaveFile = false;
  if (m_csvFile)
    filename = (m_name == "") ? "clusters_evt" : (m_name+"_clusters_evt");  
  auto cartesion_to_spherical = [](const Amg::Vector3D &xyzVec, float &eta_, float &phi_) {
    float r3 = 0;
    for (int idx = 0; idx < 3; ++idx) {
      r3 += xyzVec[idx] * xyzVec[idx];
    }
    r3 = sqrt(r3);
    phi_ = atan2(xyzVec[1], xyzVec[0]);
    float theta_ = acos(xyzVec[2] / r3);
    eta_ = log(tan(0.5 * theta_));
  };


  ///////////////////////////////////////////////////////////////////////////
  //////////////////////////////PIXEL CONTAINER//////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  m_nCL=0;
  if (m_rootFile) {
    (*m_CLhardware).clear();
    (*m_CLparticleLink_eventIndex).clear();
    (*m_CLparticleLink_barcode).clear();
    (*m_CLbarcodesLinked).clear();
    (*m_CLphis).clear();
    (*m_CLetas).clear();
    (*m_CLtots).clear();
    (*m_CLlocal_cov).clear();
  }

  if (PixelClusterContainer->size()>0) {

    const InDetSimDataCollection* sdoCollection = 0;
    if (evtStore()->retrieve(sdoCollection, "PixelSDO_Map").isFailure() ){
      ATH_MSG_ERROR( "Could not find the data object PixelSDO_Map !" );
      return StatusCode::FAILURE;
    } 

    for (const auto& clusterCollection : *PixelClusterContainer){
      // skip empty collections
      if (clusterCollection->empty()) continue;

      if (not doSaveFile) {
        if (m_csvFile)
          m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
        doSaveFile = true;
      }

      int barrel_endcap = m_pixelID->barrel_ec(clusterCollection->identify());
      int layer_disk    = m_pixelID->layer_disk(clusterCollection->identify());
      int eta_module    = m_pixelID->eta_module(clusterCollection->identify());
      int phi_module    = m_pixelID->phi_module(clusterCollection->identify());

      InDetDD::SiDetectorElement* element = m_pixelManager->getDetectorElement(clusterCollection->identify());

      std::string ModulesInfo("Modules_geo.txt");
      std::fstream outfile;
      outfile.open(ModulesInfo, std::ios_base::app);
      outfile<< "PIXEL"<<" "<< barrel_endcap <<" " << layer_disk << " " << eta_module << " " << phi_module << " " 
                                             <<element->center().z() <<" "
                                             <<element->center().x()<<" "
                                             <<element->center().y()<<" "
                                             <<clusterCollection->identify().get_compact()<<" "
                                             <<0
                                             <<std::endl;
      outfile.close();

      Amg::Vector3D my_normal = element->normal();
      float norm_x = fabs(my_normal.x())>1e-5 ? my_normal.x() : 0.;
      float norm_y = fabs(my_normal.y())>1e-5 ? my_normal.y() : 0.;
      float norm_z = fabs(my_normal.z())>1e-5 ? my_normal.z() : 0.;          

      const InDetDD::PixelModuleDesign* design (dynamic_cast<const InDetDD::PixelModuleDesign*>(&element->design()));

      if (not design) {
        ATH_MSG_ERROR("Dynamic cast failed at "<<__LINE__<<" of MergedPixelsTool.cxx.");
        return StatusCode::FAILURE;
      }

      // loop over collection
      for (InDet::PixelCluster* const & cluster : *clusterCollection ) {
        Identifier clusterId = cluster->identify();
        if ( !clusterId.is_valid() ) {
          ATH_MSG_WARNING("Pixel cluster identifier is not valid");
        }

        const Amg::MatrixX& local_cov = cluster->localCovariance();

        std::vector < std::pair<int,int> > barcodes = {};
        std::vector<int> particleLink_eventIndex = {};
        std::vector<int> particleLink_barcode = {};
	      std::vector<bool> barcodesLinked = {};
        std::vector<int> phis = {};
        std::vector<int> etas = {};
        std::vector<int> tots = {};
        int min_eta = 999;
        int min_phi = 999;
        int max_eta = -999;
        int max_phi = -999;

        float charge_count = 0;
        int pixel_count = 0;

        for (unsigned int rdo = 0; rdo<cluster->rdoList().size(); rdo++ ){
          const auto& rdoID = cluster->rdoList().at(rdo);
          int phi = m_pixelID->phi_index(rdoID);
          int eta = m_pixelID->eta_index(rdoID);
          if(min_eta > eta) min_eta = eta;
          if(min_phi > phi) min_phi = phi;
          if(max_eta < eta) max_eta = eta;
          if(max_phi < phi) max_phi = phi;

          ++pixel_count;
          charge_count += cluster->totList().at(rdo);

          phis.push_back(phi);
          etas.push_back(eta);
          tots.push_back(cluster->totList().at(rdo));

          auto pos = sdoCollection->find(rdoID);
          if (pos != sdoCollection->end()) {
            for (auto deposit : pos->second.getdeposits()){
              const HepMcParticleLink& particleLink = deposit.first;
    	        std::pair<int,int> barcode(particleLink.eventIndex(),particleLink.barcode());
              if (particleLink.isValid()) allTruthParticles.at(barcode).second++;
              if (std::find(barcodes.begin(), barcodes.end(),barcode)==barcodes.end()) {
                barcodes.push_back(barcode);
                particleLink_eventIndex.push_back(particleLink.eventIndex());
                particleLink_barcode.push_back(particleLink.barcode());
            		barcodesLinked.push_back(particleLink.isValid());
              }
            }
          }
        }

      	InDetDD::SiLocalPosition localPos_entry = design->localPositionOfCell(InDetDD::SiCellId(min_phi,min_eta));
	      InDetDD::SiLocalPosition localPos_exit  = design->localPositionOfCell(InDetDD::SiCellId(max_phi,max_eta));

        Amg::Vector3D localStartPosition(localPos_entry.xEta()-0.5*element->etaPitch(), localPos_entry.xPhi()-0.5*element->phiPitch(), -0.5*element->thickness());
        Amg::Vector3D localEndPosition  (localPos_exit.xEta() +0.5*element->etaPitch(), localPos_exit.xPhi() +0.5*element->phiPitch(),  0.5*element->thickness());

        // local direction in local coordinates
        // clusterShape: [lx, ly, lz]
        Amg::Vector3D localDirection = localEndPosition - localStartPosition;

        float loc_eta = 0, loc_phi = 0; // clusterShape: [leta, lphi]
        cartesion_to_spherical(localDirection, loc_eta, loc_phi);

        Amg::Vector3D globalStartPosition = element->globalPosition(localStartPosition);
        Amg::Vector3D globalEndPosition   = element->globalPosition(localEndPosition  );

        Amg::Vector3D direction = globalEndPosition-globalStartPosition;
        float glob_eta = 0, glob_phi = 0; // clusterShape: [geta, gphi]
        cartesion_to_spherical(direction, glob_eta, glob_phi); 

        Amg::Vector3D my_phiax = element->phiAxis();
        Amg::Vector3D my_etaax = element->etaAxis();

        float trkphicomp = direction.dot(my_phiax);
        float trketacomp = direction.dot(my_etaax);
        float trknormcomp = direction.dot(my_normal);
        double phi_angle = atan2(trknormcomp,trkphicomp);
        double eta_angle = atan2(trknormcomp,trketacomp);
        // now dumping all the values now
      	cluster->setIndex(m_selected++);
	      if (doSaveFile) {
	        // fill TXT file
          if (m_csvFile) {
            *m_outfile << m_selected <<",PIXEL,"<< cluster->globalPosition().x() << "," << cluster->globalPosition().y() << "," << cluster->globalPosition().z() << ",#," << barrel_endcap <<"," << layer_disk << "," << eta_module << "," << phi_module << ",0" << ",#,";
            for (unsigned truth = 0 ; truth<barcodes.size(); truth++)
              *m_outfile<<"("<<barcodes.at(truth).first<<","<<barcodes.at(truth).second<<","<<barcodesLinked.at(truth)<<"),";
            *m_outfile<<"#,";
            for (unsigned pixel = 0 ; pixel<tots.size(); pixel++)
              *m_outfile << "(" << etas.at(pixel) << "," << phis.at(pixel) << "," << tots.at(pixel) << "),";
            *m_outfile << "#," << pixel_count << "," << charge_count
                               << "," << loc_eta << "," << loc_phi
                               << "," << localDirection[0] << "," << localDirection[1]<< "," << localDirection[2]
                               << ","<<0<<","<<0<<","<<0
                               <<","<< glob_eta << "," << glob_phi
                               << "," << eta_angle << "," << phi_angle
                               << ",#," << norm_x << "," << norm_y << "," << norm_z;

            *m_outfile<<",#";
          }
          std::vector<double> v_local_cov;
          if (local_cov.size() > 0) {
            for (size_t i=0, nRows = local_cov.rows(), nCols = local_cov.cols(); i < nRows; i++) {
              for (size_t j=0; j < nCols; ++j){
                if (m_csvFile)
                  *m_outfile << "," << local_cov(i, j);
                v_local_cov.push_back(local_cov(i,j));
              }
            }
          }
          else if(m_csvFile) *m_outfile << ",0";
          if(m_csvFile) *m_outfile << std::endl;

          if (m_rootFile){
        	  // fill TTree
	          m_CLindex[m_nCL]=m_selected;
	          (*m_CLhardware).push_back("PIXEL");
    	      m_CLx[m_nCL]=cluster->globalPosition().x();
      	    m_CLy[m_nCL]=cluster->globalPosition().y();
      	    m_CLz[m_nCL]=cluster->globalPosition().z();
            m_CLbarrel_endcap[m_nCL] = barrel_endcap;
            m_CLlayer_disk[m_nCL] = layer_disk;
            m_CLeta_module[m_nCL] = eta_module;
            m_CLphi_module[m_nCL] = phi_module;
            m_CLside[m_nCL] = 0;
            m_CLmoduleID[m_nCL] = clusterCollection->identify().get_compact();
            (*m_CLparticleLink_eventIndex).push_back(particleLink_eventIndex);
            (*m_CLparticleLink_barcode).push_back(particleLink_barcode);
            (*m_CLbarcodesLinked).push_back(barcodesLinked);
            (*m_CLetas).push_back(etas);
            (*m_CLphis).push_back(phis);
            (*m_CLtots).push_back(tots);
            m_CLloc_direction1[m_nCL] = localDirection[0];
            m_CLloc_direction2[m_nCL] = localDirection[1];
            m_CLloc_direction3[m_nCL] = localDirection[2];
            m_CLJan_loc_direction1[m_nCL] = 0;
            m_CLJan_loc_direction2[m_nCL] = 0;
            m_CLJan_loc_direction3[m_nCL] = 0;
            m_CLpixel_count[m_nCL] = pixel_count;
            m_CLcharge_count[m_nCL] = charge_count;
            m_CLloc_eta[m_nCL] = loc_eta;
            m_CLloc_phi[m_nCL] = loc_phi;
            m_CLglob_eta[m_nCL] = glob_eta;
            m_CLglob_phi[m_nCL] = glob_phi;
            m_CLeta_angle[m_nCL] = eta_angle;
            m_CLphi_angle[m_nCL] = phi_angle;
            m_CLnorm_x[m_nCL] = norm_x;
            m_CLnorm_y[m_nCL] = norm_y;
            m_CLnorm_z[m_nCL] = norm_z;
            (*m_CLlocal_cov).push_back(v_local_cov);
          }
          m_nCL++;
      	  if (m_nCL==m_maxCL) {
            std::cout<<"================================="<<std::endl;
	          std::cout<<"DUMP : hit max number of clusters"<<std::endl;
            std::cout<<"================================="<<std::endl;
      	    break;
	        }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////////SCT CONTAINER//////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  if (SCT_ClusterContainer->size()>0) {
    const InDetSimDataCollection* sdoCollection = 0;
    if (evtStore()->retrieve(sdoCollection, "SCT_SDO_Map").isFailure()) {
      ATH_MSG_ERROR ("Could not find the data object SCT_SDO_Map !");
      return StatusCode::FAILURE;
    } 

    for (const auto& clusterCollection : *SCT_ClusterContainer){
      // skip empty collections
      if (clusterCollection->empty()) continue;

      if (not doSaveFile) {
        if(m_csvFile)
          m_outfile = new std::ofstream(("./cluster_file_evt"+std::to_string(m_event)+".txt").c_str());
        doSaveFile = true;
      }

      int barrel_endcap = m_SCT_ID->barrel_ec(clusterCollection->identify());
      int layer_disk    = m_SCT_ID->layer_disk(clusterCollection->identify());
      int eta_module    = m_SCT_ID->eta_module(clusterCollection->identify());
      int phi_module    = m_SCT_ID->phi_module(clusterCollection->identify());
      int side          = m_SCT_ID->side(clusterCollection->identify());

      InDetDD::SiDetectorElement* element = m_SCT_Manager->getDetectorElement(clusterCollection->identify());

      std::string ModulesInfo("Modules_geo.txt");
      std::fstream outfile;
      outfile.open(ModulesInfo, std::ios_base::app);

      outfile<< "STRIP"<<" "<< barrel_endcap <<" " << layer_disk << " " << eta_module << " " << phi_module << " " 
                                             <<element->center().z() << " "
                                             <<element->center().x() << " "
                                             <<element->center().y() << " "
                                             <<clusterCollection->identify().get_compact() << " "
                                             << side
                                             << std::endl;
      outfile.close();

      Amg::Vector3D my_normal = element->normal();
      float norm_x = fabs(my_normal.x())>1e-5 ? my_normal.x() : 0.;
      float norm_y = fabs(my_normal.y())>1e-5 ? my_normal.y() : 0.;
      float norm_z = fabs(my_normal.z())>1e-5 ? my_normal.z() : 0.;

      // loop over collection
      for (InDet::SCT_Cluster* const & cluster : *clusterCollection) {
        Identifier clusterId = cluster->identify();
        if ( !clusterId.is_valid() ) {
          ATH_MSG_WARNING("SCT cluster identifier is not valid");
        }

        const Amg::MatrixX& local_cov = cluster->localCovariance();

        std::vector < std::pair<int, int> > barcodes = {};
        std::vector<int> particleLink_eventIndex = {};
        std::vector<int> particleLink_barcode = {};
        std::vector < bool > barcodesLinked = {};

        std::vector<int> tots = {};
        std::vector<int> strip_ids = {};
        int min_strip = 999;
        int max_strip = -999;

        float charge_count = 0;
        int pixel_count = 0;

        for (unsigned int rdo = 0; rdo<cluster->rdoList().size(); rdo++) {
          const auto& rdoID = cluster->rdoList().at(rdo);

          int strip = m_SCT_ID->strip(rdoID);

          if (min_strip > strip) min_strip = strip;
          if (max_strip < strip) max_strip = strip;
          strip_ids.push_back(strip);
          // tots.push_back(cluster->totList().at(rdo));
          tots.push_back(0);  // FIXME
          ++pixel_count;
          // find barcodes of the truth particles
          auto pos = sdoCollection->find(rdoID);
          if (pos != sdoCollection->end()) {
            for (auto deposit : pos->second.getdeposits()) {
              const HepMcParticleLink& particleLink = deposit.first;
              std::pair<int,int> barcode(particleLink.eventIndex(),particleLink.barcode());
              // note that we are not filling the map allTruthParticles here - OK, we are not using this map for anything
              if (std::find(barcodes.begin(), barcodes.end(),barcode)==barcodes.end()) {
                barcodes.push_back(barcode);
                particleLink_eventIndex.push_back(particleLink.eventIndex());
                particleLink_barcode.push_back(particleLink.barcode());
                barcodesLinked.push_back(particleLink.isValid());
              }
            }
          }
        }

        // retrieve cluster shape
        const InDetDD::SCT_ModuleSideDesign* design (dynamic_cast<const InDetDD::SCT_ModuleSideDesign*>(&element->design()));
        if (not design) {
          ATH_MSG_ERROR("Failed at "<< __LINE__ <<" of accessing SCT ModuleSide Design");
          return StatusCode::FAILURE;
        }

        Amg::Vector2D locpos = cluster->localPosition(); 
        std::pair<Amg::Vector3D, Amg::Vector3D > ends(element->endsOfStrip(InDetDD::SiLocalPosition(locpos.y(), locpos.x(), 0)));

        Amg::Vector3D JanDirection = ends.second - ends.first;
        //std::cout<<JanDirection[0]<<" "<<JanDirection[1]<<" "<<JanDirection[2]<<std::endl;

        InDetDD::SiLocalPosition localPos_entry = design->localPositionOfCell(InDetDD::SiCellId(min_strip));
        InDetDD::SiLocalPosition localPos_exit  = design->localPositionOfCell(InDetDD::SiCellId(max_strip));

        Amg::Vector3D localStartPosition(localPos_entry.xEta()-0.5*element->etaPitch(), localPos_entry.xPhi()-0.5*element->phiPitch(), -0.5*element->thickness());
        Amg::Vector3D localEndPosition  (localPos_exit.xEta() +0.5*element->etaPitch(), localPos_exit.xPhi() +0.5*element->phiPitch(),  0.5*element->thickness());

        Amg::Vector3D localDirection = localEndPosition - localStartPosition;
        float loc_eta = 0, loc_phi = 0; // clusterShape: [leta, lphi]
        cartesion_to_spherical(localDirection, loc_eta, loc_phi);

        Amg::Vector3D globalStartPosition = element->globalPosition(localStartPosition);
        Amg::Vector3D globalEndPosition   = element->globalPosition(localEndPosition  );

        Amg::Vector3D direction = globalEndPosition-globalStartPosition;
        float glob_eta = 0, glob_phi = 0; // clusterShape: [geta, gphi]
        cartesion_to_spherical(direction, glob_eta, glob_phi); 

        Amg::Vector3D my_phiax = element->phiAxis();
        Amg::Vector3D my_etaax = element->etaAxis();

        float trkphicomp = direction.dot(my_phiax);
        float trketacomp = direction.dot(my_etaax);
        float trknormcomp = direction.dot(my_normal);
        double phi_angle = atan2(trknormcomp,trkphicomp);
        double eta_angle = atan2(trknormcomp,trketacomp);

        // now dumping all the values now
        if (doSaveFile) {
          cluster->setIndex(m_selected++);
          if (m_csvFile) {
            *m_outfile << m_selected <<",STRIP,"<< cluster->globalPosition().x() << "," << cluster->globalPosition().y() << "," << cluster->globalPosition().z() << ",#," << barrel_endcap <<"," << layer_disk << "," << eta_module << "," << phi_module << "," << side << ",#,";
            for (unsigned truth = 0 ; truth<barcodes.size(); truth++)
              *m_outfile<<"("<<barcodes.at(truth).first<<","<<barcodes.at(truth).second<<","<<barcodesLinked.at(truth)<<"),";
            *m_outfile<<"#,";
          }
          // cluster shape
          std::vector<int> cst;
          for (unsigned strip=0; strip<strip_ids.size(); strip++) {
            if (m_csvFile)
              *m_outfile << "(" << strip_ids.at(strip) << ",-1," << tots.at(strip) << "),";
            cst.push_back(-1);
          }
          if (m_csvFile) {
            *m_outfile << "#," << pixel_count << "," << charge_count 
                               << "," << loc_eta << "," << loc_phi 
                               << "," << localDirection[0] << "," << localDirection[1]<< "," << localDirection[2] 
                               << ","<<JanDirection[0]<<","<<JanDirection[1]<<","<<JanDirection[2]
                               <<","<< glob_eta << "," << glob_phi
                               << "," << eta_angle << "," << phi_angle 
                               << ",#," << norm_x << "," << norm_y << "," << norm_z;

            //*m_outfile<<"#,"<<eta_angle<<","<<phi_angle<<",#,"<< norm_x << "," << norm_y << "," << norm_z;
            // local covariance matrix
            *m_outfile<<",#";
          }
          std::vector<double> v_local_cov;
          if (local_cov.size() > 0) {
            for(size_t i=0, nRows = local_cov.rows(), nCols = local_cov.cols(); i < nRows; i++) {
              for (size_t j=0; j < nCols; ++j) {
                if (m_csvFile)
                  *m_outfile << "," << local_cov(i, j);
                v_local_cov.push_back(local_cov(i,j));
              }
            }
          }
          else if (m_csvFile)
            *m_outfile << ",0";
          if (m_csvFile) *m_outfile << std::endl;

          if (m_rootFile) {
            m_CLindex[m_nCL]=m_selected;
           	(*m_CLhardware).push_back("STRIP");
//    std::cout << std::endl << boost::typeindex::type_id_with_cvr<decltype(cluster->globalPosition().x())>() << std::endl;
//    exit(0);
  	        m_CLx[m_nCL]=cluster->globalPosition().x();
          	m_CLy[m_nCL]=cluster->globalPosition().y();
          	m_CLz[m_nCL]=cluster->globalPosition().z();
            m_CLbarrel_endcap[m_nCL] = barrel_endcap;
            m_CLlayer_disk[m_nCL] = layer_disk;
            m_CLeta_module[m_nCL] = eta_module;
            m_CLphi_module[m_nCL] = phi_module;
            m_CLside[m_nCL] = side;
            m_CLmoduleID[m_nCL] = clusterCollection->identify().get_compact();
            (*m_CLparticleLink_eventIndex).push_back(particleLink_eventIndex);
            (*m_CLparticleLink_barcode).push_back(particleLink_barcode);
            (*m_CLbarcodesLinked).push_back(barcodesLinked);
            (*m_CLetas).push_back(strip_ids);
            (*m_CLphis).push_back(cst);
            (*m_CLtots).push_back(tots);
            m_CLloc_direction1[m_nCL] = localDirection[0];
            m_CLloc_direction2[m_nCL] = localDirection[1];
            m_CLloc_direction3[m_nCL] = localDirection[2];
            m_CLJan_loc_direction1[m_nCL] = JanDirection[0];
            m_CLJan_loc_direction2[m_nCL] = JanDirection[1];
            m_CLJan_loc_direction3[m_nCL] = JanDirection[2];
            m_CLpixel_count[m_nCL] = pixel_count;
            m_CLcharge_count[m_nCL] = charge_count;
            m_CLloc_eta[m_nCL] = loc_eta;
            m_CLloc_phi[m_nCL] = loc_phi;
            m_CLglob_eta[m_nCL] = glob_eta;
            m_CLglob_phi[m_nCL] = glob_phi;
            m_CLeta_angle[m_nCL] = eta_angle;
            m_CLphi_angle[m_nCL] = phi_angle;
            m_CLnorm_x[m_nCL] = norm_x;
            m_CLnorm_y[m_nCL] = norm_y;
            m_CLnorm_z[m_nCL] = norm_z;
            (*m_CLlocal_cov).push_back(v_local_cov);
          }

          m_nCL++;
      	  if (m_nCL==m_maxCL) {
            std::cout<<"================================="<<std::endl;
	          std::cout<<"DUMP : hit max number of clusters"<<std::endl;
            std::cout<<"================================="<<std::endl;
      	    break;
	        }
        }
      }
    }    
  }

  if (doSaveFile && m_csvFile) {
    m_outfile->close();
  }

  if (m_csvFile)
    filename = (m_name == "") ? "spacepoints_evt" : (m_name+"_spacepoints_evt");

  const SpacePointContainer* PixelSpacePointContainer = 0;     
  if (evtStore()->retrieve(PixelSpacePointContainer,"PixelSpacePoints").isFailure() ) {
      //ATH_MSG_ERROR("Cannot retrieve PixelSpacePoints");
      //return StatusCode::FAILURE;
  }

  const SpacePointContainer* SCT_SpacePointContainer = 0;     
  if (evtStore()->retrieve(SCT_SpacePointContainer,"SCT_SpacePoints").isFailure() ) {
      //ATH_MSG_ERROR("Cannot retrieve SCT_SpacePoints");
      //return StatusCode::FAILURE;
  }

  /*
  ////OVERLAP
  const SpacePointOverlapCollection*  overlapCollection(0);
  if( evtStore()->retrieve(overlapCollection,"OverlapSpacePoints").isFailure() ) {
      ATH_MSG_ERROR("Cannot retrieve OverlapSpacePoints");
      return StatusCode::FAILURE;
  }
  */

  doSaveFile = false;
  int sp_index = 0;

  m_nSP = 0;
  if (PixelSpacePointContainer && PixelSpacePointContainer->size()>0) {
    for( const auto& spCollection : *PixelSpacePointContainer){
      // skip empty collections
      if( spCollection->empty() ) continue;

      if (not doSaveFile) {
        if (m_csvFile)
          m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
        doSaveFile = true;
      }

      // loop over collection
      for( auto& sp : *spCollection ) {
        // save sp x, y, z and the index of the cluster associated to that one
        if (doSaveFile) {
          const InDet::SiCluster * cl = static_cast<const InDet::SiCluster*>(sp->clusterList().first);
          if (m_csvFile)
            *m_outfile<<sp_index++<<","<<sp->globalPosition().x()<<","<<sp->globalPosition().y()<<","<<sp->globalPosition().z()<<"," <<cl->getIndex()<<std::endl;
          if (m_rootFile) {
            m_SPindex[m_nSP] = sp_index-1;
            m_SPx[m_nSP] = sp->globalPosition().x();
            m_SPy[m_nSP] = sp->globalPosition().y();
            m_SPz[m_nSP] = sp->globalPosition().z();
            m_SPCL1_index[m_nSP] = cl->getIndex();
            m_SPCL2_index[m_nSP] = -1;
          }

          m_nSP++;
      	  if (m_nSP==m_maxCL) {
            std::cout<<"====================================="<<std::endl;
	          std::cout<<"DUMP : hit max number of space points"<<std::endl;
            std::cout<<"====================================="<<std::endl;
      	    break;
          }
        }
      }
    }
  }

  if (SCT_SpacePointContainer && SCT_SpacePointContainer->size()>0) {
    for( const auto& spCollection : *SCT_SpacePointContainer){
      // skip empty collections
      if( spCollection->empty() ) continue;
      
      if (not doSaveFile) {
        if (m_csvFile)
          m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
        doSaveFile = true;
      }

      // loop over collection
      for( auto& sp : *spCollection ) {
        // save sp x, y, z and the index of the cluster associated to that one
        if (doSaveFile) {
          const InDet::SiCluster * cl_1 = static_cast<const InDet::SiCluster*>(sp->clusterList().first);
          const InDet::SiCluster * cl_2 = static_cast<const InDet::SiCluster*>(sp->clusterList().second);
          if (m_csvFile)
            *m_outfile<<sp_index++<<","<<sp->globalPosition().x()<<","<<sp->globalPosition().y()<<","<<sp->globalPosition().z()<<"," <<cl_1->getIndex()<<"," <<cl_2->getIndex()<<std::endl;
          if (m_rootFile) {
            m_SPindex[m_nSP] = sp_index-1;
            m_SPx[m_nSP] = sp->globalPosition().x();
            m_SPy[m_nSP] = sp->globalPosition().y();
            m_SPz[m_nSP] = sp->globalPosition().z();
            m_SPCL1_index[m_nSP] = cl_1->getIndex();
            m_SPCL2_index[m_nSP] = cl_2->getIndex();
          }

          m_nSP++;
      	  if (m_nSP==m_maxCL) {
            std::cout<<"====================================="<<std::endl;
	          std::cout<<"DUMP : hit max number of space points"<<std::endl;
            std::cout<<"====================================="<<std::endl;
      	    break;
          }
        }
      }
    }
  }

    /*
    /////overlap
    if (overlapCollection && overlapCollection->size()>0) {
        //for( const auto& spCollection : *overlapCollection){

            if (not doSaveFile && m_csvFile) {
                m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
                doSaveFile = true;
            }

            for( auto& sp : *overlapCollection ) {
                // save sp x, y, z and the index of the cluster associated to that one
                if (doSaveFile) {          
                    const InDet::SiCluster * cl_1 = static_cast<const InDet::SiCluster*>(sp->clusterList().first);
                    const InDet::SiCluster * cl_2 = static_cast<const InDet::SiCluster*>(sp->clusterList().second);
                    *m_outfile<<sp_index++<<","<<sp->globalPosition().x()<<","<<sp->globalPosition().y()<<","<<sp->globalPosition().z()<<"," <<cl_1->getIndex()<<"," <<cl_2->getIndex()<<std::endl;
                }
            }
    }*/

  if (doSaveFile && m_csvFile)
    m_outfile->close();

  const TrackCollection* trackCollection = 0;
  if (evtStore()->retrieve(trackCollection, "CombinedInDetTracks").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object CombinedInDetTracks !" );
    return StatusCode::FAILURE;
  }

  const TrackTruthCollection* trackTruthCollection = 0;
  trackTruthCollection = 0;
  if (evtStore()->retrieve(trackTruthCollection, "TrackTruthCollection").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object TrackTruthCollection !" );
    return StatusCode::FAILURE;
  }

  if (m_csvFile)
    filename = (m_name == "") ? "tracks_evt" : (m_name+"_tracks_evt");
  doSaveFile = false;
  int trk_index = 0;

  if (not doSaveFile) {
    if (m_csvFile)
      m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
    doSaveFile = true;
  }

  // loop over tracks (and track truth) objects
  TrackCollection::const_iterator trackIterator = (*trackCollection).begin();
  m_nTRK = 0;
  if (m_rootFile) {
    (*m_TRKproperties).clear();
    (*m_TRKpattern).clear();
    (*m_TRKperigee_position).clear();
    (*m_TRKperigee_momentum).clear();
    (*m_TRKmeasurementsOnTrack_pixcl_sctcl_index).clear();
    (*m_TRKoutliersOnTrack_pixcl_sctcl_index).clear();
  }

  for ( ; trackIterator < (*trackCollection).end(); ++trackIterator) {
    if (!((*trackIterator))) {
      std::cout <<"TrackCollection contains empty entries" << std::endl;
      continue;
    }
    const Trk::TrackInfo& info = (*trackIterator)->info();
    const Trk::FitQuality* fitQuality = (*trackIterator)->fitQuality();
    const Trk::Perigee* perigeeParameters = (*trackIterator)->perigeeParameters ();
    const DataVector<const Trk::MeasurementBase>* measurementsOnTrack = (*trackIterator)->measurementsOnTrack();
    const DataVector<const Trk::MeasurementBase>* outliersOnTrack = (*trackIterator)->outliersOnTrack();

    ElementLink<TrackCollection> tracklink;
    tracklink.setElement(const_cast<Trk::Track*>(*trackIterator));
    tracklink.setStorableObject(*trackCollection);
    const ElementLink<TrackCollection> tracklink2=tracklink;
    TrackTruthCollection::const_iterator found = trackTruthCollection->find(tracklink2);

    if (doSaveFile) {
      if (m_csvFile)
        *m_outfile << trk_index++ << "," << info.trackFitter() << "," << info.particleHypothesis() << ",#,";
      const std::bitset<Trk::TrackInfo::NumberOfTrackProperties>& properties = info.properties();
      std::vector<int> v_properties;
      for (std::size_t i = 0; i < properties.size(); i++) {
	      if (properties[i]){
          if (m_csvFile)
             *m_outfile << i << ",";
           v_properties.push_back(i);
        }
      }

      if (m_csvFile) *m_outfile << "#,#,";
      const std::bitset<Trk::TrackInfo::NumberOfTrackRecoInfo>& pattern = info.patternRecognition();
      std::vector<int> v_pattern;
      for (std::size_t i = 0; i < pattern.size(); i++) {
        if (pattern[i]) {
          if (m_csvFile)
            *m_outfile << i << ",";
          v_pattern.push_back(i);
        }
      }

      if (m_csvFile) *m_outfile << "#,";
      int ndof=-1;
      float chiSq=0;
      if (fitQuality) {
      	ndof=fitQuality->numberDoF();
	      chiSq=fitQuality->chiSquared();
      }
      if (m_csvFile) *m_outfile << ndof << "," << chiSq << ",";
      std::vector<double> position, momentum;
      int charge = 0;
      if (perigeeParameters) {
        if (m_csvFile) {
        	*m_outfile << perigeeParameters->charge() <<","
	  	               << perigeeParameters->position()[0] << "," << perigeeParameters->position()[1] << "," << perigeeParameters->position()[2] << ","
              		   << perigeeParameters->momentum()[0] << "," << perigeeParameters->momentum()[1] << "," << perigeeParameters->momentum()[2] << ",";
        }
        position.push_back(perigeeParameters->position()[0]);
        position.push_back(perigeeParameters->position()[1]);
        position.push_back(perigeeParameters->position()[2]);
        momentum.push_back(perigeeParameters->momentum()[0]);
        momentum.push_back(perigeeParameters->momentum()[1]);
        momentum.push_back(perigeeParameters->momentum()[2]);
        charge = perigeeParameters->charge();
      }
      else {
      	if(m_csvFile) *m_outfile << "0,0,0,0,0,0,0,";
        position.push_back(0);
        position.push_back(0);
        position.push_back(0);
        momentum.push_back(0);
        momentum.push_back(0);
        momentum.push_back(0);
      }
      int mot=0;
      int oot=0;
      if (measurementsOnTrack) mot=measurementsOnTrack->size();
      if (outliersOnTrack) oot=outliersOnTrack->size();
      if (m_csvFile) *m_outfile << mot <<"," << oot << ",#,";
      std::vector<int> measurementsOnTrack_pixcl_sctcl_index, outliersOnTrack_pixcl_sctcl_index;
      int TTCindex, TTCevent_index, TTCparticle_link;
      float TTCprobability;
      if (measurementsOnTrack)
	      for (size_t i=0; i<measurementsOnTrack->size(); i++) {
	        const Trk::MeasurementBase* mb = (*measurementsOnTrack)[i];
	        const InDet::PixelClusterOnTrack* pixcl = dynamic_cast<const InDet::PixelClusterOnTrack*>(mb);
	        const InDet::SCT_ClusterOnTrack* sctcl = dynamic_cast<const InDet::SCT_ClusterOnTrack*>(mb);
	        if (pixcl) {
            if (m_csvFile)
  	          *m_outfile << pixcl->prepRawData()->getIndex() << ",";
            measurementsOnTrack_pixcl_sctcl_index.push_back(pixcl->prepRawData()->getIndex());
          }

	        else
            if (sctcl) {
              if (m_csvFile)
          	    *m_outfile << sctcl->prepRawData()->getIndex() << ",";
              measurementsOnTrack_pixcl_sctcl_index.push_back(sctcl->prepRawData()->getIndex());
            }
            else {
              if (m_csvFile)
                *m_outfile << "-1,";
              measurementsOnTrack_pixcl_sctcl_index.push_back(-1);
            }
     	}

      if (m_csvFile) *m_outfile << "#,";
      if (outliersOnTrack)
      	for (size_t i=0; i<outliersOnTrack->size(); i++) {
      	  const Trk::MeasurementBase* mb = (*outliersOnTrack)[i];
	        const InDet::PixelClusterOnTrack* pixcl = dynamic_cast<const InDet::PixelClusterOnTrack*>(mb);
	        const InDet::SCT_ClusterOnTrack* sctcl = dynamic_cast<const InDet::SCT_ClusterOnTrack*>(mb);
	        if (pixcl) {
            if (m_csvFile)
        	    *m_outfile << pixcl->prepRawData()->getIndex() << ",";
            outliersOnTrack_pixcl_sctcl_index.push_back(pixcl->prepRawData()->getIndex());
          }
	        else
            if (sctcl) {
              if (m_csvFile)
          	    *m_outfile << sctcl->prepRawData()->getIndex() << ",";
              outliersOnTrack_pixcl_sctcl_index.push_back(sctcl->prepRawData()->getIndex());
            }
            else {
              if (m_csvFile)
          	    *m_outfile << "-1,";
              outliersOnTrack_pixcl_sctcl_index.push_back(-1);
            }
      	}
 
      if (found != trackTruthCollection->end()) {
        if (m_csvFile)
  	      *m_outfile << "#," << found->first.index() << "," << found->second.particleLink().eventIndex() << "," <<found->second.particleLink().barcode() << "," << found->second.probability() << std::endl;
        TTCindex = found->first.index();
        TTCevent_index = found->second.particleLink().eventIndex();
        TTCparticle_link = found->second.particleLink().barcode();
        TTCprobability = found->second.probability();
      }
      else {
        if (m_csvFile)
          *m_outfile << "#,-999,-999,-999,-1" << std::endl;
        TTCindex = TTCevent_index = TTCparticle_link = -999;
        TTCprobability = -1;
      }

      if (m_rootFile) {
        m_TRKindex[m_nTRK] = trk_index-1;
        m_TRKtrack_fitter[m_nTRK] = info.trackFitter();
        m_TRKndof[m_nTRK] = info.trackFitter();
        m_TRKparticle_hypothesis[m_nTRK] = info.particleHypothesis();
        (*m_TRKproperties).push_back(v_properties);
        (*m_TRKpattern).push_back(v_pattern);
        m_TRKndof[m_nTRK] = ndof;
        m_TRKchiSq[m_nTRK] = chiSq;
        (*m_TRKmeasurementsOnTrack_pixcl_sctcl_index).push_back(measurementsOnTrack_pixcl_sctcl_index);
        (*m_TRKoutliersOnTrack_pixcl_sctcl_index).push_back(outliersOnTrack_pixcl_sctcl_index);
        m_TRKcharge[m_nTRK] = charge;
        (*m_TRKperigee_position).push_back(position);
        (*m_TRKperigee_momentum).push_back(momentum);
        m_TRKmot[m_nTRK] = mot;
        m_TRKoot[m_nTRK] = oot;
        m_TTCindex[m_nTRK] = TTCindex;
        m_TTCevent_index[m_nTRK] = TTCevent_index;
        m_TTCparticle_link[m_nTRK] = TTCparticle_link;
        m_TTCprobability[m_nTRK] = TTCprobability;
      }

      // index
      m_nTRK++;
  	  if (m_nTRK==m_maxCL) {
        std::cout<<"====================================="<<std::endl;
        std::cout<<"DUMP : hit max number of track events"<<std::endl;
        std::cout<<"====================================="<<std::endl;
   	    break;
      }
    }
  }

  if (doSaveFile && m_csvFile) {
    m_outfile->close();
  }

  const DetailedTrackTruthCollection* detailedTrackTruthCollection = 0;
  detailedTrackTruthCollection = 0;
  if (evtStore()->retrieve(detailedTrackTruthCollection, "DetailedTrackTruth").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object DetailedTrackTruth !" );
    return StatusCode::FAILURE;
  }

  filename = (m_name == "") ? "detailedtracktruth_evt" : (m_name+"_detailedtracktruth_evt");
  doSaveFile = false;

  m_nDTT = 0;
  if (m_rootFile) {
    (*m_DTTtrajectory_eventindex).clear();
    (*m_DTTtrajectory_barcode).clear();
    (*m_DTTstTruth_subDetType).clear();
    (*m_DTTstTrack_subDetType).clear();
    (*m_DTTstCommon_subDetType).clear();
  }

  if (not doSaveFile) {
    if (m_csvFile)
      m_outfile = new std::ofstream(("./"+filename+std::to_string(m_event)+".txt").c_str());
    doSaveFile = true;
  }

  // loop over DetailedTrackTruth objects
  DetailedTrackTruthCollection::const_iterator detailedTrackTruthIterator = (*detailedTrackTruthCollection).begin();
  for ( ; detailedTrackTruthIterator != (*detailedTrackTruthCollection).end(); ++detailedTrackTruthIterator) {
    std::vector<int> DTTtrajectory_eventindex, DTTtrajectory_barcode, DTTstTruth_subDetType, DTTstTrack_subDetType, DTTstCommon_subDetType;
    if (m_csvFile)
      *m_outfile << detailedTrackTruthIterator->first.index() << ",";
    const TruthTrajectory& traj=detailedTrackTruthIterator->second.trajectory();
    if (m_csvFile) *m_outfile << traj.size() << ",#,";
    for (size_t j=0; j<traj.size(); j++) {
      if (m_csvFile)
        *m_outfile << traj[j].eventIndex() << "," <<traj[j].barcode()<<",";
      DTTtrajectory_eventindex.push_back(traj[j].eventIndex());
      DTTtrajectory_barcode.push_back(traj[j].barcode());
    }
    if (m_csvFile) *m_outfile << "#,";
    const SubDetHitStatistics& stTruth = detailedTrackTruthIterator->second.statsTruth();
    const SubDetHitStatistics& stTrack = detailedTrackTruthIterator->second.statsTrack();
    const SubDetHitStatistics& stCommon = detailedTrackTruthIterator->second.statsCommon();
    for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) {
      if (m_csvFile)
        *m_outfile << stTruth[SubDetHitStatistics::SubDetType(j)] << ",";
      DTTstTruth_subDetType.push_back(stTruth[SubDetHitStatistics::SubDetType(j)]);
    }
    for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) {
      if (m_csvFile)
        *m_outfile << stTrack[SubDetHitStatistics::SubDetType(j)] << ",";
      DTTstTrack_subDetType.push_back(stTrack[SubDetHitStatistics::SubDetType(j)]);
    }
    for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) {
      if (m_csvFile)
       *m_outfile << stCommon[SubDetHitStatistics::SubDetType(j)];
      DTTstCommon_subDetType.push_back(stCommon[SubDetHitStatistics::SubDetType(j)]);
      if (j<SubDetHitStatistics::NUM_SUBDETECTORS-1)
        if (m_csvFile) *m_outfile << ",";
    }
    if (m_csvFile) *m_outfile << std::endl;

    if (m_rootFile) {
      m_DTTindex[m_nDTT] = detailedTrackTruthIterator->first.index();
       m_DTTsize[m_nDTT] = traj.size();
      (*m_DTTtrajectory_eventindex).push_back(DTTtrajectory_eventindex);
      (*m_DTTtrajectory_barcode).push_back(DTTtrajectory_barcode);
      (*m_DTTstTruth_subDetType).push_back(DTTstTruth_subDetType);
      (*m_DTTstTrack_subDetType).push_back(DTTstTrack_subDetType);
      (*m_DTTstCommon_subDetType).push_back(DTTstCommon_subDetType);
    }

    m_nDTT++;
  }

  if (doSaveFile && m_csvFile)
    m_outfile->close();

  // Once all the information for this event has been filled in the arrays,
  // copy content of the arrays to the TTree
  if (m_rootFile)
    m_nt->Fill();

  return StatusCode::SUCCESS;

  // For the production version of the code, simply delete everything below
  //


  //const TrackCollection* trackCollection = 0;
  trackCollection = 0;
  if (evtStore()->retrieve(trackCollection, "CombinedInDetTracks").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object CombinedInDetTracks !" );
    return StatusCode::FAILURE;
  }

  //const TrackTruthCollection* trackTruthCollection = 0;
  trackTruthCollection = 0;
  if (evtStore()->retrieve(trackTruthCollection, "TrackTruthCollection").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object TrackTruthCollection !" );
    return StatusCode::FAILURE;
  }

  //const DetailedTrackTruthCollection* detailedTrackTruthCollection = 0;
  detailedTrackTruthCollection = 0;
  if (evtStore()->retrieve(detailedTrackTruthCollection, "DetailedTrackTruth").isFailure() ){
    ATH_MSG_ERROR( "Could not find the data object DetailedTrackTruth !" );
    return StatusCode::FAILURE;
  }

  std::cout<<"JAN : "<<(*trackCollection).size()<<" "<<(*trackTruthCollection).size()<<" "<<(*detailedTrackTruthCollection).size()<<std::endl;

  TrackTruthCollection::const_iterator trackTruthIterator = (*trackTruthCollection).begin();
  //trackTruthIterator = (*trackTruthCollection).begin();
  for ( ; trackTruthIterator != (*trackTruthCollection).end(); ++trackTruthIterator) {
    std::cout<<"JAN track truth : "<<trackTruthIterator->second.particleLink().eventIndex()<<" "<<trackTruthIterator->second.particleLink().barcode()<<" "<<trackTruthIterator->second.probability()<<std::endl;
    auto result = detailedTrackTruthCollection->equal_range(trackTruthIterator->first);
    //for (auto i = detailedTrackTruthCollection->find(trackTruthIterator->first); i != detailedTrackTruthCollection->end(); i++) {
    for (auto i=result.first; i != result.second; i++) {
      const TruthTrajectory& traj=i->second.trajectory();
      std::cout<<"JAN trajectory : "<<traj.size()<<std::endl;
      for (size_t j=0; j<traj.size(); j++) {
	std::cout<<"JAN "<<traj[j].eventIndex()<<" "<<traj[j].barcode()<<" "<<std::endl;
      }
      const SubDetHitStatistics& stTruth = i->second.statsTruth();
      const SubDetHitStatistics& stTrack = i->second.statsTrack();
      const SubDetHitStatistics& stCommon = i->second.statsCommon();
     std::cout<<"JAN : ";
      for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) std::cout<<stTruth[SubDetHitStatistics::SubDetType(j)]<<",";
      std::cout<<std::endl;
      std::cout<<"JAN : ";
      for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) std::cout<<stTrack[SubDetHitStatistics::SubDetType(j)]<<",";
      std::cout<<std::endl;
      std::cout<<"JAN : ";
      for (unsigned j=0; j<SubDetHitStatistics::NUM_SUBDETECTORS; j++) std::cout<<stCommon[SubDetHitStatistics::SubDetType(j)]<<",";
      std::cout<<std::endl;
    }
  }

  //DetailedTrackTruthCollection::const_iterator dttIterator = (*detailedTrackTruthCollection).begin();
  //for ( ; dttIterator < (*detailedTrackTruthCollection).end(); ++dttIterator) {
  //TruthTrajectory& traj=(*dttIterator)->trajectory();
  //std::cout<<"JAN : "<<traj.size()<<" ";
  //for (size_t i=0; i<traj.size(); i++) {
  //  std::cout<<traj[i].eventIndex()<<" "<<traj[i].barcode()<<" "<<std::endl;
  //}
  //const SubDetHitStatistics& stTruth = (*dttIterator)->statsTruth();
  //std::cout<<"JAN : ";
  //for (int i=0; i<SubDetHitStatistics::NUM_SUBDETECTORS; i++) std::cout<<stTruth[i]<<",";
  //std::cout<<std::endl;
  //}
  
  return StatusCode::SUCCESS;

  // loop over tracks
  //TrackCollection::const_iterator trackIterator = (*trackCollection).begin();
  trackIterator = (*trackCollection).begin();
  for ( ; trackIterator < (*trackCollection).end(); ++trackIterator) {
    if (!((*trackIterator))) {
      std::cout <<"TrackCollection contains empty entries" << std::endl;
      continue;
    }
    const Trk::TrackInfo& info = (*trackIterator)->info();
    const Trk::FitQuality* fitQuality = (*trackIterator)->fitQuality();
    const Trk::Perigee* perigeeParameters = (*trackIterator)->perigeeParameters ();
    const DataVector<const Trk::MeasurementBase>* measurementsOnTrack = (*trackIterator)->measurementsOnTrack();
    const DataVector<const Trk::MeasurementBase>* outliersOnTrack = (*trackIterator)->outliersOnTrack();
    const Trk::TrackSummary* trackSummary = (*trackIterator)->trackSummary();
    std::cout<<"Track: "<<info.dumpInfo()<<std::endl;
    if (fitQuality) {
      std::cout<<" "<<fitQuality->chiSquared()<<" "<<fitQuality->numberDoF()<<std::endl;
    } else {
      std::cout<<" no fit quality"<<std::endl;
    }
    if (perigeeParameters) {
      std::cout<<" "<<perigeeParameters->charge()<<"   "
	       <<perigeeParameters->position()[0]<<" "<<perigeeParameters->position()[1]<<" "<<perigeeParameters->position()[2]<<"   "
	       <<perigeeParameters->momentum()[0]<<" "<<perigeeParameters->momentum()[1]<<" "<<perigeeParameters->momentum()[2]<<std::endl;
    } else {
      std::cout<<" no perigee parameters"<<std::endl;
    }
    if (measurementsOnTrack) {
      std::cout<<" "<<measurementsOnTrack->size()<<std::endl;
      for (size_t i=0; i<measurementsOnTrack->size(); i++) {
	const Trk::MeasurementBase* mb = (*measurementsOnTrack)[i];
	const Trk::SpacePoint* sp = dynamic_cast<const Trk::SpacePoint*>(mb);
	const Trk::RIO_OnTrack* rio = dynamic_cast<const Trk::RIO_OnTrack*>(mb);
	const InDet::SiClusterOnTrack* sicl = dynamic_cast<const InDet::SiClusterOnTrack*>(mb);
	const InDet::PixelClusterOnTrack* pixcl = dynamic_cast<const InDet::PixelClusterOnTrack*>(mb);
	const InDet::SCT_ClusterOnTrack* sctcl = dynamic_cast<const InDet::SCT_ClusterOnTrack*>(mb);
	std::cout<<mb<<" "<<sp<<" "<<rio<<" "<<sicl<<" "<<pixcl<<" "<<sctcl<<std::endl;
	if (pixcl) {
	  std::cout<<pixcl->prepRawData()->getIndex()<<std::endl;
	} else if (sctcl) {
	  std::cout<<sctcl->prepRawData()->getIndex()<<std::endl;
	}
      }
    } else {
      std::cout<<" no measurements on track"<<std::endl;
    }
    if (outliersOnTrack) {
      std::cout<<" "<<outliersOnTrack->size()<<std::endl;
    } else {
      std::cout<<" no outliers on track"<<std::endl;
    }
    if (trackSummary) {
      std::cout<<" do have track summary"<<std::endl;
    } else {
      std::cout<<" no track summary"<<std::endl;
    }

    //virtual Amg::Vector3D position() const override final;


  }

  return StatusCode::SUCCESS;
  
}


//--------------------------------
StatusCode DumpObjects::finalize() {
//--------------------------------
  if (m_rootFile) {
    delete[] m_CLID;

    delete[] m_CLindex;
    delete m_CLhardware;
    delete[] m_CLx;
    delete[] m_CLy;
    delete[] m_CLz;
    delete[] m_CLbarrel_endcap;
    delete[] m_CLlayer_disk;
    delete[] m_CLeta_module;
    delete[] m_CLphi_module;
    delete[] m_CLside;
    delete[] m_CLmoduleID;
    delete m_CLparticleLink_eventIndex;
    delete m_CLparticleLink_barcode;
    delete m_CLbarcodesLinked;
    delete m_CLphis;
    delete m_CLetas;
    delete m_CLtots;
    delete[] m_CLloc_direction1;
    delete[] m_CLloc_direction2;
    delete[] m_CLloc_direction3;
    delete[] m_CLJan_loc_direction1;
    delete[] m_CLJan_loc_direction2;
    delete[] m_CLJan_loc_direction3;
    delete[] m_CLpixel_count;
    delete[] m_CLcharge_count;
    delete[] m_CLloc_eta;
    delete[] m_CLloc_phi;
    delete[] m_CLglob_eta;
    delete[] m_CLglob_phi;
    delete[] m_CLeta_angle;
    delete[] m_CLphi_angle;
    delete[] m_CLnorm_x;
    delete[] m_CLnorm_y;
    delete[] m_CLnorm_z;
    delete m_CLlocal_cov;

    delete[] m_CLevent_number;
    delete[] m_CLbarcode;
    delete[] m_CLpx;
    delete[] m_CLpy;
    delete[] m_CLpz;
    delete[] m_CLpt;
    delete[] m_CLeta;
    delete[] m_CLvx;
    delete[] m_CLvy;
    delete[] m_CLvz;
    delete[] m_CLradius;
    delete[] m_CLstatus;
    delete[] m_CLcharge;
    delete[] m_CLpdg_id;
    delete[] m_CLpassed;
//  if ((*m_CLvParentID).size()) delete[] m_CLvParentID;
//  if ((*m_CLvParentBarcode).size()) delete[] m_CLvParentBarcode;

    delete[] m_SPindex;
    delete[] m_SPx;
    delete[] m_SPy;
    delete[] m_SPz;
    delete[] m_SPCL1_index;
    delete[] m_SPCL2_index;

    delete[] m_TRKindex;
    delete[] m_TRKtrack_fitter;
    delete[] m_TRKparticle_hypothesis;
    delete m_TRKproperties;
    delete m_TRKpattern;
    delete[] m_TRKndof;
    delete[] m_TRKmot;
    delete[] m_TRKoot;
    delete[] m_TRKchiSq;
    delete m_TRKmeasurementsOnTrack_pixcl_sctcl_index;
    delete m_TRKoutliersOnTrack_pixcl_sctcl_index;
    delete[] m_TRKcharge;
    delete m_TRKperigee_position;
    delete m_TRKperigee_momentum;
    delete[] m_TTCindex;
    delete[] m_TTCevent_index;
    delete[] m_TTCparticle_link;
    delete[] m_TTCprobability;

    delete[] m_DTTindex;
    delete[] m_DTTsize;
    delete m_DTTtrajectory_eventindex;
    delete m_DTTtrajectory_barcode;
    delete m_DTTstTruth_subDetType;
    delete m_DTTstTrack_subDetType;
    delete m_DTTstCommon_subDetType;
  }

  return StatusCode::SUCCESS;

}

//--------------------------------------------------------------------------------------------
bool DumpObjects::isPassed(const HepMC::GenParticle * particle,
			   float& px, float& py, float& pz, float& pt, float& eta, 
			   float& vx, float& vy, float& vz, float& radius, 
			   float& status, float& charge,
			   std::vector<int> &vParentID, std::vector<int> &vParentBarcode,
                           int &vProdNin, int &vProdNout, int &vProdStatus, int &vProdBarcode) {
//--------------------------------------------------------------------------------------------
  
  px = particle->momentum().px();
  py = particle->momentum().py();
  pz = particle->momentum().pz();
  
  pt = std::sqrt(px*px+py*py);
  eta = particle->momentum().eta();
  
  int pdgCode = particle->pdg_id();
  
  int absPdgCode = std::abs(pdgCode);
  // get the charge: ap->charge() is used later, DOES NOT WORK RIGHT NOW
  const HepPDT::ParticleData* ap = m_particleDataTable->particle(absPdgCode);
  charge = 1.;
  if (ap) charge = ap->charge();
  // since the PDT table only has abs(PID) values for the charge
  charge *= (pdgCode > 0.) ?  1. : -1.;
  
  status = particle->status();
      
  if (particle->production_vertex()) {
    vx = particle->production_vertex()->position().x();
    vy = particle->production_vertex()->position().y();
    vz = particle->production_vertex()->position().z();
    radius = particle->production_vertex()->position().perp();
  } else {
    vx = vy = vz = -1; 
    radius = 999;    
    if (status==1)
      std::cout << "no vertex for particle with status 1" << std::endl;
  }

  if (particle->production_vertex()) {
    vProdNin=particle->production_vertex()->particles_in_size();
    vProdNout=particle->production_vertex()->particles_out_size();
    vProdStatus=particle->production_vertex()->id();
    vProdBarcode=particle->production_vertex()->barcode();
    for ( HepMC::GenVertex::particles_in_const_iterator  p = particle->production_vertex()->particles_in_const_begin(); p != particle->production_vertex()->particles_in_const_end(); ++p ) {
      vParentID.push_back((*p)->pdg_id());
      vParentBarcode.push_back((*p)->barcode());
    }
  } else {
    vProdNin=0;
    vProdNout=0;
    vProdStatus=-999;
    vProdBarcode=999;
  }
  
  bool passEta = (pt > 0.1) ? (std::abs(eta) < m_max_eta) : false;
  if (not passEta) return false;  

  bool passPt = (pt > m_min_pt);
  if (not passPt) return false;

  bool passBarcode = (particle->barcode() < m_max_barcode);
  if (not passBarcode) return false;

  bool passCharge = not (charge == 0.);
  if (not passCharge) return false;

  bool passStatus = (status==1);
  if (not passStatus) return false;

  bool passProdRadius = (radius < m_maxProdVertex);
  if (not passProdRadius) return false;
  
  return true;
}
