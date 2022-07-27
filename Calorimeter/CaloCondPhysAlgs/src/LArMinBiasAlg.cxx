/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "LArMinBiasAlg.h"


#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/LArEM_ID.h"
#include "CaloIdentifier/CaloCell_ID.h"

#include "LArSimEvent/LArHit.h"
#include "LArSimEvent/LArHitContainer.h"
#include "TTree.h"
#include "GaudiKernel/ITHistSvc.h"
#include "CaloDetDescr/CaloDetDescrElement.h"



  //Constructor
  LArMinBiasAlg:: LArMinBiasAlg(const std::string& name, ISvcLocator* pSvcLocator):
    AthAlgorithm(name,pSvcLocator),
    m_datasetID_lowPt(119995),
    m_datasetID_highPt(119996),
    m_weight_lowPt(39.8606),
    m_weight_highPt(0.138128)
  {
     declareProperty("datasetID_lowPt",m_datasetID_lowPt);
     declareProperty("datasetID_highPt",m_datasetID_highPt);
     declareProperty("weight_highPt",m_weight_highPt);
     declareProperty("weight_lowPt",m_weight_lowPt);
     m_first=true;
  }
  
  //__________________________________________________________________________
  //Destructor
  LArMinBiasAlg::~LArMinBiasAlg()
  {
  }
  //__________________________________________________________________________
  StatusCode LArMinBiasAlg::initialize()
  {
    
    ATH_MSG_INFO ( " LArMinBiasAlg initialize() " );

    ATH_CHECK( service("THistSvc",m_thistSvc) );


    m_tree = new TTree("m_tree","Offset ntuple");
    m_tree->Branch("ncell",&m_nsymcell,"ncell/I");
    m_tree->Branch("nevt_total",&m_nevt_total,"nevt_total/I");
    m_tree->Branch("identifier",m_identifier,"identifier[ncell]/I");
    m_tree->Branch("region",m_region,"region[ncell]/I");
    m_tree->Branch("ieta",m_ieta,"ieta[ncell]/I");
    m_tree->Branch("layer",m_layer,"layer[ncell]/I");
    m_tree->Branch("region",m_region,"region[ncell]/I");
    m_tree->Branch("ieta",m_ieta,"ieta[ncell]/I");
    m_tree->Branch("eta",m_eta,"eta[ncell]/F");
    m_tree->Branch("phi",m_phi,"phi[ncell]/F");
    m_tree->Branch("nevt",m_nevt,"nevt[ncell]/D");
    m_tree->Branch("average",m_average,"average[ncell]/D");
    m_tree->Branch("rms",m_rms,"rms[ncell]/D");
    m_tree->Branch("reference",m_offset,"offset[ncell]/D");
   
    if( m_thistSvc->regTree("/file1/offset",m_tree).isFailure()) {
       ATH_MSG_WARNING(" cannot register ntuple " );
       return StatusCode::SUCCESS;
     }

    const CaloIdManager* mgr = nullptr;
    ATH_CHECK( detStore()->retrieve( mgr ) );
    m_larem_id   = mgr->getEM_ID();
    m_calo_id    = mgr->getCaloCell_ID();


    ATH_CHECK(m_mcSymKey.initialize());

    ATH_CHECK(m_cablingKey.initialize());

    ATH_CHECK(m_caloMgrKey.initialize());

    ATH_CHECK(m_eventInfoKey.initialize());

    m_n1=0;
    m_n2=0;

    m_nevt_total=0;

    return StatusCode::SUCCESS; 

  }
  //__________________________________________________________________________
  StatusCode LArMinBiasAlg::stop()
  {
    ATH_MSG_INFO("number of events in the two samples  " << m_n1 << " " << m_n2);
    this->fillNtuple();
    ATH_MSG_INFO(" stop after fill ntuple");
    return StatusCode::SUCCESS;
  }
  //__________________________________________________________________________
  StatusCode LArMinBiasAlg::finalize()
  {
    ATH_MSG_INFO(" finalize()");
    return StatusCode::SUCCESS; 
  }
  
  //__________________________________________________________________________
  StatusCode LArMinBiasAlg::execute()
  {
    //.............................................
    
    ATH_MSG_DEBUG(" LArMinBiasAlg execute()");

    const EventContext& ctx = Gaudi::Hive::currentContext();

    if (m_first) {

      SG::ReadCondHandle<CaloDetDescrManager> caloMgrHandle{m_caloMgrKey};
      ATH_CHECK(caloMgrHandle.isValid());

      SG::ReadCondHandle<LArMCSym>  mcsym      (m_mcSymKey, ctx);
      SG::ReadCondHandle<LArOnOffIdMapping> cablingHdl (m_cablingKey, ctx);
      const LArOnOffIdMapping* cabling{*cablingHdl};
      if(!cabling) {
         ATH_MSG_ERROR( "Do not have cabling mapping from key " << m_cablingKey.key() );
         return StatusCode::FAILURE;
      }

      m_ncell = m_calo_id->calo_cell_hash_max();

      ATH_MSG_INFO(" --- first event " << m_ncell);
      m_symCellIndex.resize(m_ncell,-1);
      std::vector<int> doneCell;
      doneCell.resize(m_ncell,-1);

      //m_readCell.resize(m_ncell,0);

      m_eCell.resize(m_ncell,0.);

      m_CellList.reserve(MAX_SYM_CELLS);
      int nsym=0; 
      // loop over cells
      // and find symmetry cells
      for (unsigned int i=0;i<((unsigned int)(m_ncell));i++) {
        IdentifierHash idHash=i;
        Identifier id=m_calo_id->cell_id(idHash);
        if (m_calo_id->is_tile(id)) continue;
        // convert cell id to symmetric identifier
        HWIdentifier hwid2 = mcsym->ZPhiSymOfl(id);
        Identifier id2 = cabling->cnvToIdentifier(hwid2);
        int i2 = (int) (m_calo_id->calo_cell_hash(id2));
        if(i2>=m_ncell) {
           ATH_MSG_WARNING("problem: i2: "<<i2<<" for id: "<<m_calo_id->print_to_string(id)<<" symmetrized: "<<m_calo_id->print_to_string(id2));
        }
        // we have already processed this hash => just need to associate cell i to the same symmetric cell
        if (doneCell[i2]>=0) {
           m_symCellIndex[i]=doneCell[i2];
        }
        // we have not already processed this hash, add an entry for this new symmetric cell
        else {  
           doneCell[i2]=nsym;
           m_symCellIndex[i] = nsym; 
           CellInfo cell;
           const CaloDetDescrElement* calodde = (*caloMgrHandle)->get_element(id);
           cell.eta =  calodde->eta();
           cell.phi = calodde->phi();
           cell.region = m_calo_id->region(id);
           cell.ieta = m_calo_id->eta(id);
           cell.layer = m_calo_id->calo_sample(id);
           cell.region = m_calo_id->region(id);
           cell.ieta = m_calo_id->eta(id);
           //cell.identifier = id2.get_identifier32().get_compact();
           cell.identifier = id2;
           cell.average=0.;
           cell.offset=0.;
           cell.rms=0.;
           cell.nevt=0.;
           m_CellList.push_back(cell);
           nsym++;
        }
      }
      ATH_MSG_INFO(" --- number of symmetric cells found " << nsym << " " << m_CellList.size());
      if (nsym>=MAX_SYM_CELLS) ATH_MSG_ERROR(" More than "<<MAX_SYM_CELLS<<" number of symmetric cells... Fix array size for ntuple writing !!!! ");
      m_nsymcell=nsym;
      m_first=false;
    }

    SG::ReadHandle<xAOD::EventInfo> eventInfo(m_eventInfoKey);
    if (!eventInfo.isValid()) {
      ATH_MSG_ERROR ("Could not retrieve EventInfo");
      return StatusCode::FAILURE;
    }
    int channelNumber = eventInfo->mcChannelNumber();

  m_nevt_total++;

  double weight=1.;

// Dataset ID for lowPt MinBias
  if (channelNumber==m_datasetID_lowPt) {
     weight = m_weight_lowPt;
     m_n1+=1;
  }
//  Dataset ID for highPt MinBias
  else if (channelNumber==m_datasetID_highPt) {
     weight = m_weight_highPt;
     m_n2+=1;
  }
  else {
     ATH_MSG_WARNING(" Neither low nor high Pt MinBias sample " << channelNumber << "  set weight to 1.0 ");
     weight=1.;
  } 


  if ((m_nevt_total%100)==1) ATH_MSG_INFO(" ---- process event number " << m_nevt_total <<  " " << channelNumber << "  weight " << weight);

   for (int i=0;i<m_ncell;i++) m_eCell[i]=0.;

    std::vector <std::string> HitContainer;
    HitContainer.emplace_back("LArHitEMB");
    HitContainer.emplace_back("LArHitEMEC");
    HitContainer.emplace_back("LArHitHEC");
    HitContainer.emplace_back("LArHitFCAL");
    for (unsigned int iHitContainer=0;iHitContainer<HitContainer.size();iHitContainer++)
    {
      const LArHitContainer* hit_container ;
      ATH_CHECK(evtStore()->retrieve(hit_container,HitContainer[iHitContainer]));
      int ihit = 0;
      for (const LArHit* hit : *hit_container)
      {       
          ihit++; 
          Identifier cellID=hit->cellID();
          double energy = hit->energy(); 
          double time =hit->time();
          int index = (int) (m_calo_id->calo_cell_hash(cellID));
          if (index < m_ncell && index>=0 && fabs(time)<25.) {
             m_eCell[index] += energy;
        }

      }  // loop over hits in container
    }  // loop over containers

    for (int i=0;i<m_ncell;i++) {
          addCell(i,m_eCell[i],0.,weight);
    }


  return StatusCode::SUCCESS;
 }

 void LArMinBiasAlg::addCell(int index, double  energy, double eshift, double weight) 
 {
          if (index < m_ncell && index>=0) {
            int index2= m_symCellIndex[index];
            if (index2<0) return;
            if (index2 >= ((int)(m_CellList.size())) ) {
	      ATH_MSG_INFO(" LArMinBiasAlg::addCell: for " << index << ", " << index2 << " is out-of bounds for list of size " << m_CellList.size());
	      return;
            }
            double oldN =  m_CellList[index2].nevt;
            double oldAverage = m_CellList[index2].average;
            double oldRMS     = m_CellList[index2].rms;
            double oldAverage2 = m_CellList[index2].offset;

            double frac = oldN/(weight+oldN);
            double Anew = weight+oldN;
            double newAverage = frac*oldAverage + weight*energy/Anew;
            double deltaE = energy-newAverage;
            double newRMS     = frac*(oldRMS + (newAverage-oldAverage)*(newAverage-oldAverage)) + weight*deltaE*deltaE/Anew;

            double newAverage2 = frac*oldAverage2 + weight*eshift/Anew;

            m_CellList[index2].nevt = Anew;
            m_CellList[index2].average = newAverage;
            m_CellList[index2].rms = newRMS; 
            m_CellList[index2].offset = newAverage2; 

          }
 }

 void LArMinBiasAlg::fillNtuple()
 {

   ATH_MSG_INFO(" in fillNtuple " << m_nsymcell);
   for (int i=0;i<m_nsymcell;i++) {
     m_identifier[i] = m_CellList[i].identifier.get_identifier32().get_compact();
     m_layer[i] = m_CellList[i].layer;
     m_region[i] = m_CellList[i].region;
     m_ieta[i] = m_CellList[i].ieta;
     m_eta[i] = m_CellList[i].eta;
     m_phi[i] = m_CellList[i].phi;
     m_nevt[i] = m_CellList[i].nevt;
     m_offset[i] = (float) (m_CellList[i].offset);
     m_average[i] = (float) (m_CellList[i].average);
     m_rms[i] = (float) (sqrt(m_CellList[i].rms));
   }
   m_tree->Fill();
   ATH_MSG_INFO(" after tree fill ");

   for (int i=0;i<m_nsymcell;i++) {
     m_CellList[i].nevt=0;
     m_CellList[i].offset=0.;
     m_CellList[i].average=0;
     m_CellList[i].rms=0;
   }
   ATH_MSG_INFO(" end of fillNtuple ");
 
 }
