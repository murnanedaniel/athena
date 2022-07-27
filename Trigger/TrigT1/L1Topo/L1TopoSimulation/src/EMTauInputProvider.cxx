/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "./EMTauInputProvider.h"

#include <math.h>

#include "GaudiKernel/ITHistSvc.h"

#include "TrigT1CaloEvent/EmTauROI_ClassDEF.h"


#include "TrigT1Interfaces/CPRoIDecoder.h"

#include "L1TopoEvent/ClusterTOB.h"
#include "L1TopoEvent/TopoInputEvent.h"


using namespace std;
using namespace LVL1;

EMTauInputProvider::EMTauInputProvider(const std::string& type, const std::string& name, 
                                       const IInterface* parent) :
   base_class(type, name, parent),
   m_histSvc("THistSvc", name),
   m_emTauLocation( TrigT1CaloDefs::EmTauTopoTobLocation )
{
   declareInterface<LVL1::IInputTOBConverter>( this );
   declareProperty( "EmTauROILocation", m_emTauLocation, "Storegate key for the EMTAU info from CMX" );
}

EMTauInputProvider::~EMTauInputProvider()
{}

StatusCode
EMTauInputProvider::initialize() {

   CHECK(m_histSvc.retrieve());

   ServiceHandle<IIncidentSvc> incidentSvc("IncidentSvc", "EnergyInputProvider");
   CHECK(incidentSvc.retrieve());
   incidentSvc->addListener(this,"BeginRun", 100);
   incidentSvc.release().ignore();

   CHECK(m_emTauLocation.initialize());

   return StatusCode::SUCCESS;
}


void
EMTauInputProvider::handle(const Incident& incident) {
   if (incident.type()!="BeginRun") return;
   ATH_MSG_DEBUG( "In BeginRun incident");

   string histPath = "/EXPERT/" + name() + "/";
   replace( histPath.begin(), histPath.end(), '.', '/'); 

   auto hEMEt = std::make_unique<TH1I>( "EMTOBEt", "EM TOB Et", 200, 0, 400);
   hEMEt->SetXTitle("E_{T} [GeV]");
   auto hEMPhiEta = std::make_unique<TH2I>( "EMTOBPhiEta", "EM TOB Location", 200, -50, 50, 128, 0, 64);
   hEMPhiEta->SetXTitle("#eta#times10");
   hEMPhiEta->SetYTitle("#phi#times10");
   auto hEMEtEta = std::make_unique<TH2I>( "EMTOBEtEta", "Et vs eta", 200, -50, 50, 200, 0, 400);
   hEMEtEta->SetXTitle("#eta#times10");
   hEMEtEta->SetYTitle("E_{t} [GeV]");
   auto hEMEtPhi = std::make_unique<TH2I>( "EMTOBEtPhi", "Et vs phi", 128, 0, 64, 200, 0, 400);
   hEMEtPhi->SetXTitle("#phi#times10");
   hEMEtPhi->SetYTitle("E_{t} [GeV]");

   auto hTauEt = std::make_unique<TH1I>( "TauTOBEt", "Tau TOB Et", 200, 0, 400);
   hTauEt->SetXTitle("E_{T} [GeV]");
   auto hTauPhiEta = std::make_unique<TH2I>( "TauTOBPhiEta", "Tau TOB Location", 200, -50, 50, 128, 0, 64);
   hTauPhiEta->SetXTitle("#eta#times10");
   hTauPhiEta->SetYTitle("#phi#times10");
   auto hTauEtEta = std::make_unique<TH2I>( "TauTOBEtEta", "Et vs eta", 200, -50, 50, 200, 0, 400);
   hTauEtEta->SetXTitle("#eta#times10");
   hTauEtEta->SetYTitle("E_{t} [GeV]");
   auto hTauEtPhi = std::make_unique<TH2I>( "TauTOBEtPhi", "Et vs phi", 128, 0, 64, 200, 0, 400);
   hTauEtPhi->SetXTitle("#phi#times10");
   hTauEtPhi->SetYTitle("E_{t} [GeV]");


   if (m_histSvc->regShared( histPath + "EMTOBEt", std::move(hEMEt), m_hEMEt ).isSuccess()){
     ATH_MSG_DEBUG("EMTOBEt histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register EMTOBEt histogram for EMTauProvider");
   }
   if (m_histSvc->regShared( histPath + "EMTOBPhiEta", std::move(hEMPhiEta), m_hEMPhiEta ).isSuccess()){
     ATH_MSG_DEBUG("EMTOBPhiEta histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register EMTOBPhiEta histogram for EMTauProvider");
   }
   if (m_histSvc->regShared( histPath + "EMTOBEtEta", std::move(hEMEtEta), m_hEMEtEta ).isSuccess()){
     ATH_MSG_DEBUG("EMTOBEtEta histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register EMTOBEtEta histogram for EMTauProvider");
   }
   if (m_histSvc->regShared( histPath + "EMTOBEtPhi", std::move(hEMEtPhi), m_hEMEtPhi ).isSuccess()){
     ATH_MSG_DEBUG("EMTOBEtPhi histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register EMTOBEtPhi histogram for EMTauProvider");
   }

   if (m_histSvc->regShared( histPath + "TauTOBEt", std::move(hTauEt), m_hTauEt ).isSuccess()){
     ATH_MSG_DEBUG("TauTOBEt histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register TauTOBEt histogram for EMTauProvider");
   }
   if (m_histSvc->regShared( histPath + "TauTOBPhiEta", std::move(hTauPhiEta), m_hTauPhiEta ).isSuccess()){
     ATH_MSG_DEBUG("TauTOBPhiEta histogram has been registered successfully for EMTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register TauTOBPhiEta histogram for EMTauProvider");
   }
   if (m_histSvc->regShared( histPath + "TauTOBEtEta", std::move(hTauEtEta), m_hTauEtEta ).isSuccess()){
     ATH_MSG_DEBUG("TauTOBEtEta histogram has been registered successfully for TauTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register TauTOBEtEta histogram for TauTauProvider");
   }
   if (m_histSvc->regShared( histPath + "TauTOBEtPhi", std::move(hTauEtPhi), m_hTauEtPhi ).isSuccess()){
     ATH_MSG_DEBUG("TauTOBEtPhi histogram has been registered successfully for TauTauProvider.");
   }
   else{
     ATH_MSG_WARNING("Could not register TauTOBEtPhi histogram for TauTauProvider");
   }
   
}



StatusCode
EMTauInputProvider::fillTopoInputEvent(TCS::TopoInputEvent& inputEvent) const {

   /** this will be the new format
       https://indico.cern.ch/conferenceDisplay.py?confId=284687
       Electron ROI:
       | 0 0 1 0 | 2b Crate | 4b CPM Num | 3b CPChip | 3b Local coords | 0 0 0 | 5b electron isolation/veto | 8b electron energy |
       Tau ROI:
       | 0 0 1 1 | 2b Crate | 4b CPM Num | 3b CPChip | 3b Local coords | 0 0 0 | 5b tau isolation/veto      | 8b tau energy      |
   */

   // Retrieve EMTAU RoIs (they are built by EMTAUTrigger)

   SG::ReadHandle<DataVector < CPCMXTopoData> > emtau(m_emTauLocation);
   if( !emtau.isValid() ) {
      ATH_MSG_WARNING("No CPCMXTopoDataCollection with SG key '" << m_emTauLocation.key() << "' found in the event. No EM or TAU input for the L1Topo simulation.");
      return StatusCode::RECOVERABLE;
   }

   ATH_MSG_DEBUG("Filling the input event. Number of emtau topo data objects: " << emtau->size());
   // topoData is read in in reverse in order to obtain the same crate order as in the hardware (3->0, not 0->3)
   for(auto iTopoData = emtau->rbegin(); iTopoData != emtau->rend(); ++iTopoData) {
      const CPCMXTopoData *topoData = *iTopoData;

      // fill the vector of TOBs
      std::vector< CPTopoTOB > tobs;
      topoData->tobs(tobs);
      ATH_MSG_DEBUG("Emtau topo data object has # TOBs: " << tobs.size());
      for(const CPTopoTOB & tob : tobs ) {
         ATH_MSG_DEBUG( "EMTAU TOB with cmx = " << tob.cmx() << "[" << (tob.cmx()==0?"EM":"TAU") << "]"
                        << " : e = " << setw(3) << tob.et() << ", isolation " << tob.isolation()
                        << ", eta = " << setw(2) << tob.eta() << ", phi = " << tob.phi()
                        << ", ieta = " << setw(2) << tob.ieta() << ", iphi = " << tob.iphi()
                        << ", word = " << hex << tob.roiWord() << dec
                        );

         TCS::ClusterTOB cl(tob.et(), tob.isolation(), tob.ieta(), tob.iphi(), tob.cmx()==0 ? TCS::CLUSTER : TCS::TAU, tob.roiWord() );
         cl.setEtaDouble( tob.eta() );
         cl.setPhiDouble( tob.phi() );
         
         if(tob.cmx()==0) {
            inputEvent.addCluster( cl );
            m_hEMEt->Fill(cl.Et());
            m_hEMPhiEta->Fill(cl.eta(),cl.phi());
            m_hEMEtEta->Fill(cl.eta(),cl.Et());
            m_hEMEtPhi->Fill(cl.phi(),cl.Et());
         } else {
            inputEvent.addTau( cl );            
            m_hTauEt->Fill(cl.Et());
            m_hTauPhiEta->Fill(cl.eta(),cl.phi());
            m_hTauEtEta->Fill(cl.eta(),cl.Et());
            m_hTauEtPhi->Fill(cl.phi(),cl.Et());
         }
      }
      if(topoData->overflow()){
          inputEvent.setOverflowFromEmtauInput(true);
          ATH_MSG_DEBUG("setOverflowFromEmtauInput : true");
      }
   }
   return StatusCode::SUCCESS;
}

void 
EMTauInputProvider::CalculateCoordinates(int32_t roiWord, double & eta, double & phi) const {
   CPRoIDecoder get;
   constexpr double TwoPI = 2 * M_PI;
   CoordinateRange coordRange = get.coordinate( roiWord );
   
   eta = coordRange.eta();
   phi = coordRange.phi();
   if( phi > M_PI ) phi -= TwoPI;
}

