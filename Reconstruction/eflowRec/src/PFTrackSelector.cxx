#include "eflowRec/eflowTrackExtrapolatorBaseAlgTool.h"
#include "eflowRec/PFTrackSelector.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "xAODEgamma/ElectronxAODHelpers.h"


PFTrackSelector::PFTrackSelector(const std::string& name, ISvcLocator* pSvcLocator):
  AthAlgorithm(name, pSvcLocator),
  m_tracksReadHandle("InDetTrackParticles"),
  m_electronsReadHandle("eflowRec_selectedElectrons"),
  m_muonsReadHandle("eflowRec_selectedMuons"),
  m_eflowRecTracksWriteHandle("eflowRecTracks"),
  m_theTrackExtrapolatorTool("Trk::ParticleCaloExtensionTool",this),
  m_upperTrackPtCut(100.0)
{
  declareProperty("tracksName", m_tracksReadHandle);
  declareProperty("electronsName", m_electronsReadHandle);  
  declareProperty("muonsName",  m_muonsReadHandle);
  declareProperty("eflowRecTracksOutputName",  m_eflowRecTracksWriteHandle);
  declareProperty("trackExtrapolatorTool", m_theTrackExtrapolatorTool, "AlgTool to use for track extrapolation");
  declareProperty("trackSelectionTool", m_trackSelectorTool);
  declareProperty("upperTrackPtCut",m_upperTrackPtCut);
}

StatusCode PFTrackSelector::initialize(){

  ATH_CHECK(m_theTrackExtrapolatorTool.retrieve());  
  ATH_CHECK(m_trackSelectorTool.retrieve());

  ATH_CHECK(m_tracksReadHandle.initialize());
  ATH_CHECK(m_electronsReadHandle.initialize());
  ATH_CHECK(m_muonsReadHandle.initialize());

  ATH_CHECK(m_eflowRecTracksWriteHandle.initialize());
  
  return StatusCode::SUCCESS;

}

StatusCode PFTrackSelector::execute(){

  ATH_CHECK(m_eflowRecTracksWriteHandle.record(std::make_unique<eflowRecTrackContainer>()));

  /* Verify the read handle has a valid pointer, and if not return */
  if (!m_tracksReadHandle.isValid()){
    if (msgLvl(MSG::WARNING)) { msg(MSG::WARNING) << "Can not retrieve xAOD::TrackParticleContainer with name: " << m_tracksReadHandle.key() << endmsg; }
    return StatusCode::FAILURE;
  }

  /* Do the track selection for tracks to be used in all of the following steps: */
  xAOD::TrackParticleContainer::const_iterator itTrackParticle = m_tracksReadHandle->begin();
  int trackIndex = 0;
  for (; itTrackParticle != m_tracksReadHandle->end(); ++itTrackParticle, ++trackIndex) {
    const xAOD::TrackParticle* track = (*itTrackParticle);
    if (!track){
      if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "Have invalid pointer to xAOD::TrackParticle " << endmsg;
      continue;	
    }

    bool rejectTrack(!selectTrack(*track));

    bool isElectron = this->isElectron(track);
    bool isMuon = this->isMuon(track);
    if (true == isElectron || true == isMuon) rejectTrack = true;

    if (!rejectTrack) {
      /* Create the eflowRecCluster and put it in the container */
      std::unique_ptr<eflowRecTrack> thisEFRecTrack  = std::make_unique<eflowRecTrack>(ElementLink<xAOD::TrackParticleContainer>(*m_tracksReadHandle, trackIndex), m_theTrackExtrapolatorTool);
      thisEFRecTrack->setTrackId(trackIndex);
      m_eflowRecTracksWriteHandle->push_back(std::move(thisEFRecTrack));
    }
  }

  std::sort(m_eflowRecTracksWriteHandle->begin(), m_eflowRecTracksWriteHandle->end(), eflowRecTrack::SortDescendingPt());

  return StatusCode::SUCCESS;
}

StatusCode PFTrackSelector::finalize(){return StatusCode::SUCCESS;}

bool PFTrackSelector::selectTrack(const xAOD::TrackParticle& track) {
  if (track.pt()*0.001 < m_upperTrackPtCut) return m_trackSelectorTool->accept(track, track.vertex());
  else return false;
}

bool PFTrackSelector::isElectron(const xAOD::TrackParticle* track){

  if (m_electronsReadHandle.isValid()){

    xAOD::ElectronContainer::const_iterator firstElectron = m_electronsReadHandle->begin();
    xAOD::ElectronContainer::const_iterator lastElectron = m_electronsReadHandle->end();
    
    for (; firstElectron != lastElectron; ++firstElectron){
      const xAOD::Electron* this_egamma = *firstElectron;
      if (this_egamma){
	unsigned int nTrack = this_egamma->nTrackParticles();
	
	if (0 != nTrack){	  
	  const xAOD::TrackParticle* origTrack = xAOD::EgammaHelpers::getOriginalTrackParticle(this_egamma);	  
	  if (origTrack){
	    if (track == origTrack) {
	      return true;
	    }
	  }//if valid track pointer
	  else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "Electron object map has NULL pointer to original TrackParticle " << endmsg;
	}//if has a track
	else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "Electron object has " << nTrack << " tracks " << endmsg;
      }//if valid pointer
      else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "Electron is a NULL pointer " << endmsg;
    }//electron loop    
  }
  else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << " Invalid ReadHandle for electrons with key: " << m_electronsReadHandle.key() << endmsg;

  return false;

}

bool PFTrackSelector::isMuon(const xAOD::TrackParticle* track){

  if (m_muonsReadHandle.isValid()){

    xAOD::MuonContainer::const_iterator firstMuon = m_muonsReadHandle->begin();
    xAOD::MuonContainer::const_iterator lastMuon = m_muonsReadHandle->end();
    
    for (; firstMuon != lastMuon ; ++firstMuon){
      const xAOD::Muon* theMuon = *firstMuon;
      if (theMuon){
	const ElementLink< xAOD::TrackParticleContainer > theLink = theMuon->inDetTrackParticleLink();
	if (theLink.isValid()){
	  const xAOD::TrackParticle* ID_track = *theLink;
	  if (ID_track){
	    if (track == ID_track) return true;
	    return false;
	  }
	  else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "This muon has a NULL pointer to the track " << endmsg;
	}
	else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "This muon has an invalid link to the track " << endmsg;
      }//if muon pointer is valid
      else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << "This muon is a NULL pointer " << endmsg;
    }//muon loop
  }
   else if (msgLvl(MSG::WARNING)) msg(MSG::WARNING) << " Invalid ReadHandle for muons with key: " << m_muonsReadHandle.key() << endmsg;

  return false;
}
