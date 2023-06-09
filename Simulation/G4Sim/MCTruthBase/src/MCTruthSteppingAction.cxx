/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "MCTruthSteppingAction.h"
#include "MCTruth/TrackHelper.h"
#include "RecordingEnvelope.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include <map>
#include <iostream>

MCTruthSteppingAction::MCTruthSteppingAction(const std::string& type,
                                             const std::string& name,
                                             const IInterface* parent)
  : UserActionBase(type, name, parent), m_needToInitialize(true)
{
  declareProperty("VolumeCollectionMap", m_volumeCollectionMap,
                  "Map of volume name to output collection name");
}

void MCTruthSteppingAction::BeginOfEvent (const G4Event*)
{
  if(m_needToInitialize)
    {
      // Fills m_recordingEnvelopes from ToolHandleArray
      if(this->FillRecordingEnvelopeVector().isFailure())
        {
          //bail out. TODO should we throw a G4Exception here perhaps?
          return;
        }
      if (m_recordingEnvelopes.size() == 0)
        {
          ATH_MSG_WARNING("No recording envelopes found!");
        }
    }
  for (auto& recEnvelope : m_recordingEnvelopes)
    {
      recEnvelope.BeginOfEvent();
    }
  return;
}

StatusCode MCTruthSteppingAction::FillRecordingEnvelopeVector()
{
  ATH_MSG_DEBUG("Setting up " << m_volumeCollectionMap.size() <<
                " recording envelopes:");
  m_recordingEnvelopes.reserve( m_volumeCollectionMap.size() );
  for(const auto& volCollPair : m_volumeCollectionMap) {
    ATH_MSG_DEBUG("  " << volCollPair.first << ", " << volCollPair.second);

    // Construct the helper in place on the vector
    m_recordingEnvelopes.emplace_back(volCollPair.first, volCollPair.second);
    RecordingEnvelope& recEnvelope = m_recordingEnvelopes.back();

    // Make sure the RecEnvelope can initialize properly
    if(!recEnvelope.Initialize()) {
      //FIXME - should this be an error?
      ATH_MSG_WARNING("Envelope volume " << recEnvelope.GetVolumeName() <<
                      " not found in geometry!");
      ATH_MSG_WARNING("TrackRecordCollection " <<
                      recEnvelope.GetTrackRecordCollectionName() <<
                      " will NOT be recorded");
      // Throw away uninitialized RecordingEnvelope
      m_recordingEnvelopes.pop_back();
    }
  }
  m_needToInitialize=false;
  return StatusCode::SUCCESS;
}

void MCTruthSteppingAction::Step(const G4Step* aStep)
{
  if (m_recordingEnvelopes.size() == 0) return;
  TrackHelper trackHelper(aStep->GetTrack());
  if (trackHelper.IsSecondary()) return;

  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4StepPoint* postStep = aStep->GetPostStepPoint();

  G4VPhysicalVolume* preVol = preStep->GetPhysicalVolume();
  G4VPhysicalVolume* postVol = postStep->GetPhysicalVolume();

  if (preVol==postVol) return;

  G4TouchableHistory* preTHist = (G4TouchableHistory*) preStep->GetTouchable();
  G4TouchableHistory* postTHist = (G4TouchableHistory*) postStep->GetTouchable();
  const int preStepVolumeDepth = preTHist->GetHistoryDepth();
  const int postStepVolumeDepth = postTHist->GetHistoryDepth();

  for (auto& recEnvelope : m_recordingEnvelopes)
    {
      const int envelopeLevel = recEnvelope.GetLevel();
      if (envelopeLevel <= preStepVolumeDepth)
        {
          //NB  preTHist->GetVolume(preStepVolumeDepth) would just give us the World volume.
          const G4LogicalVolume* logicalvolume1 =
            preTHist->GetVolume(preStepVolumeDepth-envelopeLevel)->GetLogicalVolume();
          if (logicalvolume1 != recEnvelope.GetLogicalVolume()) continue;

          if (envelopeLevel <= postStepVolumeDepth)
            {
              if (logicalvolume1 ==
                  postTHist->GetVolume(postStepVolumeDepth-envelopeLevel)->GetLogicalVolume())
                continue;
            }
          // We have a track crossing a recording envelope volume boundary, so make a TrackRecord
          recEnvelope.AddTrackRecord(aStep);
          // Done with this volume.
          break;
        }
    }
  return;
}

// TODO: drop in favor of declareInterface
StatusCode MCTruthSteppingAction::queryInterface(const InterfaceID& riid,
                                                 void** ppvInterface)
{
  if ( IUserAction::interfaceID().versionMatch(riid) ) {
    *ppvInterface = dynamic_cast<IUserAction*>(this);
    addRef();
  } else {
    // Interface is not directly available : try out a base class
    return UserActionBase::queryInterface(riid, ppvInterface);
  }
  return StatusCode::SUCCESS;
}


//=============================================================================
// Version2 of the MCTruthSteppingAction for multi-threading follows
//=============================================================================

namespace G4UA
{

  //---------------------------------------------------------------------------
  // Constructor
  //---------------------------------------------------------------------------
  MCTruthSteppingAction::
  MCTruthSteppingAction(const VolumeCollectionMap_t& volCollMap,
                        IMessageSvc* msgSvc, MSG::Level level)
    : AthMessaging(msgSvc, "MCTruthSteppingAction"),
      m_isInitialized(false),
      m_volumeCollectionMap(volCollMap)
  {
    msg().setLevel(level);
  }

  //---------------------------------------------------------------------------
  // Setup the recording envelopes
  //---------------------------------------------------------------------------
  void MCTruthSteppingAction::setupRecEnvelopes()
  {
    ATH_MSG_DEBUG("Setting up " << m_volumeCollectionMap.size() <<
                  " recording envelopes:");
    m_recordingEnvelopes.reserve( m_volumeCollectionMap.size() );
    for(const auto& volCollPair : m_volumeCollectionMap) {
      ATH_MSG_DEBUG("  " << volCollPair.first << ", " << volCollPair.second);

      // Construct the helper in place on the vector
      m_recordingEnvelopes.emplace_back(volCollPair.first, volCollPair.second);
      RecordingEnvelope& recEnvelope = m_recordingEnvelopes.back();

      // Make sure the RecEnvelope can initialize properly
      if(!recEnvelope.Initialize()) {
        //FIXME - should this be an error?
        ATH_MSG_WARNING("Envelope volume " << recEnvelope.GetVolumeName() <<
                        " not found in geometry!");
        ATH_MSG_WARNING("TrackRecordCollection " <<
                        recEnvelope.GetTrackRecordCollectionName() <<
                        " will NOT be recorded");
        // Throw away uninitialized RecordingEnvelope
        m_recordingEnvelopes.pop_back();
      }
    }
    if (m_recordingEnvelopes.size() == 0) {
      ATH_MSG_WARNING("No recording envelopes found!");
    }
  }

  //---------------------------------------------------------------------------
  // Beginning of event
  //---------------------------------------------------------------------------
  void MCTruthSteppingAction::beginOfEvent(const G4Event*)
  {
    // First time initialization
    if(!m_isInitialized) {
      setupRecEnvelopes();
      if (m_recordingEnvelopes.size() == 0) {
        ATH_MSG_WARNING("No recording envelopes found!");
      }
      m_isInitialized = true;
    }
    // Every event initialization
    for (auto& recEnvelope : m_recordingEnvelopes) {
      recEnvelope.BeginOfEvent();
    }
  }

  //---------------------------------------------------------------------------
  // Process one tracking step
  //---------------------------------------------------------------------------
  void MCTruthSteppingAction::processStep(const G4Step* aStep)
  {
    if (m_recordingEnvelopes.size() == 0) return;
    TrackHelper trackHelper(aStep->GetTrack());
    if (trackHelper.IsSecondary()) return;

    G4StepPoint* preStep = aStep->GetPreStepPoint();
    G4StepPoint* postStep = aStep->GetPostStepPoint();

    G4VPhysicalVolume* preVol = preStep->GetPhysicalVolume();
    G4VPhysicalVolume* postVol = postStep->GetPhysicalVolume();

    if (preVol == postVol) return;

    G4TouchableHistory* preTHist = (G4TouchableHistory*) preStep->GetTouchable();
    G4TouchableHistory* postTHist = (G4TouchableHistory*) postStep->GetTouchable();
    const int preStepVolDepth = preTHist->GetHistoryDepth();
    const int postStepVolDepth = postTHist->GetHistoryDepth();

    for (auto& recEnvelope : m_recordingEnvelopes)
    {
      const int envelopeLevel = recEnvelope.GetLevel();
      if (envelopeLevel <= preStepVolDepth)
      {
        //NB  preTHist->GetVolume(preStepVolDepth) would just give us the World volume.
        const G4LogicalVolume* logicalVolume1 =
          preTHist->GetVolume(preStepVolDepth-envelopeLevel)->GetLogicalVolume();
        if (logicalVolume1 != recEnvelope.GetLogicalVolume()) continue;

        if (envelopeLevel <= postStepVolDepth &&
            logicalVolume1 == postTHist->GetVolume(postStepVolDepth-envelopeLevel)
                              ->GetLogicalVolume())
        {
          continue;
        }

        // We have a track crossing a recording envelope
        // volume boundary, so make a TrackRecord
        recEnvelope.AddTrackRecord(aStep);

        // Done with this volume.
        break;
      }
    }
  }

} // namespace G4UA
