/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#include "SimHelpers/SecondaryTracksHelper.h"



const G4TrackingManager* SecondaryTracksHelper::theTrackingManager=0;

SecondaryTracksHelper::SecondaryTracksHelper()
{
	currentNrOfSecondaries=0;
}
SecondaryTracksHelper::SecondaryTracksHelper(const G4TrackingManager *tm)
{
	theTrackingManager=tm;
	currentNrOfSecondaries=0;
}
void SecondaryTracksHelper::SetTrackingManager(const G4TrackingManager *tm)
{
	theTrackingManager=tm;
}
const G4TrackingManager *SecondaryTracksHelper::GetTrackingManager()
{
	return theTrackingManager;
}
int SecondaryTracksHelper::NrOfNewSecondaries()
{
	int value=0;
//	std::cout<<" NrOfNewSecondaries "<<currentNrOfSecondaries << 
//	         " "<<CurrentNrOfSecondaries()<<std::endl;
	int curr=CurrentNrOfSecondaries();
	if (curr>currentNrOfSecondaries) 
	{
		value=curr-currentNrOfSecondaries;
		currentNrOfSecondaries=curr;
	}
	return value;
}
std::vector<G4Track*> SecondaryTracksHelper::GetSecondaries(int ntracks) const
{
	std::vector<G4Track*> vTrack;
	int iSize=theTrackingManager->GimmeSecondaries()->size();
	for (int i=iSize-1;i>iSize-ntracks-1;i--) 
		vTrack.push_back((*(theTrackingManager->GimmeSecondaries()))[i]);
	return vTrack;
}
