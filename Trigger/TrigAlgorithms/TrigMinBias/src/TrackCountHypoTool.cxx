
/*
Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigCompositeUtils/HLTIdentifier.h"
#include "AthenaMonitoringKernel/Monitored.h"

#include "TrackCountHypoTool.h"

using namespace TrigCompositeUtils;

TrackCountHypoTool::TrackCountHypoTool(const std::string &type, const std::string &name, const IInterface *parent)
	: AthCheckedComponent<AthAlgTool>(type, name, parent),
	  m_decisionId(HLT::Identifier::fromToolName(name)) {}


StatusCode TrackCountHypoTool::decide(TrkCountsInfo &trkinfo) const
{
	if (trkinfo.previousDecisionIDs.count(m_decisionId.numeric()) == 0)
	{
		ATH_MSG_DEBUG("Already rejected");
		return StatusCode::SUCCESS;
	}
	if ( m_acceptAll ) {		
		addDecisionID(m_decisionId.numeric(), trkinfo.decision);
		ATH_MSG_DEBUG("REGTEST event accepted because of acceptAll");
	}

	std::vector<int> counts;
	trkinfo.counts->getDetail<std::vector<int>>("counts", counts);
	if ( m_exclusive ) {
		if ( counts[0] > m_exclusivityThreshold ) {
			ATH_MSG_DEBUG("Lowest pt tracks count " << counts[0] << " exceeds exclusivity cut, " <<  m_exclusivityThreshold<<" rejecting");
			return StatusCode::SUCCESS;
		}
	}
	std::vector<float> pTcuts;
	std::vector<float> z0cuts;
	std::vector<float> vertexZcuts;
	trkinfo.counts->getDetail<std::vector<float>>("pTcuts", pTcuts);
	trkinfo.counts->getDetail<std::vector<float>>("z0cuts", z0cuts);
	trkinfo.counts->getDetail<std::vector<float>>("vertexZcuts", vertexZcuts);


	float countForConfiguredPtThreshold{};
	bool found{false};
	// find the right counter
	for (size_t i = 0; i < counts.size(); ++i)
	{
		if (std::abs(pTcuts[i] - m_minPt) < 0.001 
			&& std::abs(z0cuts[i] - m_maxZ0) < 0.001
			&& std::abs(vertexZcuts[i] - m_maxVertexZ) < 0.001)
		{
			found = true;
			countForConfiguredPtThreshold = counts[i];
		}
	}

	if (!found)
	{
		ATH_MSG_ERROR("Unable to find tracks count for requested: "  
					<< " min pT " << m_minPt 
					<< " nax z0  " << m_maxZ0
					<< " vertex Z " << m_maxVertexZ 
					<< " need to fix hypo tool configuration or add new threshold in tracks counting");
		for (size_t i = 0; i < counts.size(); ++i)
		{
			ATH_MSG_ERROR("Count of tracks of pTcuts " << pTcuts[i] << " z0Cuts "  << z0cuts[i] << " vertexZcuts " <<  vertexZcuts[i] << " that are available");
		}
		return StatusCode::FAILURE;
	}
	else
	{
		ATH_MSG_DEBUG("REGTEST found " << countForConfiguredPtThreshold << " tracks for " << m_minPt);
	}



	const bool minTrkPassed = (m_minNtrks == -1) or (countForConfiguredPtThreshold >= m_minNtrks);
	const bool maxTrkPassed = (m_maxNtrks == -1) or (countForConfiguredPtThreshold < m_maxNtrks);

	if ( minTrkPassed and maxTrkPassed ) {
		addDecisionID(m_decisionId.numeric(), trkinfo.decision);
		ATH_MSG_DEBUG("REGTEST event accepted");
	}
	return StatusCode::SUCCESS;
}
