/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

// ****************************************************************************
// ----------------------------------------------------------------------------
// FourMuonTool header file
//
// James Catmore <James.Catmore@cern.ch>

// ----------------------------------------------------------------------------
// ****************************************************************************
#ifndef BPHY4TOOL_H
#define BPHY4TOOL_H
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "InDetConversionFinderTools/InDetConversionFinderTools.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"

#include <vector>
#include <string>
/////////////////////////////////////////////////////////////////////////////

namespace Trk {
    class IVertexFitter;
    //    class VxCandidate;
    //    class TrackParticleBase;
    //    class VxTrackAtVertex;
    //    class RecVertex;
    class TrkV0VertexFitter;
    class ITrackSelectorTool;
    class V0Tools;
    //    class ExtendedVxCandidate;
}

namespace InDet { class VertexPointEstimator; }

namespace DerivationFramework {
    
    static const InterfaceID IID_FourMuonTool("FourMuonTool", 1, 0);
    
    // Struct and enum to associate muon pairs with track pairs
    // and make the program flow more straightforward
    struct Combination
    {
        std::vector<const xAOD::Muon*> muons;
        std::vector<unsigned int> quadIndices;
        std::pair<unsigned int, unsigned int> pairIndices;
        
        std::string combinationCharges() {
            std::string chargeStr = "";
            if (muons.at(0)->charge() > 0) {chargeStr += "+";}
            else {chargeStr += "-";}
            if (muons.at(1)->charge() > 0) {chargeStr += "+";}
            else {chargeStr += "-";}
            if (muons.size()==4) {
                if (muons.at(2)->charge() > 0) {chargeStr += "+";}
                else {chargeStr += "-";}
                if (muons.at(3)->charge() > 0) {chargeStr += "+";}
                else {chargeStr += "-";}
            }
            return chargeStr;
        }
        
        std::string combinationIndices() {
            std::string indexStr = "";
            std::stringstream ss;
            if (muons.size()==2) {
                ss.str(""); ss.clear();
                ss << pairIndices.first;
                indexStr+=ss.str();
                ss.str(""); ss.clear();
                ss << pairIndices.second;
                indexStr+=ss.str();
            }
            if (muons.size()==4) {
                for (unsigned int i=0; i<4; ++i) {
                    ss.str(""); ss.clear();
                    ss << quadIndices[i];
                    indexStr+=ss.str();
                }
            }
            return indexStr;
        }

        const xAOD::TrackParticle* GetMuonTrack(const xAOD::Muon* mu) const{
            auto& link = mu->inDetTrackParticleLink();
            return link.isValid() ? *link : nullptr;
        }

        std::vector<const xAOD::TrackParticle*> trackParticles(const std::string& specify) {
            std::vector<const xAOD::TrackParticle*> theTracks;
            bool oppCh(false);
            if (muons.at(0)->charge()*muons.at(1)->charge() < 0) oppCh=true;
            if (specify=="pair1") {
                theTracks.push_back(GetMuonTrack(muons.at(0)));
                theTracks.push_back(GetMuonTrack(muons.at(1)));
            }
            if (specify=="pair2") {
                theTracks.push_back(GetMuonTrack(muons.at(2)));
                theTracks.push_back(GetMuonTrack(muons.at(3)));
            }
            if (specify=="DC") {
                if (oppCh) {
                    theTracks.push_back(GetMuonTrack(muons.at(0)));
                    theTracks.push_back(GetMuonTrack(muons.at(1)));
                    theTracks.push_back(GetMuonTrack(muons.at(2)));
                    theTracks.push_back(GetMuonTrack(muons.at(3)));
                } else {
                    theTracks.push_back(GetMuonTrack(muons.at(0)));
                    theTracks.push_back(GetMuonTrack(muons.at(2)));
                    theTracks.push_back(GetMuonTrack(muons.at(1)));
                    theTracks.push_back(GetMuonTrack(muons.at(3)));
                }
            }
            if (specify=="AC") {
                theTracks.push_back(GetMuonTrack(muons.at(0)));
                theTracks.push_back(GetMuonTrack(muons.at(3)));
                theTracks.push_back(GetMuonTrack(muons.at(1)));
                theTracks.push_back(GetMuonTrack(muons.at(2)));
            }
            if (specify=="SS") {
                if (oppCh) {
                    theTracks.push_back(GetMuonTrack(muons.at(0)));
                    theTracks.push_back(GetMuonTrack(muons.at(2)));
                    theTracks.push_back(GetMuonTrack(muons.at(1)));
                    theTracks.push_back(GetMuonTrack(muons.at(3)));
                } else {
                    theTracks.push_back(GetMuonTrack(muons.at(0)));
                    theTracks.push_back(GetMuonTrack(muons.at(1)));
                    theTracks.push_back(GetMuonTrack(muons.at(2)));
                    theTracks.push_back(GetMuonTrack(muons.at(3)));
                }
            }
            return theTracks;
        }
        
    };
    
    class FourMuonTool:  virtual public AthAlgTool
    {
    public:
        FourMuonTool(const std::string& t, const std::string& n, const IInterface*  p);
        ~FourMuonTool();
        StatusCode initialize();
        
        static const InterfaceID& interfaceID() { return IID_FourMuonTool;}
        
        //-------------------------------------------------------------------------------------
        //Doing Calculation and inline functions
        StatusCode performSearch(xAOD::VertexContainer*& pairVxContainer, xAOD::VertexAuxContainer*& pairVxAuxContainer,
                                 xAOD::VertexContainer*& quadVxContainer, xAOD::VertexAuxContainer*& quadVxAuxContainer, bool &acceptEvent) const;
        xAOD::Vertex* fit(const std::vector<const xAOD::TrackParticle*>& ,const xAOD::TrackParticleContainer* importedTrackCollection, const Amg::Vector3D &beamSpot) const;
        static std::vector<std::vector<unsigned int> > getQuadIndices(unsigned int length);
        static std::vector<std::pair<unsigned int, unsigned int> > getPairIndices(unsigned int length);
        static std::vector<std::vector<unsigned int> > mFromN(unsigned int m, unsigned int n);
        static void combinatorics(unsigned int offset,
                           unsigned int k,
                           std::vector<unsigned int> &combination,
                           std::vector<unsigned int> &mainList,
                           std::vector<std::vector<unsigned int> > &allCombinations);
        static void buildCombinations(const std::vector<const xAOD::Muon*> &muonsIn,
                               std::vector<Combination> &pairs,
                               std::vector<Combination> &quadruplets,
                               unsigned int nSelectedMuons);
        static bool passesQuadSelection(const std::vector<const xAOD::Muon*> &muonsIn);
        //-------------------------------------------------------------------------------------
        
    private:
        double m_ptCut;
        double m_etaCut;
        bool m_useV0Fitter;
        SG::ReadHandleKey<xAOD::MuonContainer>          m_muonCollectionKey;
        SG::ReadHandleKey<xAOD::TrackParticleContainer> m_TrkParticleCollection;
        ToolHandle < Trk::IVertexFitter > m_iVertexFitter;
        ToolHandle < Trk::IVertexFitter > m_iV0VertexFitter;
        ToolHandle < Trk::V0Tools > m_V0Tools;
        ToolHandle < Trk::ITrackSelectorTool > m_trkSelector;
        ToolHandle < InDet::VertexPointEstimator > m_vertexEstimator;
        SG::ReadHandleKey<xAOD::EventInfo> m_eventInfo_key{this, "EventInfo", "EventInfo", "Input event information"};
        SG::WriteDecorHandleKey<xAOD::MuonContainer> m_muonIndex{this, "muonIndexDec", "Muons.BPHY4MuonIndex"};

    };
} // end of namespace
#endif
