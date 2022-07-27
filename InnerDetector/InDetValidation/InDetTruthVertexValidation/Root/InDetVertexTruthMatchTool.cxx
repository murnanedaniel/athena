/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetTruthVertexValidation/InDetVertexTruthMatchTool.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"

#include "InDetTruthVertexValidation/InDetVertexTruthMatchUtils.h"
using namespace InDetVertexTruthMatchUtils;

InDetVertexTruthMatchTool::InDetVertexTruthMatchTool( const std::string & name ) : asg::AsgTool(name) {
  declareProperty("trackMatchProb", m_trkMatchProb = 0.5 );
  declareProperty("vertexMatchWeight", m_vxMatchWeight = 0.7 );
  declareProperty("trackPtCut", m_trkPtCut = 500. );
}

StatusCode InDetVertexTruthMatchTool::initialize() {
  ATH_MSG_INFO("Initializing");

  return StatusCode::SUCCESS;
}

StatusCode InDetVertexTruthMatchTool::finalize()
{
  if (m_nBadLinks>0) {
     ATH_MSG_INFO("Vertex fraction with invalid track particle links: " << std::setw(8) << static_cast<double>(m_nVtxWithBadLinks) / m_nVtx << "; invalid track particle link fraction " << static_cast<double>(m_nBadLinks) / m_nLinks << ".");
  }
  return StatusCode::SUCCESS;
}

namespace {
//Helper methods for this file only

//output typedef
// this is defined like this in InDetTruthMatchUtils:
// typedef std::pair<ElementLink<xAOD::TruthEventBaseContainer>, float> VertexTruthMatchInfo;
// std::pair<ElementLink<>, T> is special templated pair in ElementLink.h even
// pair of link to a truth event with relative weight of matched tracks


// Create a truth map by decoarting truth particles with a back link to the truth event they live in
// Needed because the track->truth assoc gives us the particles but they don't store event normally
// Add as decoration to avoid full loop for every track ( this time only once per event )
// Use a vector so any number of truth event collections can be used at once -- but the pointers need to be valid
void createTruthMap(const std::vector<const xAOD::TruthEventBaseContainer *> & truthEventContainers ) {

  static const xAOD::TruthParticle::Decorator<ElementLink<xAOD::TruthEventBaseContainer> > backLinkDecor("TruthEventLink");

  for ( auto cit : truthEventContainers ) {

    const xAOD::TruthEventBaseContainer & truthEvents = *cit;

    for ( size_t i = 0; i < truthEvents.size(); ++i) {

      for ( auto & tkit : truthEvents[i]->truthParticleLinks() ) { //std::vector<ElementLink...

        const ElementLink<xAOD::TruthEventBaseContainer> elLink = ElementLink<xAOD::TruthEventBaseContainer>( truthEvents, i );

        if (elLink.isValid() && tkit.isValid()) {
          backLinkDecor(**tkit) = elLink;
        }

      }

    }

  }

}

void createTrackTruthMap(const std::vector<const xAOD::TruthEventBaseContainer *> & truthEventContainers,
                         const xAOD::TrackParticleContainer & trackParticleContainer,
                         float matchCut)
{

  createTruthMap(truthEventContainers);

  xAOD::TrackParticle::ConstAccessor<ElementLink<xAOD::TruthParticleContainer> > trk_truthPartAcc("truthParticleLink");
  xAOD::TrackParticle::ConstAccessor<float> trk_truthProbAcc("truthMatchProbability");
  static const xAOD::TruthParticle::Decorator<ElementLink<xAOD::TruthEventBaseContainer> > backLinkDecor("TruthEventLink");
  static const xAOD::TrackParticle::Decorator<ElementLink<xAOD::TruthEventBaseContainer> > trackLinkDecor("TrackEventLink");

  int nGood = 0;
  int nMatch = 0;
  int nLink = 0;
  for (auto trk : trackParticleContainer)
  {
    {
      nGood++;
      if (trk_truthPartAcc.isAvailable(*trk) && trk_truthProbAcc.isAvailable(*trk)
          && trk_truthPartAcc(*trk).isValid() && trk_truthProbAcc(*trk) >= matchCut)
      {
        nMatch++;
        const auto& truthParticle = trk_truthPartAcc(*trk);
        if (backLinkDecor.isAvailable(**truthParticle) && backLinkDecor(**truthParticle).isValid())
        {
          nLink++;
          trackLinkDecor(*trk) = backLinkDecor(**truthParticle);
        }
      }
    }
  }
  // won't compile, no idea why
  //ATH_MSG_DEBUG("Linked/Matched/Good/All: " << nLink << " / " << nMatch << " / " << nGood << " / " << trackParticleContainer.size());
}

//In the vector of match info, find the element corresponding to link and return its index; create a new one if necessary
size_t indexOfMatchInfo( std::vector<VertexTruthMatchInfo> & matches, ElementLink<xAOD::TruthEventBaseContainer> & link ) {
  for ( size_t i = 0; i < matches.size(); ++i ) {
    if ( link.key() == std::get<0>(matches[i]).key() && link.index() == std::get<0>(matches[i]).index() )
      return i;
  }
  // This is the first time we've seen this truth vertex, so make a new entry
  matches.emplace_back( link, 0., 0. );
  return matches.size() - 1;
}


//for sorting the container -> highest relative match weight first
bool compareMatchPair(const VertexTruthMatchInfo& a, const VertexTruthMatchInfo& b ) { return std::get<1>(a) > std::get<1>(b); }

}


const xAOD::TrackParticleContainer*
InDetVertexTruthMatchTool::findTrackParticleContainer( const xAOD::VertexContainer& vxContainer ) 
{
  for (auto vtx : vxContainer)
  {
    for (const ElementLink<xAOD::TrackParticleContainer>& tpLink : vtx->trackParticleLinks())
    {
      if (tpLink.isValid())
      {
        return tpLink.getStorableObjectPointer();
      }
    }
  }
  return nullptr;
}

StatusCode InDetVertexTruthMatchTool::matchVertices( const xAOD::VertexContainer & vxContainer ) const {

  ATH_MSG_DEBUG("Start vertex matching");
  if (vxContainer.empty() ||   // reject empty vertex containers
       (vxContainer.size() == 1 && vxContainer.at(0)->vertexType() == xAOD::VxType::NoVtx)){  // as well as containers containing only a dummy vertex
    ATH_MSG_DEBUG("No vertices to match.");
    return StatusCode::SUCCESS;
  }
  // Identify MC vertices to match to -- this is the collection for hard scatter
  const xAOD::TruthEventBaseContainer * truthEvents = nullptr;
  if ( evtStore()->contains<xAOD::TruthEventBaseContainer>( "TruthEvents" ) )
    ATH_CHECK( evtStore()->retrieve( truthEvents, "TruthEvents" ) );
  else
    ATH_CHECK( evtStore()->retrieve( truthEvents, "TruthEvent" ) );

  std::vector<const xAOD::TruthEventBaseContainer *> truthContainers;
  truthContainers.push_back( truthEvents );

  ATH_MSG_DEBUG("Found Hard Scatter collection");

  // These are the pile-up truth -- don't want to fail if they don't exist
  const xAOD::TruthEventBaseContainer * truthPileup = nullptr;
  if ( evtStore()->contains<xAOD::TruthEventBaseContainer>( "TruthPileupEvents" ) )
    ATH_CHECK( evtStore()->retrieve( truthPileup, "TruthPileupEvents" ) );
  if (truthPileup)
    truthContainers.push_back( truthPileup );

  ATH_MSG_DEBUG("Found Pileup collection");

  // Find the trackParticle container associated with our reconstructed vertices
  // We could pass this, but it would break the original interface...
  const xAOD::TrackParticleContainer* tkContainer = findTrackParticleContainer(vxContainer);
  if (!tkContainer)
  {
    ATH_MSG_WARNING("Vertex container has no vertices with valid TrackParticle links");
    return StatusCode::SUCCESS;
  }

  ATH_MSG_DEBUG("Found track collection");

  // create the particle links to events to avoid excessive looping
  // also decorate reconstructed tracks passing selection with truthEvent links
  createTrackTruthMap( truthContainers, *tkContainer, m_trkMatchProb );

  // Accessor for the links we just created
  xAOD::TruthParticle::ConstAccessor<ElementLink<xAOD::TruthEventBaseContainer> > backLinkDecor("TruthEventLink");

  //setup decorators for truth matching info
  static const xAOD::Vertex::Decorator<std::vector<VertexTruthMatchInfo> > matchInfoDecor("TruthEventMatchingInfos");
  static const xAOD::Vertex::Decorator<std::vector<VertexTruthMatchInfo> > rawMatchInfoDecor("TruthEventRawMatchingInfos");
  static const xAOD::Vertex::Decorator<VertexMatchType> matchTypeDecor("VertexMatchType");
  static const xAOD::Vertex::Decorator<std::vector<ElementLink<xAOD::VertexContainer> > > splitPartnerDecor("SplitPartners");
  static const xAOD::Vertex::Decorator<int> nHSTrkDecor("nHSTrk");

  //setup accessors
  // can switch to built in method in xAOD::Vertex once don't have to deal with changing names anymore
  xAOD::Vertex::ConstAccessor<xAOD::Vertex::TrackParticleLinks_t> trkAcc("trackParticleLinks");
  xAOD::Vertex::ConstAccessor<std::vector<float> > weightAcc("trackWeights");

  xAOD::TrackParticle::ConstAccessor<ElementLink<xAOD::TruthParticleContainer> > trk_truthPartAcc("truthParticleLink");
  xAOD::TrackParticle::ConstAccessor<float> trk_truthProbAcc("truthMatchProbability");

  static const xAOD::TrackParticle::Decorator<ElementLink<xAOD::VertexContainer> > trk_recoVtx("RecoVertex");
  static const xAOD::TrackParticle::Decorator<float> trk_wtVtx("WeightVertex");

  //some variables to store
  size_t ntracks;
  xAOD::VxType::VertexType vxType;

  ATH_MSG_DEBUG("Starting Loop on Vertices");

  //=============================================================================
  //First loop over vertices: get tracks, then TruthParticles, and store relative
  //weights of contribution from each TruthEvent
  //=============================================================================
  size_t vxEntry = 0;
  unsigned int n_bad_links = 0;
  unsigned int n_links = 0;
  unsigned int n_vx_with_bad_links = 0;

  for ( auto vxit : vxContainer.stdcont() ) {
    vxEntry++;
    vxType = static_cast<xAOD::VxType::VertexType>( vxit->vertexType() );
    if (vxType == xAOD::VxType::NoVtx) {
      //skip dummy vertices -> match info will be empty vector if someone tries to access later
      //type will be set to dummy
      ATH_MSG_DEBUG("FOUND xAOD::VxType::NoVtx");
      continue;
    }


    //create the vector we will add as matching info decoration later
    std::vector<VertexTruthMatchInfo> matchinfo;
    std::vector<VertexTruthMatchInfo> rawMatchinfo; //not normalized to one for each reco vertex

    //if don't have track particles
    if (!trkAcc.isAvailable(*vxit) || !weightAcc.isAvailable(*vxit) ) {
      ATH_MSG_DEBUG("trackParticles or trackWeights not available, setting fake");
      // Add invalid link for fakes
      matchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), 1., 0. );
      matchInfoDecor( *vxit ) = matchinfo;
      rawMatchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), 1., 0. );
      rawMatchInfoDecor( *vxit ) = rawMatchinfo;
      nHSTrkDecor( *vxit ) = 0;
      continue;
    }

    //things we need to do the matching
    const xAOD::Vertex::TrackParticleLinks_t & trkParts = trkAcc( *vxit );
    ntracks = trkParts.size();
    const std::vector<float> & trkWeights = weightAcc( *vxit );

    //double check
    if ( trkWeights.size() != ntracks ) {
      ATH_MSG_DEBUG("Vertex without same number of tracks and trackWeights, setting fake");
      matchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), 1., 0. );
      matchInfoDecor( *vxit ) = matchinfo;
      rawMatchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), 1., 0. );
      rawMatchInfoDecor( *vxit ) = rawMatchinfo;
      nHSTrkDecor( *vxit ) = 0;
      continue;
    }

    ATH_MSG_DEBUG("Matching new vertex at (" << vxit->x() << ", " << vxit->y() << ", " << vxit->z() << ")" << " with " << ntracks << " tracks, at index: " << vxit->index());

    float totalWeight = 0.;
    float totalFake = 0.;
    int nHSTrk = 0;

    unsigned vx_n_bad_links = 0;
    //loop element link to track particle
    for ( size_t t = 0; t < ntracks; ++t ) {
      if (!trkParts[t].isValid()) {
         ++vx_n_bad_links;
         continue;
      }
      const xAOD::TrackParticle & trk = **trkParts[t];

      totalWeight += trkWeights[t];
      trk_recoVtx(trk) = ElementLink<xAOD::VertexContainer>(vxContainer, vxEntry - 1);
      trk_wtVtx(trk) = trkWeights[t];

      const ElementLink<xAOD::TruthParticleContainer> & truthPartLink = trk_truthPartAcc( trk );
      float prob = trk_truthProbAcc( trk );

      if (!truthPartLink.isValid()) continue;

      if (prob > m_trkMatchProb) {
        const xAOD::TruthParticle & truthPart = **truthPartLink;
        //check if the truth particle is "good"
        if ( pass( truthPart) ) {
          ElementLink<xAOD::TruthEventBaseContainer> match = backLinkDecor( truthPart );
          //check we have an actual link
          if ( match.isValid() ) {
            size_t matchIdx = indexOfMatchInfo( matchinfo, match );
            std::get<1>(matchinfo[matchIdx]) += trkWeights[t];
            std::get<2>(matchinfo[matchIdx]) += (trk.pt()/1000.) * (trk.pt()/1000.) * trkWeights[t];
            matchIdx = indexOfMatchInfo( rawMatchinfo, match );
            std::get<1>(rawMatchinfo[matchIdx]) += trkWeights[t];
            std::get<2>(rawMatchinfo[matchIdx]) += (trk.pt()/1000.) * (trk.pt()/1000.) * trkWeights[t];
            if((*match)->type() == xAOD::Type::TruthEvent && match.index() == 0) nHSTrk++;
          } else {
            totalFake += trkWeights[t];
          }

        } else {
          //truth particle failed cuts -> add to fakes
          totalFake += trkWeights[t];
        }
      } else {
        //not valid or low matching probability -> add to fakes
        totalFake += trkWeights[t];
      }
    }//end loop over tracks in vertex
    n_links     += ntracks;
    n_bad_links += vx_n_bad_links;
    if (vx_n_bad_links>0) {
       ++n_vx_with_bad_links;
    }

    //finalize the match info vector
    if ( totalWeight < 1e-12 ) { // in case we don't have any passing tracks we want to make sure labeled fake
      ATH_MSG_DEBUG("   Declaring vertex fully fake (no passing tracks included)");
      totalWeight = 1.;
      totalFake = 1.;
    }
    if ( totalFake > 0. )
    {
      matchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), totalFake, 0. );
      rawMatchinfo.emplace_back( ElementLink<xAOD::TruthEventBaseContainer>(), totalFake, 0. );
    }

    for ( auto & mit : matchinfo ) {
      std::get<1>(mit) /= totalWeight;
    }
    std::sort( matchinfo.begin(), matchinfo.end(), compareMatchPair );
    std::sort( rawMatchinfo.begin(), rawMatchinfo.end(), compareMatchPair );
    matchInfoDecor( *vxit ) = matchinfo;
    rawMatchInfoDecor( *vxit ) = rawMatchinfo;
    nHSTrkDecor( *vxit ) = nHSTrk;
  }
  m_nVtx             += vxContainer.stdcont().size();
  m_nVtxWithBadLinks += n_vx_with_bad_links;
  m_nBadLinks        += n_bad_links;
  m_nLinks           += n_links;

  //After first loop, all vertices have been decorated with their vector of match info (link to TruthEvent paired with weight)
  //now we want to use that information from the whole collection to assign types

  //keep track of whether a type is assigned
  //useful since looking for splits involves a double loop, and then setting types ahead in the collection
  std::vector<bool> assignedType( vxContainer.size(), false );

  for ( size_t i = 0; i < vxContainer.size(); ++i ) {

    if ( assignedType[i] ) continue; // make sure we don't reclassify vertices already found in the split loop below

    std::vector<VertexTruthMatchInfo> & info = matchInfoDecor( *vxContainer[i] );
    if (info.empty()) {
      matchTypeDecor( *vxContainer[i] ) = DUMMY;
    } else if ( !std::get<0>(info[0]).isValid() ) {
      matchTypeDecor( *vxContainer[i] ) = FAKE;
    } else if ( std::get<1>(info[0]) > m_vxMatchWeight ) {
      matchTypeDecor( *vxContainer[i] ) = MATCHED;
    } else {
      matchTypeDecor( *vxContainer[i] ) = MERGED;
    }

    //check for splitting
    if ( matchTypeDecor( *vxContainer[i] ) == MATCHED || matchTypeDecor( *vxContainer[i] ) == MERGED ) {
      std::vector<size_t> foundSplits;
      for ( size_t j = i + 1; j < vxContainer.size(); ++j ) {
        std::vector<VertexTruthMatchInfo> & info2 = matchInfoDecor( *vxContainer[j] );
        //check second vertex is not dummy or fake, and that it has same elementlink as first vertex
        //equality test is in code but doesnt seem to work for ElementLinks that I have?
        //so i am just checking that the contianer key hash and the index are the same
        if (matchTypeDecor( *vxContainer[j] ) == FAKE || matchTypeDecor( *vxContainer[j] ) == DUMMY) continue;
        if (!info2.empty() && std::get<0>(info2[0]).isValid() && std::get<0>(info[0]).key() == std::get<0>(info2[0]).key() && std::get<0>(info[0]).index() == std::get<0>(info2[0]).index() ) {
          //add split links; first between first one found and newest one
          splitPartnerDecor( *vxContainer[i] ).emplace_back( vxContainer, j );
          splitPartnerDecor( *vxContainer[j] ).emplace_back( vxContainer, i );
          //then between any others we found along the way
          for ( auto k : foundSplits ) { //k is a size_t in the vector of splits
            splitPartnerDecor( *vxContainer[k] ).emplace_back( vxContainer, j );
            splitPartnerDecor( *vxContainer[j] ).emplace_back( vxContainer, k );
          }
          //then keep track that we found this one
          foundSplits.push_back(j);
        } //if the two vertices match to same TruthEvent
      }//inner loop over vertices

      // Correct labelling of split vertices - keep highest sumpt2 vertex labelled as matched/merged
      float maxSumpT2 = std::get<2>( matchInfoDecor( *vxContainer[i] )[0] );
      size_t indexOfMax = i;
      for ( auto l : foundSplits ) {
        if ( std::get<2>( matchInfoDecor( *vxContainer[l] )[0] ) > maxSumpT2 ){
          maxSumpT2 = std::get<2>( matchInfoDecor( *vxContainer[l] )[0] );
          indexOfMax = l;
        } else {
          matchTypeDecor( *vxContainer[l] ) = SPLIT;
          assignedType[l] = true;
        }
      }
      if ( indexOfMax!=i ) matchTypeDecor( *vxContainer[i] ) = SPLIT;
    } //if matched or merged
  } //outer loop

  //DEBUG MATCHING
  if (msgLvl(MSG::DEBUG)) {
    for (const auto &vxit : vxContainer.stdcont() ) {
      ATH_MSG_DEBUG("Matched vertex (index " << (*vxit).index() << ") to type " << matchTypeDecor(*vxit) << " with following info of size " << matchInfoDecor(*vxit).size() << ":");
      for (const auto &vit : matchInfoDecor(*vxit) ) {
        if ( std::get<0>(vit).isValid() ) {
          ATH_MSG_DEBUG("  GenEvent type " << (* std::get<0>(vit))->type() << ", index " << std::get<0>(vit).index() << " with relative weight " << std::get<1>(vit) );
        } else {
          ATH_MSG_DEBUG("  Fakes with relative weight " << std::get<1>(vit) );
        }
      }
      if (matchTypeDecor(*vxit) == SPLIT) {
        ATH_MSG_DEBUG("  Split partners are:");
        for (const auto &split : splitPartnerDecor( *vxit ) ) {
          if ( split.isValid() )
            ATH_MSG_DEBUG("    Vertex " << split.index());
          else
            ATH_MSG_DEBUG("    ERROR");
        }
      }
    }
  }

  return StatusCode::SUCCESS;

}


//Set up any cuts on either the tracks or truth particles to allow here
//A failing track is removed from consideration entirely
//If a passing track matches to a failing truth particle it will be considered "fake"

bool InDetVertexTruthMatchTool::pass( const xAOD::TruthParticle & truthPart ) const {

  //remove the registered secondaries
  return truthPart.pt() >= m_trkPtCut;

}
