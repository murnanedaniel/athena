/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonLayerAmbiguitySolverTool.h"

#include "MuonLayerEvent/MuonCandidate.h"
#include "MuonSegment/MuonSegment.h"
#include "TrkEventPrimitives/FitQuality.h"

namespace Muon {

    MuonLayerAmbiguitySolverTool::MuonLayerAmbiguitySolverTool(const std::string& type, const std::string& name, const IInterface* parent) :
        AthAlgTool(type, name, parent) {
        declareInterface<IMuonLayerAmbiguitySolverTool>(this);
    }

    StatusCode MuonLayerAmbiguitySolverTool::initialize() {
        ATH_CHECK(m_segmentSelector.retrieve());
        ATH_CHECK(m_segmentMatchingTool.retrieve());
        ATH_CHECK(m_muonTrackBuilder.retrieve());
        ATH_CHECK(m_printer.retrieve());

        return StatusCode::SUCCESS;
    }

    void MuonLayerAmbiguitySolverTool::resolveOverlaps(const EventContext& ctx, const std::vector<MuonLayerRecoData>& allLayers,
                                                       std::vector<MuonCandidate>& resolvedCandidates) const {
        // re-organise data to allow hash based access, resolve small large overlaps
        std::vector<std::vector<MuonLayerIntersection> > muonLayerDataHashVec;
        buildLayerVec(ctx, allLayers, muonLayerDataHashVec);

        // build candidate by selecting seeds and extending them
        unsigned int nseeds = 0;                    // counter for number of seeds up to now
        std::set<const MuonSegment*> usedSegments;  // keep track of the segments already used
        std::vector<MuonStationIndex::StIndex> inverseSeedLayerOrder = {MuonStationIndex::BO, MuonStationIndex::BI, MuonStationIndex::BM,
                                                                        MuonStationIndex::EO, MuonStationIndex::EE, MuonStationIndex::EI,
                                                                        MuonStationIndex::EM};
        while (nseeds < m_maxSeeds) {
            // first get a seed
            MuonLayerIntersection layerIntersection;
            if (!getNextSeed(muonLayerDataHashVec, usedSegments, inverseSeedLayerOrder, layerIntersection)) {
                ATH_MSG_VERBOSE("No more seeds, total used seeds " << nseeds);
                break;
            }

            // create first candidate from seed and extend it
            std::vector<MuonLayerIntersection> layerIntersections = {layerIntersection};
            std::vector<MuonCandidate> candidates = {MuonCandidate(std::move(layerIntersections))};
            if (extendCandidatesWithLayers(ctx, candidates, muonLayerDataHashVec, inverseSeedLayerOrder)) {
                // add candidates to output list
                ATH_MSG_DEBUG(" Completed seed extension " << candidates.size());
                if (msgLvl(MSG::VERBOSE)) {
                    for (const auto& candidate : candidates) {
                        msg(MSG::VERBOSE) << " Candidate with layers " << candidate.layerIntersections.size();
                        for (const auto& entry : candidate.layerIntersections) {
                            msg(MSG::VERBOSE) << std::endl << "  " << m_printer->print(*entry.segment);
                        }
                    }
                    msg(MSG::VERBOSE) << endmsg;
                }

                // add all segments on the candidates to the exlusion list
                for (const auto& candidate : candidates) {
                    for (const auto& layer : candidate.layerIntersections) { usedSegments.insert(layer.segment.get()); }
                }
                resolvedCandidates.insert(resolvedCandidates.end(), std::make_move_iterator(candidates.begin()),
                                          std::make_move_iterator(candidates.end()));
            }
            ++nseeds;
        }

        ATH_MSG_DEBUG("Completed ambiguity solving using " << nseeds << " seeds, resulting in " << resolvedCandidates.size()
                                                           << " track candidates ");
    }

    bool MuonLayerAmbiguitySolverTool::extendCandidatesWithLayers(
        const EventContext& ctx, std::vector<MuonCandidate>& candidates,
        const std::vector<std::vector<MuonLayerIntersection> >& muonLayerDataHashVec,
        const std::vector<MuonStationIndex::StIndex>& inverseSeedLayerOrder) const {
        // break recursive call chain once we processed all layers
        if (inverseSeedLayerOrder.empty()) return true;

        ATH_MSG_VERBOSE("extendCandidates " << candidates.size() << " remaining layers " << inverseSeedLayerOrder.size());

        // get data in current layer
        MuonStationIndex::StIndex currentStIndex = inverseSeedLayerOrder.back();
        const std::vector<MuonLayerIntersection>& layerIntersections = muonLayerDataHashVec[currentStIndex];
        if (!layerIntersections.empty()) {
            // store new MuonCandidates
            std::vector<MuonCandidate> newCandidates;

            // loop over candidates
            for (MuonCandidate& candidate : candidates) {
                // if more than one segment is selected in the layer, create a new candidate
                unsigned int selectedSegmentsInLayer = 0;
                // loop over data in layer
                for (const MuonLayerIntersection& layerIntersection : layerIntersections) {
                    // match segment to candidate. Segment pairs with the same identifier are
                    // excluded from the beginning
                    if (match(ctx, candidate, layerIntersection)) {
                        // if first add to existing candidate, else create a new candidate
                        if (selectedSegmentsInLayer == 0) {
                            candidate.layerIntersections.push_back(layerIntersection);
                        } else {
                            MuonCandidate newCandidate = candidate;
                            /// replace the last intersection on the original candidate
                            newCandidate.layerIntersections.back() = layerIntersection;
                            newCandidates.emplace_back(std::move(newCandidate));
                        }
                        ++selectedSegmentsInLayer;
                    }
                }
            }
            // add new candidates to list
            if (!newCandidates.empty()) {
                ATH_MSG_VERBOSE("Found multiple solutions, add new candidates " << newCandidates.size());
                candidates.insert(candidates.end(), std::make_move_iterator(newCandidates.begin()),
                                  std::make_move_iterator(newCandidates.end()));
            }
        }

        // remove the current layer and call extendCandidatesWithLayers for the next layer
        std::vector<MuonStationIndex::StIndex> newInverseSeedLayerOrder = inverseSeedLayerOrder;
        newInverseSeedLayerOrder.pop_back();
        return extendCandidatesWithLayers(ctx, candidates, muonLayerDataHashVec, newInverseSeedLayerOrder);
    }

    bool MuonLayerAmbiguitySolverTool::match(const EventContext& ctx, const MuonCandidate& candidate,
                                             const MuonLayerIntersection& layerIntersection) const {
        // loop over layers and match each segment to the new one, if any fails, fail the combination
        for (const Muon::MuonLayerIntersection& layer : candidate.layerIntersections) {
            if (!m_segmentMatchingTool->match(ctx, *layer.segment, *layerIntersection.segment)) return false;
        }
        return true;
    }

    bool MuonLayerAmbiguitySolverTool::getNextSeed(const std::vector<std::vector<MuonLayerIntersection> >& muonLayerDataHashVec,
                                                   std::set<const MuonSegment*>& usedSegments,
                                                   std::vector<MuonStationIndex::StIndex>& inverseSeedLayerOrder,
                                                   MuonLayerIntersection& layerIntersection) const {
        ATH_MSG_VERBOSE("getNextSeed, remaining layers " << inverseSeedLayerOrder.size());
        // loop over the inverse seed layers
        std::vector<MuonStationIndex::StIndex>::const_reverse_iterator rit = inverseSeedLayerOrder.rbegin();
        std::vector<MuonStationIndex::StIndex>::const_reverse_iterator rit_end = inverseSeedLayerOrder.rend();
        for (; rit != rit_end; ++rit) {
            // loop over segments and find the next 'good' one that was not used yet
            for (const MuonLayerIntersection& muonLayerIntersection : muonLayerDataHashVec[*rit]) {
                /// select segment
                if (muonLayerIntersection.quality < m_seedQualityThreshold) continue;
                // only consider once
                const MuonSegment* segment = muonLayerIntersection.segment.get();
                if (usedSegments.count(segment)) continue;
                usedSegments.insert(segment);

                // return result
                layerIntersection = muonLayerIntersection;
                ATH_MSG_VERBOSE("Selected seed " << m_printer->print(*segment));
                return true;
            }
            // if we get here, we processed all the possible seeds in the layer so we can remove it from the list
            inverseSeedLayerOrder.pop_back();
        }
        return false;
    }

    void MuonLayerAmbiguitySolverTool::buildLayerVec(const EventContext& ctx, const std::vector<MuonLayerRecoData>& allLayers,
                                                     std::vector<std::vector<MuonLayerIntersection> >& muonLayerDataHashVec) const {
        // clear and resize hash vector, initialize with null_ptr
        muonLayerDataHashVec.clear();
        muonLayerDataHashVec.resize(MuonStationIndex::StIndexMax);

        // loop over layers
        for (const MuonLayerRecoData& layer : allLayers) {
            const MuonLayerSurface& layerSurface = layer.intersection.layerSurface;
            MuonStationIndex::StIndex stIndex = MuonStationIndex::toStationIndex(layerSurface.regionIndex, layerSurface.layerIndex);

            // create layer intersections
            std::vector<MuonLayerIntersection> layerIntersections;
            layerIntersections.reserve(layer.segments.size());
            for (const std::shared_ptr<const MuonSegment>& segment : layer.segments) {
                // get quality and skip bad ones
                int quality = m_segmentSelector->quality(*segment);
                if (quality < m_minSegmentQuality) continue;
                layerIntersections.emplace_back(layer.intersection, segment, quality);
            }

            // if there are no segments yet in the layer, directly add them
            if (muonLayerDataHashVec[stIndex].empty()) {
                muonLayerDataHashVec[stIndex] = std::move(layerIntersections);
            } else {
                // there are already segment, try resolving small/large overlaps
                resolveSmallLargeOverlaps(ctx, muonLayerDataHashVec[stIndex], layerIntersections);
            }

            // finally sort the segments
            std::stable_sort(
                muonLayerDataHashVec[stIndex].begin(), muonLayerDataHashVec[stIndex].end(),
                [](const Muon::MuonLayerIntersection& a, const Muon::MuonLayerIntersection& b) { return a.quality > b.quality; });
        }

        if (msgLvl(MSG::DEBUG)) {
            msg(MSG::DEBUG) << " Done building segment vector ";
            for (const auto& vec : muonLayerDataHashVec) {
                for (const auto& entry : vec) { msg(MSG::DEBUG) << std::endl << "  " << m_printer->print(*entry.segment); }
            }
            msg(MSG::DEBUG) << endmsg;
        }
    }

    void MuonLayerAmbiguitySolverTool::resolveSmallLargeOverlaps(const EventContext& ctx,
                                                                 std::vector<MuonLayerIntersection>& existingLayerIntersections,
                                                                 std::vector<MuonLayerIntersection>& newLayerIntersections) const {
        ATH_MSG_VERBOSE(" resolveSmallLargeOverlaps: existing " << existingLayerIntersections.size() << " new "
                                                                << newLayerIntersections.size());

        // keep track of segments that have been merged
        std::set<const MuonSegment*> combinedSegments;
        std::vector<MuonLayerIntersection> combinedIntersections;

        // loop over all permutations
        for (const MuonLayerIntersection& layerIntersection1 : existingLayerIntersections) {
            // get quality and skip bad ones
            if (layerIntersection1.quality < m_minSegmentQuality) continue;
            for (const MuonLayerIntersection& layerIntersection2 : newLayerIntersections) {
                // get quality and skip bad ones
                if (layerIntersection2.quality < m_minSegmentQuality) continue;

                // require at least one of the segments to be above seeding threshold
                if (layerIntersection1.quality < m_seedQualityThreshold && layerIntersection2.quality < m_seedQualityThreshold) continue;

                // match segments
                if (!m_segmentMatchingTool->match(ctx, *layerIntersection1.segment, *layerIntersection2.segment)) continue;

                // build new segment
                std::shared_ptr<const MuonSegment> newseg{
                    m_muonTrackBuilder->combineToSegment(ctx, *layerIntersection1.segment, *layerIntersection2.segment, nullptr)};
                if (!newseg) {
                    ATH_MSG_DEBUG(" Fit of combination of segments failed ");
                    continue;
                }

                // check fit quality
                const Trk::FitQuality* fq = newseg->fitQuality();
                if (!fq || fq->numberDoF() == 0) {
                    ATH_MSG_WARNING(" No fit quality, dropping segment ");
                    continue;
                }
                if (fq->chiSquared() / fq->numberDoF() > 2.5) {
                    ATH_MSG_DEBUG("Bad fit quality, dropping segment " << fq->chiSquared() / fq->numberDoF());
                    continue;
                }

                // check that quality of the combined segment is not worse that the original ones
                int qualitynew = m_segmentSelector->quality(*newseg);
                if (qualitynew < layerIntersection1.quality || qualitynew < layerIntersection2.quality) {
                    ATH_MSG_DEBUG("Quality got worse after combination: new " << qualitynew << " q1 " << layerIntersection1.quality
                                                                              << " q2 " << layerIntersection2.quality);
                    continue;
                }

                // select intersection closest to the IP and create new MuonLayerIntersection
                // get line from combined segment to the IP
                Amg::Vector3D direction = newseg->globalPosition().unit();
                // lambda to project the intersection positions on the line, treats the case where pointer is null
                auto getDistance = [](const MuonLayerIntersection& layerIntersection, const Amg::Vector3D& direction) {
                    if (!layerIntersection.intersection.trackParameters) return 1e9;
                    return layerIntersection.intersection.trackParameters->position().dot(direction);
                };
                double dist1 = getDistance(layerIntersection1, direction);
                double dist2 = getDistance(layerIntersection2, direction);
                if (dist1 < dist2)
                    combinedIntersections.emplace_back(layerIntersection1.intersection, newseg, qualitynew);
                else
                    combinedIntersections.emplace_back(layerIntersection2.intersection, newseg, qualitynew);

                ATH_MSG_DEBUG(" Combined segments " << std::endl
                                                    << " first   " << m_printer->print(*layerIntersection1.segment) << std::endl
                                                    << " second  " << m_printer->print(*layerIntersection2.segment) << std::endl
                                                    << " combined " << m_printer->print(*newseg));

                // add segments to exclusion list
                combinedSegments.insert(layerIntersection1.segment.get());
                combinedSegments.insert(layerIntersection2.segment.get());
            }
        }

        // lambda to loop over the input intersections and add them to the new list
        auto insert_intersection = [&combinedSegments](const Muon::MuonLayerIntersection& inter_sect) {
            return !combinedSegments.count(inter_sect.segment.get());
        };
        // do the insertion and swap the new vector with the existingLayerIntersections
        combinedIntersections.reserve(existingLayerIntersections.size() + newLayerIntersections.size());
        std::copy_if(std::make_move_iterator(existingLayerIntersections.begin()), std::make_move_iterator(existingLayerIntersections.end()),
                     std::back_inserter(combinedIntersections), insert_intersection);
        std::copy_if(std::make_move_iterator(newLayerIntersections.begin()), std::make_move_iterator(newLayerIntersections.end()),
                     std::back_inserter(combinedIntersections), insert_intersection);
        existingLayerIntersections = std::move(combinedIntersections);
    }
}  // namespace Muon
