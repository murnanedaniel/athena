/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//////////////////////////////////////////////////////////////////////////////
// MuonCombinedFitTagTool
//  AlgTool performing combined fit of ID and MS tracks (Muid)
//  A CombinedFitTag is added to the InDetCandidate object.
//
//////////////////////////////////////////////////////////////////////////////

#include "MuonCombinedFitTagTool.h"

#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "MuonCombinedEvent/CombinedFitTag.h"
#include "MuonCombinedEvent/InDetCandidate.h"
#include "MuonCombinedEvent/InDetCandidateToTagMap.h"
#include "MuonCombinedEvent/MuonCandidate.h"
#include "MuonRIO_OnTrack/MdtDriftCircleOnTrack.h"
#include "TrkEventUtils/IdentifierExtractor.h"
#include "TrkMaterialOnTrack/MaterialEffectsOnTrack.h"
#include "TrkMaterialOnTrack/ScatteringAngles.h"
#include "TrkTrack/TrackStateOnSurface.h"
#include "TrkTrackSummary/TrackSummary.h"
#include "muonEvent/CaloEnergy.h"
#include "xAODTracking/Vertex.h"

namespace {
    constexpr double probCut = .00001;  // cut on max probability: below this cut, we don't attempt to form a combined track unless no
                                       // combined track has yet been successfully created
    inline double chi2 (const Trk::FitQuality* fitQuality) {
        return !fitQuality || !fitQuality->numberDoF() ? 1.e12 : fitQuality->chiSquared() / fitQuality->doubleNumberDoF();
    } 

}
namespace MuonCombined {

    MuonCombinedFitTagTool::MuonCombinedFitTagTool(const std::string& type, const std::string& name, const IInterface* parent) :
        AthAlgTool(type, name, parent) {
        declareInterface<IMuonCombinedTagTool>(this);
    }

    StatusCode MuonCombinedFitTagTool::initialize() {
        ATH_MSG_INFO("Initializing MuonCombinedFitTagTool");

        ATH_CHECK(m_printer.retrieve());
        ATH_CHECK(m_trackBuilder.retrieve());
        if (!m_outwardsBuilder.empty()) ATH_CHECK(m_outwardsBuilder.retrieve());
        ATH_CHECK(m_trackQuery.retrieve());
        ATH_CHECK(m_momentumBalanceTool.retrieve());
        if (!m_muonRecovery.empty()) ATH_CHECK(m_muonRecovery.retrieve());
        ATH_CHECK(m_matchQuality.retrieve());
        ATH_CHECK(m_trackScoringTool.retrieve());
        /// handle to the magnetic field cache
        ATH_CHECK(m_fieldCacheCondObjInputKey.initialize());

        // The trigger doesn't use the vertex information
        if (!m_vertexKey.empty()) ATH_CHECK(m_vertexKey.initialize());

        return StatusCode::SUCCESS;
    }

    void MuonCombinedFitTagTool::combine(const MuonCandidate& muonCandidate, const std::vector<const InDetCandidate*>& indetCandidates,
                                         InDetCandidateToTagMap& tagMap, TrackCollection* combTracks, TrackCollection* METracks,
                                         const EventContext& ctx) const {
        ATH_MSG_DEBUG("muon candidate: " << muonCandidate.toString());

        std::unique_ptr<CombinedFitTag> bestTag, currentTag;
        std::unique_ptr<Trk::Track> bestCombTrack, bestMETrack, combinedTrack, METrack;
        const InDetCandidate* bestCandidate = nullptr;

        // map of ID candidates by max probability of match (based on match chi2 at IP and MS entrance)
        using InDetProbMatch = std::pair<double, const InDetCandidate*>;
        std::vector<InDetProbMatch> sortedInDetCandidates;
        sortedInDetCandidates.reserve(indetCandidates.size());
        // loop over ID candidates
        for (const MuonCombined::InDetCandidate* idTP : indetCandidates) {
            const Trk::Track* id_track = idTP->indetTrackParticle().track();
            double outerMatchProb = m_matchQuality->outerMatchProbability(*id_track, muonCandidate.muonSpectrometerTrack(), ctx);
            double innerMatchProb = -1;
            if (muonCandidate.extrapolatedTrack())
                innerMatchProb = m_matchQuality->innerMatchProbability(*id_track, *muonCandidate.extrapolatedTrack(), ctx);
            const double maxProb = std::max(outerMatchProb, innerMatchProb);
            sortedInDetCandidates.emplace_back(maxProb, idTP);
        }
        /// Sort such that the track with the largest probability comes first
        std::sort(sortedInDetCandidates.begin(), sortedInDetCandidates.end(),
                  [](const InDetProbMatch& a, const InDetProbMatch& b) { return a.first > b.first; });

        bool fitBadMatches = false;
        for (const InDetProbMatch& cand_prob : sortedInDetCandidates) {
            ATH_MSG_DEBUG("in det candidate prob: " << cand_prob.first);
            if (cand_prob.first < probCut && !fitBadMatches) {
                if (!bestCombTrack) {
                    ATH_MSG_DEBUG("no combined track yet, keep fitting");
                    fitBadMatches = true;
                } else {
                    ATH_MSG_DEBUG("combined track found, we're done here");
                    break;
                }
            }
            const Trk::Track* id_track = cand_prob.second->indetTrackParticle().track();
            ATH_MSG_DEBUG("Doing combined fit with ID track " << cand_prob.second->toString());
            ATH_MSG_DEBUG("Doing combined fit with MS track " << muonCandidate.toString());

            // fit the combined ID-MS track
            combinedTrack = buildCombinedTrack(ctx, *id_track, muonCandidate.muonSpectrometerTrack(), muonCandidate.extrapolatedTrack());
            if (!combinedTrack) {
                ATH_MSG_DEBUG("Combination fit failed");
                continue;
            }

            if (msgLevel() >= MSG::DEBUG) {
                dumpCaloEloss(ctx, combinedTrack.get(), "Combined Track ");
                dumpCaloEloss(ctx, muonCandidate.extrapolatedTrack(), "Extrapolated Track ");
            }

            // calculate track score
            Trk::TrackScore score = m_trackScoringTool->score(*combinedTrack, true);

            // add fit info into tag object
            currentTag = std::make_unique<CombinedFitTag>(xAOD::Muon::MuidCo, muonCandidate, score);

            // re-fit standalone track (if needed) and store output into tag object
            METrack = evaluateMatchProperties(ctx, combinedTrack.get(), *currentTag, cand_prob.second->indetTrackParticle());

            // select the best combined track
            if (!bestCandidate || bestMatchChooser(*cand_prob.second, *currentTag, *combinedTrack, METrack.get(), *bestCandidate, *bestTag,
                                                   *bestCombTrack, bestMETrack.get())) {
                bestCandidate = cand_prob.second;
                bestTag.swap(currentTag);
                bestCombTrack.swap(combinedTrack);
                bestMETrack.swap(METrack);
            }
        }
        /// try recovery
        if (!bestCandidate && !m_muonRecovery.empty()) {
            for (const InDetProbMatch& cand_prob : sortedInDetCandidates) {
                const Trk::Track* id_track = cand_prob.second->indetTrackParticle().track();
                combinedTrack = m_muonRecovery->recoverableMatch(*id_track, muonCandidate.muonSpectrometerTrack(), ctx);
                if (combinedTrack && combinedTrackQualityCheck(ctx, *combinedTrack, *id_track)) {
                    combinedTrack->info().addPatternReco(id_track->info());
                    combinedTrack->info().addPatternReco(muonCandidate.muonSpectrometerTrack().info());
                    combinedTrack->info().setParticleHypothesis(Trk::muon);
                    combinedTrack->info().setPatternRecognitionInfo(Trk::TrackInfo::MuidCombined);
                    // calculate track score
                    Trk::TrackScore score = m_trackScoringTool->score(*combinedTrack, true);

                    // add fit info into tag object
                    currentTag = std::make_unique<CombinedFitTag>(xAOD::Muon::MuidCo, muonCandidate, score);

                    if (msgLevel() >= MSG::DEBUG) {
                        dumpCaloEloss(ctx, combinedTrack.get(), "Recovery Combined Track ");
                        dumpCaloEloss(ctx, muonCandidate.extrapolatedTrack(), "Recovery Extrapolated Track ");
                    }

                    // re-fit standalone track (if needed) and store output into tag object
                    METrack = evaluateMatchProperties(ctx, combinedTrack.get(), *currentTag, cand_prob.second->indetTrackParticle());

                    // select the best combined track
                    if (!bestCandidate || bestMatchChooser(*cand_prob.second, *currentTag, *combinedTrack, METrack.get(), *bestCandidate,
                                                           *bestTag, *bestCombTrack, bestMETrack.get())) {
                        bestCandidate = cand_prob.second;
                        bestTag.swap(currentTag);
                        bestCombTrack.swap(combinedTrack);
                        bestMETrack.swap(METrack);
                    }
                }
            }
        }

        if (bestCandidate) {
            // take the best MS Track, first the update extrapolated, than the extrapolated, last the spectrometer track
            if (msgLevel() >= MSG::DEBUG && bestMETrack) {
                dumpCaloEloss(ctx, bestCombTrack.get(), " bestCandidate Combined Track ");
                dumpCaloEloss(ctx, bestMETrack.get(), " bestCandidate Extrapolated Track ");
            }
            ATH_MSG_DEBUG("Final combined muon: " << m_printer->print(*bestCombTrack));
            ATH_MSG_DEBUG(m_printer->printStations(*bestCombTrack));
            ATH_MSG_DEBUG("Combined Muon with ID " << m_printer->print(bestCandidate->indetTrackParticle().perigeeParameters())
                                                   << " match chi2 " << bestTag->matchChi2());
            combTracks->push_back(std::move(bestCombTrack));
            ElementLink<TrackCollection> comblink(*combTracks, combTracks->size() - 1);
            bestTag->setCombinedTrackLink(comblink);
            if (bestMETrack) {
                METracks->push_back(std::move(bestMETrack));
                ElementLink<TrackCollection> melink(*METracks, METracks->size() - 1);
                bestTag->setUpdatedExtrapolatedTrackLink(melink);
            } else {
                bestTag->setUpdatedExtrapolatedTrackLink(ElementLink<TrackCollection>());
            }
            tagMap.addEntry(bestCandidate, bestTag.release());
        }
    }

    std::unique_ptr<Trk::Track> MuonCombinedFitTagTool::buildCombinedTrack(const EventContext& ctx,
                                                                           const Trk::Track& indetTrack,
                                                                           const Trk::Track& spectrometerTrack,
                                                                           const Trk::Track* extrapolatedTrack) const {
        // if no extrapolation is available
        if (!extrapolatedTrack) extrapolatedTrack = &spectrometerTrack;
        // build and fit the combined track
        std::unique_ptr<Trk::Track> combinedTrack;
        double combinedFitChi2 = 9999.;
        if (!m_trackBuilder.empty()) {
            combinedTrack = m_trackBuilder->combinedFit(ctx, indetTrack, *extrapolatedTrack, spectrometerTrack);
            if (combinedTrack && combinedTrack->fitQuality()) {
                combinedTrack->info().addPatternReco(extrapolatedTrack->info());
                combinedFitChi2 = chi2(combinedTrack->fitQuality());
            }
        }
        if (combinedFitChi2 > m_badFitChi2 && !m_outwardsBuilder.empty()) {
            std::unique_ptr<Trk::Track> outwardsTrack(
                m_outwardsBuilder->combinedFit(ctx, indetTrack, *extrapolatedTrack, spectrometerTrack));
            if (outwardsTrack && chi2(outwardsTrack->fitQuality()) < combinedFitChi2) {
                ATH_MSG_VERBOSE("buildCombinedTrack: choose outwards track");
                outwardsTrack->info().addPatternReco(spectrometerTrack.info());
                combinedTrack.swap(outwardsTrack);
            }
        }

        // filter out rubbish fits
        if (combinedTrack && combinedTrackQualityCheck(ctx, *combinedTrack, indetTrack)) {
            combinedTrack->info().addPatternReco(indetTrack.info());
            combinedTrack->info().setParticleHypothesis(Trk::muon);
            combinedTrack->info().setPatternRecognitionInfo(Trk::TrackInfo::MuidCombined);
            return combinedTrack;
        }
        return nullptr;
    }

    bool MuonCombinedFitTagTool::combinedTrackQualityCheck(const EventContext& ctx,
                                                           const Trk::Track& combinedTrack, 
                                                           const Trk::Track& indetTrack) const {
        // require calo correctly associated to track
        if (!m_trackQuery->isCaloAssociated(combinedTrack, ctx)) {
            ATH_MSG_DEBUG(" No Calorimeter CaloDeposit found on combined track ");
            return false;
        }
        // loose cut on momentumBalanceSignificance
        double significance = m_momentumBalanceTool->momentumBalanceSignificance(combinedTrack);
        if (std::abs(significance) > m_momentumBalanceCut) {
            ATH_MSG_DEBUG(" combinedTrackQualityCheck fails with momentumBalanceSignificance " << significance);
            return false;
        }

        // loose cut on indet/combined q/p pull (not applicable to indet line fit)
        if (!indetTrack.info().trackProperties(Trk::TrackInfo::StraightTrack)) {
            const Trk::Perigee* combinedPerigee = combinedTrack.perigeeParameters();
            const Trk::Perigee* indetPerigee = indetTrack.perigeeParameters();
            if (combinedPerigee->covariance() && indetPerigee->covariance()) {
                const double dpOverP2 = Amg::error(*combinedPerigee->covariance(), Trk::qOverP) * combinedPerigee->momentum().mag2();
                if (dpOverP2 < 1.E-6) {
                    // fail with unphysical momentum covariance
                    ATH_MSG_DEBUG("combinedTrackQualityCheck: fail with unphysical momentum covariance");
                    return false;
                }
                const double sigma = Amg::error(*indetPerigee->covariance(), Trk::qOverP);
                const double pull = (combinedPerigee->parameters()[Trk::qOverP] - indetPerigee->parameters()[Trk::qOverP]) / sigma;

                if (std::abs(pull) > m_indetPullCut) {
                    // fail with too high momentum pull
                    ATH_MSG_DEBUG("combinedTrackQualityCheck: fail with momentum pull above cut: "
                                  << pull << " pid " << 1. / indetPerigee->parameters()[Trk::qOverP] << " pcb "
                                  << 1. / combinedPerigee->parameters()[Trk::qOverP] << " 1./sigma " << 1. / sigma);
                    return false;
                }
            } else
                return false;
        }
        return true;
    }

    std::unique_ptr<Trk::Track> MuonCombinedFitTagTool::evaluateMatchProperties(const EventContext& ctx,
                                                                                const Trk::Track* combinedTrack, CombinedFitTag& tag,
                                                                                const xAOD::TrackParticle& idTrackParticle) const {
        const Trk::Track& idTrack = *idTrackParticle.track();
        // evaluate field integral and momentum balance significance for combined track
        tag.fieldIntegral(m_trackQuery->fieldIntegral(*combinedTrack, ctx));
        tag.momentumBalanceSignificance(m_momentumBalanceTool->momentumBalanceSignificance(*combinedTrack));

        if (tag.muonCandidate().extrapolatedTrack()) {
            std::pair<int, std::pair<double, double> > aTriad =
                m_matchQuality->innerMatchAll(idTrack, *tag.muonCandidate().extrapolatedTrack(), ctx);
            const int matchDoF = aTriad.first;
            const double matchChi2 = aTriad.second.first;
            const double matchProb = aTriad.second.second;
            // store the inner matching quality in the tag object
            tag.innerMatch(matchChi2, matchDoF, matchProb);
            ATH_MSG_DEBUG(" extrapolatedTrack innerMatch " << matchChi2);
        }

        // refit extrapolated from combined track (i.e. after cleaning)
        std::unique_ptr<Trk::Track> refittedExtrapolatedTrack;
        bool dorefit = true;

        // no SA refit for Toroid off
        MagField::AtlasFieldCache fieldCache;
        // Get field cache object

        SG::ReadCondHandle<AtlasFieldCacheCondObj> readHandle{m_fieldCacheCondObjInputKey, ctx};
        const AtlasFieldCacheCondObj* fieldCondObj{*readHandle};

        if (!fieldCondObj) {
            ATH_MSG_ERROR("Failed to retrieve AtlasFieldCacheCondObj with key " << m_fieldCacheCondObjInputKey.key());
            return nullptr;
        }
        fieldCondObj->getInitializedCache(fieldCache);
        if (!fieldCache.toroidOn()) dorefit = false;

        Amg::Vector3D origin{0., 0., 0.};

        const xAOD::Vertex* matchedVertex{nullptr};
        if (!m_vertexKey.empty()) {
            SG::ReadHandle<xAOD::VertexContainer> vertices{m_vertexKey, ctx};
            if (vertices.isValid()) {
                for (const xAOD::Vertex* vx : *vertices) {
                    for (const auto& tpLink : vx->trackParticleLinks()) {
                        if (*tpLink == &idTrackParticle) {
                            matchedVertex = vx;
                            break;
                        }
                    }
                    if (matchedVertex) break;
                }
            }
        }
        if (matchedVertex) {
            origin = Amg::Vector3D{matchedVertex->x(), matchedVertex->y(), matchedVertex->z()};
            ATH_MSG_DEBUG(" found matched vertex  bs " << origin);
        } else {
            //    take for beamspot point of closest approach of ID track in  x y z
            origin[Amg::x] = -idTrackParticle.d0() * std::sin(idTrackParticle.phi()) + idTrackParticle.vx();
            origin[Amg::y] = idTrackParticle.d0() * std::cos(idTrackParticle.phi()) + idTrackParticle.vy();
            origin[Amg::z] = idTrackParticle.z0() + idTrackParticle.vz();
            ATH_MSG_DEBUG(" NO matched vertex  take track perigee  " << origin);
        }

        ATH_MSG_DEBUG(" refit SA track " << dorefit);
        if (dorefit) {
            if (!m_trackBuilder.empty()) refittedExtrapolatedTrack = m_trackBuilder->standaloneRefit(ctx, *combinedTrack, origin);
            if (!refittedExtrapolatedTrack && !m_outwardsBuilder.empty())
                refittedExtrapolatedTrack = m_outwardsBuilder->standaloneRefit(ctx, *combinedTrack, origin);
        }
        // include vertex region pseudo for extrapolation failure
        unsigned numberPseudo =
            tag.muonCandidate().extrapolatedTrack() ? m_trackQuery->numberPseudoMeasurements(*tag.muonCandidate().extrapolatedTrack()) : 1;

        // get track quality and store
        if (refittedExtrapolatedTrack) {
            std::pair<int, std::pair<double, double> > aTriad = m_matchQuality->innerMatchAll(idTrack, *refittedExtrapolatedTrack, ctx);
            const int matchDoF = aTriad.first;
            const double matchChi2 = aTriad.second.first;
            const double matchProb = aTriad.second.second;

            // store the inner matching quality in the tag object
            tag.innerMatch(matchChi2, matchDoF, matchProb);
            ATH_MSG_DEBUG(" refittedExtrapolatedTrack innerMatch " << matchChi2);

            // print comparison with original track
            if (tag.muonCandidate().extrapolatedTrack()) {
                double oldmatchChi2 = m_matchQuality->innerMatchChi2(idTrack, *tag.muonCandidate().extrapolatedTrack(), ctx);

                ATH_MSG_VERBOSE("evaluateMatchProperties: chi2 re-evaluated from " << oldmatchChi2 << " to " << matchChi2);

                if (matchChi2 > 1.1 * oldmatchChi2)
                    ATH_MSG_DEBUG("evaluateMatchProperties: chi2 got worse: from " << oldmatchChi2 << " to " << matchChi2);
            } else
                ATH_MSG_VERBOSE("evaluateMatchProperties: added new extrapolated track with chi2 " << matchChi2);

        } else if (!numberPseudo) {
            // failed re-evaluation of match chi2
            ATH_MSG_DEBUG("evaluateMatchProperties: fail re-evaluation of match chi2");
        }
        return refittedExtrapolatedTrack;
    }

    void MuonCombinedFitTagTool::dumpCaloEloss(const EventContext& ctx, const Trk::Track* inTrack, const std::string& txt) const {
        // will refit if extrapolated track was definitely bad
        if (!inTrack) return;
        if (!m_trackQuery->isCaloAssociated(*inTrack, ctx)) {
            ATH_MSG_DEBUG(txt << " no TSOS in Calorimeter ");
            return;
        }
        const Trk::Track& originalTrack = *inTrack;
        const CaloEnergy* caloEnergy = m_trackQuery->caloEnergy(originalTrack);
        if (caloEnergy) {
            ATH_MSG_DEBUG(txt << " Calorimeter Eloss " << caloEnergy->deltaE());
        } else {
            ATH_MSG_DEBUG(txt << " No Calorimeter Eloss");
        }

        const Trk::TrackStates* trackTSOS = inTrack->trackStateOnSurfaces();

        double Eloss{0.}, idEloss{0.}, caloEloss{0.}, msEloss{0.}, deltaP{0.},
               pcalo{0.}, pstart{0.},eta{0.}, pMuonEntry{0.};
        for (const Trk::TrackStateOnSurface* m : *trackTSOS) {
            const Trk::MeasurementBase* mot = m->measurementOnTrack();
            if (m->trackParameters()) pMuonEntry = m->trackParameters()->momentum().mag();
            if (mot) {
                Identifier id = Trk::IdentifierExtractor::extract(mot);
                if (id.is_valid()) {
                    // skip after first Muon hit
                    if (m_idHelperSvc->isMuon(id)) break;
                }
            }
            if (pstart == 0 && m->trackParameters()) {
                pstart = m->trackParameters()->momentum().mag();
                eta = m->trackParameters()->momentum().eta();
                ATH_MSG_DEBUG("Start pars found eta " << eta << " r " << (m->trackParameters())->position().perp() << " z "
                                                      << (m->trackParameters())->position().z() << " pstart " << pstart);
            }
            if (m->materialEffectsOnTrack()) {
                const Trk::MaterialEffectsOnTrack* meot = dynamic_cast<const Trk::MaterialEffectsOnTrack*>(m->materialEffectsOnTrack());
                if (meot) {
                    if (meot->thicknessInX0() > 20) {
                        const Trk::ScatteringAngles* scatAngles = meot->scatteringAngles();
                        ATH_MSG_DEBUG(" Calorimeter X0  " << meot->thicknessInX0() << "  pointer scat " << scatAngles);
                        if (scatAngles) {
                            pcalo = m->trackParameters()->momentum().mag();
                            const double pullPhi = scatAngles->deltaPhi() / scatAngles->sigmaDeltaPhi();
                            const double pullTheta = scatAngles->deltaTheta() / scatAngles->sigmaDeltaTheta();
                            ATH_MSG_DEBUG(" Calorimeter scatterer deltaPhi " << scatAngles->deltaPhi() << " pull " << pullPhi
                                                                             << " deltaTheta " << scatAngles->deltaTheta() << " pull "
                                                                             << pullTheta);
                        }
                    }
                    const Trk::EnergyLoss* energyLoss = meot->energyLoss();
                    if (energyLoss) {
                        ATH_MSG_DEBUG("Eloss found r " << (m->trackParameters())->position().perp() << " z "
                                                       << (m->trackParameters())->position().z() << " value " << energyLoss->deltaE()
                                                       << " Eloss " << Eloss);
                        if (m->type(Trk::TrackStateOnSurface::CaloDeposit)) {
                            idEloss = Eloss;
                            caloEloss = std::abs(energyLoss->deltaE());
                            Eloss = 0.;
                            deltaP = m->trackParameters()->momentum().mag() - pcalo;
                            const Trk::Surface& surface = m->surface();
                            ATH_MSG_DEBUG(" Calorimeter surface " << surface);
                            ATH_MSG_DEBUG(txt << " Calorimeter delta p " << deltaP << " deltaE " << caloEloss
                                              << " delta pID = pcaloEntry-pstart " << pcalo - pstart);
                        } else {
                            Eloss += std::abs(energyLoss->deltaE());
                        }
                    }
                }
            }
        }
        msEloss = Eloss;
        Eloss = idEloss + caloEloss + msEloss;
        ATH_MSG_DEBUG(txt << " eta " << eta << " pstart " << pstart / 1000. << " Eloss on TSOS idEloss " << idEloss << " caloEloss "
                          << caloEloss << " msEloss " << msEloss << " Total " << Eloss << " pstart - pMuonEntry " << pstart - pMuonEntry);
    }

    bool MuonCombinedFitTagTool::extrapolatedNeedsRefit(const EventContext& ctx,const Trk::Track& combTrack, const Trk::Track* extrTrack) const {
        // will refit if extrapolated track was definitely bad
        if (!extrTrack) return true;
        if (!m_trackQuery->isCaloAssociated(*extrTrack, ctx)) return true;

        // otherwise will keep original SA fit if no change to MS or Calo TSOS
        const Trk::Track& originalTrack = *extrTrack;

        // refit if bad extrapolated fit - otherwise no refit if bad combined fit
        if (chi2(originalTrack.fitQuality()) > m_badFitChi2) return true;

        if (chi2(combTrack.fitQuality()) > m_badFitChi2) return true;

        // check if need to update calo association
        const CaloEnergy* caloEnergyCombined = m_trackQuery->caloEnergy(combTrack);
        const CaloEnergy* caloEnergyExtrapolated = m_trackQuery->caloEnergy(originalTrack);
        if (!caloEnergyCombined || !caloEnergyExtrapolated) {
            // no refit for combined track without CaloEnergy
            ATH_MSG_VERBOSE("extrapolatedNeedsRefit: no refit for combined track without CaloEnergy");
            return false;
        }
        double deltaE = caloEnergyExtrapolated->deltaE() - caloEnergyCombined->deltaE();
        if (std::abs(deltaE) > 0.3 * caloEnergyExtrapolated->sigmaDeltaE()) {
            ATH_MSG_VERBOSE("extrapolatedNeedsRefit: caloEnergy difference " << deltaE << "  sigma "
                                                                             << caloEnergyExtrapolated->sigmaDeltaE() << "  ratio "
                                                                             << deltaE / caloEnergyExtrapolated->sigmaDeltaE());
            return true;
        }

        Trk::TrackStates::const_reverse_iterator o = originalTrack.trackStateOnSurfaces()->rbegin();
        Trk::TrackStates::const_reverse_iterator c = combTrack.trackStateOnSurfaces()->rbegin();
        for (; o != originalTrack.trackStateOnSurfaces()->rend(); ++o) {
            if (dynamic_cast<const Trk::PerigeeSurface*>(&(**o).surface())) break;

            // compare measurements
            if ((**o).measurementOnTrack() && (**o).trackParameters()) {
                // check measurements in phase
                while (c != combTrack.trackStateOnSurfaces()->rend() && (!(**c).measurementOnTrack() || !(**c).trackParameters())) ++c;

                if (c == combTrack.trackStateOnSurfaces()->rend()) continue;

                double separation =
                    ((**o).trackParameters()->associatedSurface().center() - (**c).trackParameters()->associatedSurface().center()).mag();
                if (std::abs(separation) > 1. * CLHEP::mm) {
                    ATH_MSG_VERBOSE("extrapolatedNeedsRefit: measurement out-of-phase: "
                                    << " separation " << separation << "  extrap " << (**o).trackParameters()->associatedSurface().center()
                                    << "   comb " << (**c).trackParameters()->associatedSurface().center());
                    return true;
                }

                // different outlier
                if ((**o).type(Trk::TrackStateOnSurface::Outlier) != (**c).type(Trk::TrackStateOnSurface::Outlier)) {
                    if ((**c).type(Trk::TrackStateOnSurface::Outlier)) {
                        ATH_MSG_VERBOSE("extrapolatedNeedsRefit: outlier only on combined track ");
                    } else {
                        ATH_MSG_VERBOSE("extrapolatedNeedsRefit: outlier only on extrapolated track ");
                    }
                    return true;
                }

                // drift sign flip
                if (dynamic_cast<const Muon::MdtDriftCircleOnTrack*>((**o).measurementOnTrack())) {
                    if ((**o).measurementOnTrack()->localParameters()[Trk::driftRadius] *
                            (**c).measurementOnTrack()->localParameters()[Trk::driftRadius] <
                        0.) {
                        ATH_MSG_VERBOSE("extrapolatedNeedsRefit: drift sign flip ");
                        return true;
                    }
                }
                ++c;
            }
        }
        return false;
    }

    bool MuonCombinedFitTagTool::bestMatchChooser(const InDetCandidate& curCandidate, const CombinedFitTag& curTag,
                                                  const Trk::Track& curTrack, const Trk::Track* curMETrack,
                                                  const InDetCandidate& /*bestCandidate*/, const CombinedFitTag& bestTag,
                                                  const Trk::Track& bestTrack, const Trk::Track* bestMETrack) const

    {
        // pointers to extrapolated track
        const Trk::Track* curExtrTrack = curMETrack ? curMETrack : curTag.muonCandidate().extrapolatedTrack();
        const Trk::Track* bestExtrTrack = bestMETrack ? bestMETrack : bestTag.muonCandidate().extrapolatedTrack();

        // 1 current
        // 2 best
        // returned bool: true means current is better; false means "best is better"

        const double matchChiSq1 = curTag.matchChi2();
        const double matchChiSq2 = bestTag.matchChi2();
        const Trk::TrackSummary* summary1 = curTrack.trackSummary();
        const Trk::TrackSummary* summary2 = bestTrack.trackSummary();
        ATH_MSG_VERBOSE("bestMatchChooser: matchChiSq " << matchChiSq1 << "  " << matchChiSq2);
        if (summary1 && summary2) {
            ATH_MSG_VERBOSE("bestMatchChooser: matchChiSq " << matchChiSq1 << "  " << matchChiSq2 << "  numTRTHits "
                                                            << summary1->get(Trk::numberOfTRTHits) << "  "
                                                            << summary2->get(Trk::numberOfTRTHits) << "  field integrals: ID  "
                                                            << curTag.fieldIntegral().betweenInDetMeasurements() << "  "
                                                            << bestTag.fieldIntegral().betweenInDetMeasurements() << "  MS "
                                                            << curTag.fieldIntegral().betweenSpectrometerMeasurements() << "  "
                                                            << bestTag.fieldIntegral().betweenSpectrometerMeasurements());
        } else {
            ATH_MSG_VERBOSE("bestMatchChooser: matchChiSq " << matchChiSq1 << "  " << matchChiSq2);
        }

        // selection when only one match has a good combined fit
        const double fitChiSq1 = chi2(curTrack.fitQuality());
        const double fitChiSq2 = chi2(bestTrack.fitQuality());      
        const unsigned int numberDoF1 = curTrack.fitQuality()->numberDoF();
        const unsigned int numberDoF2 = bestTrack.fitQuality()->numberDoF();
        ATH_MSG_VERBOSE("bestMatchChooser: fitChiSq " << fitChiSq1 << "  " << fitChiSq2);
        if (std::abs(fitChiSq1 - fitChiSq2) > m_badFitChi2) {
            if (fitChiSq1 < m_badFitChi2) {
                if (matchChiSq1 > matchChiSq2 && matchChiSq2 < m_matchChiSquaredCut) {  // may want to suppress this warning!
                    ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq as better fitChiSq. "
                                    << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2 << "   fitChiSq/DoF 1,2 " << fitChiSq1
                                    << "/" << numberDoF1 << ", " << fitChiSq2 << "/" << numberDoF2);
                }
                return true;
            }
            if (fitChiSq2 < m_badFitChi2) {
                if (matchChiSq1 < matchChiSq2 && matchChiSq1 < m_matchChiSquaredCut) {  // may want to suppress this warning!
                    ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq as better fitChiSq. "
                                    << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2 << "   fitChiSq/DoF 1,2 " << fitChiSq1
                                    << "/" << numberDoF1 << ", " << fitChiSq2 << "/" << numberDoF2);
                }
                return false;
            }
        }

        // selection when only one match has a good match chi2
        if (std::abs(matchChiSq1 - matchChiSq2) > m_matchChiSquaredCut) {
            if (matchChiSq1 < m_matchChiSquaredCut) return true;
            if (matchChiSq2 < m_matchChiSquaredCut) return false;
        }

        // energy balance (absolute)
        // track length:
        //          field integral
        // 		# MS stations
        // 		pixel hits (-outliers)
        // 		trt drift hits + outliers

        // protect momentum balance and field integral when magnets off:
        if (!curCandidate.indetTrackParticle().track()->info().trackProperties(Trk::TrackInfo::StraightTrack)) {
            double cutRatio{1.5}, integral1{0.}, integral2{0.};

            if (curExtrTrack && !curExtrTrack->info().trackProperties(Trk::TrackInfo::StraightTrack) && bestExtrTrack &&
                !bestExtrTrack->info().trackProperties(Trk::TrackInfo::StraightTrack)) {
                // selection when only one match has good momentum balance or a significantly better balance
                ATH_MSG_VERBOSE("bestMatchChooser: momentumBalanceSignificance " << curTag.momentumBalanceSignificance() << "  "
                                                                                 << bestTag.momentumBalanceSignificance());
                double significanceCut = 2.0;
                const double significance1 = std::abs(curTag.momentumBalanceSignificance());
                const double significance2 = std::abs(bestTag.momentumBalanceSignificance());
                if (std::abs(significance1 - significance2) > significanceCut) {
                    if (significance1 < significanceCut) {
                        if (matchChiSq1 > matchChiSq2 && matchChiSq2 < m_matchChiSquaredCut) {
                            // NOT choosing bestMatchChi2:
                            ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq as better momentum balance. "
                                            << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2
                                            << "   momentumBalanceSignificance 1,2 " << curTag.momentumBalanceSignificance() << ", "
                                            << bestTag.momentumBalanceSignificance());
                        }
                        return true;
                    }
                    if (significance2 < significanceCut) {
                        if (matchChiSq1 < matchChiSq2 && matchChiSq1 < m_matchChiSquaredCut) {
                            // NOT choosing bestMatchChi2:
                            ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq as better momentum balance. "
                                            << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2
                                            << "   momentumBalanceSignificance 1,2 " << curTag.momentumBalanceSignificance() << ", "
                                            << bestTag.momentumBalanceSignificance());
                        }
                        return false;
                    }
                }

                // keep significantly larger measured field integral
                //   for MS
                ATH_MSG_VERBOSE("bestMatchChooser: spectrometer field integral ratio ");
                integral1 = std::abs(curTag.fieldIntegral().betweenSpectrometerMeasurements());
                integral2 = std::abs(bestTag.fieldIntegral().betweenSpectrometerMeasurements());
                if (integral1 > cutRatio * integral2) return true;
                if (integral2 > cutRatio * integral1) return false;
            }
            //   for indet
            ATH_MSG_VERBOSE("bestMatchChooser: indet field integral ratio ");
            integral1 = std::abs(curTag.fieldIntegral().betweenInDetMeasurements());
            integral2 = std::abs(bestTag.fieldIntegral().betweenInDetMeasurements());
            if (integral1 > cutRatio * integral2) return true;
            if (integral2 > cutRatio * integral1) return false;
        }

        // repeat fit/match quality selection with sharper cuts (times 2)
        ATH_MSG_VERBOSE("bestMatchChooser: sharp fit chi2 cut ");
        if (std::abs(fitChiSq1 - fitChiSq2) > 0.5 * m_badFitChi2) {
            if (fitChiSq1 < 0.5 * m_badFitChi2) {
                if (matchChiSq1 > matchChiSq2 && matchChiSq2 < m_matchChiSquaredCut) {
                    // NOT choosing bestMatchChi2:
                    ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq according to overall quality. "
                                    << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2 << "   fitChiSq/DoF 1,2 " << fitChiSq1
                                    << "/" << numberDoF1 << ", " << fitChiSq2 << "/" << numberDoF2);
                }
                return true;
            }
            if (fitChiSq2 < 0.5 * m_badFitChi2) {
                if (matchChiSq1 < matchChiSq2 && matchChiSq1 < m_matchChiSquaredCut) {
                    // NOT choosing bestMatchChi2:
                    ATH_MSG_WARNING("bestMatchChooser: choose worse, but acceptable, matchChiSq according to overall quality. "
                                    << " matchChiSq 1,2 " << matchChiSq1 << ", " << matchChiSq2 << "   fitChiSq/DoF 1,2 " << fitChiSq1
                                    << "/" << numberDoF1 << ", " << fitChiSq2 << "/" << numberDoF2);
                }
                return false;
            }
        }

        ATH_MSG_VERBOSE("bestMatchChooser: sharp match chi2 cut ");
        if (std::abs(matchChiSq1 - matchChiSq2) > 0.5 * m_matchChiSquaredCut) {
            if (matchChiSq1 < 0.5 * m_matchChiSquaredCut) return true;
            if (matchChiSq2 < 0.5 * m_matchChiSquaredCut) return false;
        }

        // track quality:
        // pixel holes (outliers)
        // silicon holes + outliers
        // trt drift hits + outliers

        // kink finding:
        // neighbour signif
        // curvature signif

        // energy balance (absolute)

        // field off protection
        if (curExtrTrack && curExtrTrack->info().trackProperties(Trk::TrackInfo::StraightTrack) && bestExtrTrack &&
            bestExtrTrack->info().trackProperties(Trk::TrackInfo::StraightTrack)) {
            // best fit chi2
            return fitChiSq1 < fitChiSq2;
        } else {
            // best match chi2
            return matchChiSq1 < matchChiSq2;
        }
    }
}  // namespace MuonCombined
