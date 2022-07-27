/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
//   get quality of combined track match
///////////////////////////////////////////////////////////////////

#ifndef MUIDTRACKBUILDER_MUONMATCHQUALITY_H
#define MUIDTRACKBUILDER_MUONMATCHQUALITY_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaKernel/Units.h"
#include "GaudiKernel/ToolHandle.h"
#include "MuidInterfaces/IMuonMatchQuality.h"
#include "MuidInterfaces/IMuonTrackQuery.h"
#include "MuonCombinedToolInterfaces/IMuonTrackTagTool.h"

namespace Units = Athena::Units;

namespace Rec {

    class MuonMatchQuality : public AthAlgTool, virtual public IMuonMatchQuality {
    public:
        MuonMatchQuality(const std::string& type, const std::string& name, const IInterface* parent);
        virtual ~MuonMatchQuality() = default;  // destructor

        StatusCode initialize() override;

        /** IMuonMatchQuality interface:
            match chiSquared between two tracks expressed at same inner (IP) surface,
            expected to handle indet with extrapolated spectrometer track or combined
            with constituent track */
        double innerMatchChi2(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const override;

        /** IMuonMatchQuality interface:
            degrees of freedom for chi2 match at IP */
        int innerMatchDOF(const Trk::Track& track1, const Trk::Track& track2) const override;

        /** IMuonMatchQuality interface:
            match probability for chi2 match at IP */
        double innerMatchProbability(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const override;

        /** IMuonMatchQuality interface:
            degrees of freedom, chi2, probability  for chi2 match at IP */
        std::pair<int, std::pair<double, double> > innerMatchAll(const Trk::Track& track1, const Trk::Track& track2,
                                                                 const EventContext& ctx) const override;

        /** IMuonMatchQuality interface:
            match chiSquared between two tracks expressed at first muon spectrometer hit,
            extrapolates indet to first hit of spectrometer track */
        double outerMatchChi2(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const override;

        /** IMuonMatchQuality interface:
            degrees of freedom for chi2 match at first MS hit */
        virtual int outerMatchDOF(const Trk::Track& track1, const Trk::Track& track2) const override;

        /** IMuonMatchQuality interface:
            match probability for chi2 match at first MS hit */
        double outerMatchProbability(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const override;

        /** IMuonMatchQuality interface:
            check the track perigee parameters are expressed at the same surface */
        virtual bool shareOrigin(const Trk::Track& track1, const Trk::Track& track2) const override;

        /** IMuonMatchQuality interface:
            as inner match chiSquared but simplified to just use diagonal errors */
        double simpleChi2(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const override;

    private:
        // cache
        struct CacheAll {
            double innerMatchChi2{FLT_MAX};
            int innerMatchDOF{0};
            double innerMatchProbability{0.};
            double simpleChi2{FLT_MAX};
        };

        CacheAll setCache(const Trk::Track& track1, const Trk::Track& track2, const EventContext& ctx) const;

        // helpers, managers, tools
        ToolHandle<MuonCombined::IMuonTrackTagTool> m_tagTool{
            this,
            "TagTool",
            "",
            "Track tag tool",
        };

        ToolHandle<IMuonTrackQuery> m_trackQuery{
            this,
            "TrackQuery",
            "Rec::MuonTrackQuery/MuonTrackQuery",
            "Track query tool",
        };

        // estimate of ID/MS alignment uncertainties
        std::unique_ptr<AmgSymMatrix(5)> m_alignmentUncertainty;
        Gaudi::Property<double> m_directionUncertainty{
            this,
            "ID_MS_DirectionUncertainty",
            0.000001,
        };  // not used anymore angle ID and MS: done by m_addIDMSerrors

        Gaudi::Property<double> m_positionUncertainty{
            this,
            "ID_MS_PositionUncertainty",
            0.01 * Units::mm,
        };  // not used anymore shift ID and MS: done by m_addIDMSerrors

    };  // end of class MuonMatchQuality

}  // end of namespace Rec

#endif  // MUIDTRACKBUILDER_MUONMATCHQUALITY_H
