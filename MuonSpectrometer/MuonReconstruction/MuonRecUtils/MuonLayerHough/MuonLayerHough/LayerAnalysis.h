/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LAYERANALYSIS_H
#define LAYERANALYSIS_H

#include "CxxUtils/checker_macros.h"
#include "MuonLayerHough/HitNtuple.h"
#include "MuonLayerHough/MuonLayerHough.h"
#include "MuonLayerHough/MuonLayerHoughSelector.h"
#include "MuonLayerHough/MuonRegionHough.h"
#include "MuonStationIndex/MuonStationIndex.h"
#include <vector>
#include <map>

class TH1F;
class TH2F;
class TTree;

namespace MuonHough {

    class MuonDetectorHough;

    struct Plots {
        Plots(const char* title, int nBinsX, float xMin, float xMax, int nBinsY, float yMin, float yMax);

        TH2F* Reco{};
        TH2F* Truth{};
        TH2F* Matched{};
        TH2F* Unmatched{};
        TH1F* Efficiency{};
        TH1F* FakeEfficiency{};
        TH2F* Diff{};
    };

    class LayerAnalysis {
    public:
        LayerAnalysis(TTree& tree) : m_tree(&tree), m_ncalls(0) {
            m_ntuple.reset();
            m_ntuple.initForRead(tree);
        }

        void initialize ATLAS_NOT_THREAD_SAFE();
        void analyse ATLAS_NOT_THREAD_SAFE();
        void finalize();

    private:
        void analysis(std::map<int, SectorData>& data);
        void drawSector(int region, int sector, SectorData& data, MuonDetectorHough& detectorHough,
                        MuonDetectorHough& detectorHoughTruth) const;
        static void calculateVariables(Plots* Plot);
        static void SetStyle ATLAS_NOT_THREAD_SAFE();
        static int getMaxLayers(Muon::MuonStationIndex::DetectorRegionIndex region, int sector);
        void finishplot(TH1F* h) const;
        static float linear_extrapolate(const MuonLayerHough::Maximum& ref, const MuonLayerHough::Maximum& ex);
        float parab_extrapolate(const MuonLayerHough::Maximum& ref, const MuonLayerHough::Maximum& ex) const;

        TTree* m_tree{};
        HitNtuple m_ntuple{};
        EventData m_event{};
        int m_ncalls{};
        std::vector<MuonDebugInfo> m_muons;
        std::vector<SegDebugInfo> m_segs;
        std::vector<MuonHough::MuonLayerHoughSelector> m_selectors;

        std::vector<Plots*> m_hMaximaHeightPerChIndex;
        bool m_DEBUG{};
        bool m_DEBUG_seg{};
        TH1F* m_h_dtheta{};
        TH1F* m_h_dtheta_truth{};
        TH1F* m_h_diff_MI_e{};
        TH1F* m_h_diff_MO_e{};
        TH1F* m_h_diff_MI_b{};
        TH1F* m_h_diff_MO_b{};
        TH1F* m_h_diff_MI_e_truth{};
        TH1F* m_h_diff_MO_e_truth{};
        TH1F* m_h_diff_MI_b_truth{};
        TH1F* m_h_diff_MO_b_truth{};
        TH1F* m_h_expt_MI_e{};
        TH1F* m_h_expt_MO_e{};
        TH1F* m_h_expt_MI_b{};
        TH1F* m_h_expt_MO_b{};
        TH1F* m_h_expt_MI_e_truth{};
        TH1F* m_h_expt_MO_e_truth{};
        TH1F* m_h_expt_MI_b_truth{};
        TH1F* m_h_expt_MO_b_truth{};
        TH1F* m_h_comp_MI_e{};
        TH1F* m_h_comp_MO_e{};
        TH1F* m_h_comp_MI_b{};
        TH1F* m_h_comp_MO_b{};
        TH1F* m_h_comp_MI_e_truth{};
        TH1F* m_h_comp_MO_e_truth{};
        TH1F* m_h_comp_MI_b_truth{};
        TH1F* m_h_comp_MO_b_truth{};

        MuonHough::MuonDetectorHough m_detectorHough{};
        MuonHough::MuonDetectorHough m_detectorHoughTruth{};
    };

}  // namespace MuonHough

#endif
