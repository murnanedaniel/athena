/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MCAST_MUONCALIBRATIONANDMEARINGTOOL_H
#define MCAST_MUONCALIBRATIONANDMEARINGTOOL_H

// Framework include(s):
#include "AsgDataHandles/ReadHandleKey.h"
#include "AsgTools/AnaToolHandle.h"
#include "AsgTools/AsgTool.h"
#include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "PATInterfaces/SystematicsCache.h"
#include "xAODEventInfo/EventInfo.h"

// ROOT include(s)
#include "TFile.h"
#include "TH3F.h"
#include "TMatrix.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include "TVectorD.h"

// C++ include(s)
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#define EPSILON 1.0E-6
#define DEFAULT_INIT_VAL -999
#define MCAST_MAX_PT 100000000
#define MZPDG 91.1876
namespace CP {

    namespace MCAST {

        namespace DataType {
            enum { Data10 = 1, Data11 = 2, Data12 = 3, Data15 = 4, Data16 = 5, Data17 = 6, Data18 = 7 };
        }
        namespace AlgoType {
            enum { Muid = 1, Staco = 2, Muons = 3 };
        }
        namespace Release {
            enum {
                Rel16_6 = 1,
                Rel17 = 2,
                Rel17_2 = 3,
                Rel17_2_Repro = 4,
                Rel17_2_Sum13 = 5,
                PreRec = 6,
                PreRec_2015_06_22 = 7,
                PreRec_2015_08_06 = 8,
                Rec_2015_11_15 = 9,
                Rec_2016_01_13 = 10,
                Rec_2016_01_19 = 11,
                PreRec_2016_05_23 = 12,
                Recs2016_08_07 = 13,
                Recs2016_15_07 = 14,
                Recs2017_08_02 = 15,
                Recs2019_05_30 = 16,
                Recs2019_10_12 = 17,
                Recs2020_03_03 = 18,
                Recs2021_07_01 = 19
            };
        }
        namespace SmearingType {
            enum { Pt = 1, QoverPt = 2 };
        }
        namespace DetectorType {
            enum { MS = 1, ID = 2, CB = 3 };
        }
        namespace SystVariation {
            enum { Default = 0, Down = -1, Up = 1 };
        }
        namespace SagittaCorType {
            enum { CB = 0, ID = 1, ME = 2, WEIGHTS = 3, AUTO = 4 };
        }
        namespace SagittaSysType {
            enum { NOMINAL = 0, RHO = 1, BIAS = 2, DATASTAT = 3 };
        }
        namespace MST_Categories {
            enum { Undefined = -1, Zero = 0, One = 1, Two = 2, Three = 3, Four = 4, Total = 5 };
        }
        namespace SagittaInputHistType {
            enum { NOMINAL = 0, SINGLE = 1 };
        }
    }  // namespace MCAST
    //  [[deprecated]] 
    class MuonCalibrationAndSmearingTool : public virtual IMuonCalibrationAndSmearingTool,
                                           public virtual ISystematicsTool,
                                           public asg::AsgTool {
        // Create a proper constructor for Athena
        ASG_TOOL_CLASS3(MuonCalibrationAndSmearingTool, CP::IMuonCalibrationAndSmearingTool, CP::ISystematicsTool,
                        CP::IReentrantSystematicsTool)

    public:
        // Interface methods that must be defined
        // Interface - Apply the correction on a modifyable object
        virtual CorrectionCode applyCorrection(xAOD::Muon& mu) const override;
        // Interface - Create a corrected copy from a constant muon
        virtual CorrectionCode correctedCopy(const xAOD::Muon& input, xAOD::Muon*& output) const override;
        // Interface - Is the tool affected by a specific systematic?
        virtual bool isAffectedBySystematic(const SystematicVariation& systematic) const override;
        // Interface - Which systematics have an effect on the tool's behaviour?
        virtual SystematicSet affectingSystematics() const override;
        // Interface - Systematics to be used for physics analysis
        virtual SystematicSet recommendedSystematics() const override;
        // Interface - Use specific systematic
        virtual StatusCode applySystematicVariation(const SystematicSet& systConfig) override;
        // Interface - get the expected resolution of the muon
        virtual double expectedResolution(const std::string& DetType, const xAOD::Muon& mu, const bool mc) const override;
        // Interface - get the expected resolution of the muon
        virtual double expectedResolution(const int& DetType, const xAOD::Muon& mu, const bool mc) const override;

    public:
        // InfoHelper is intended to be used to ease the passing of information between internal
        // methods within this class.  This is created in anticipation of usage in AthenaMT
        struct InfoHelper {
            double ptms = 0;
            double ptid = 0;
            double ptcb = 0;
            double eta = 0;
            double phi = 0;
            double sagitta_calibrated_ptcb = 0;
            double sagitta_calibrated_ptid = 0;
            double sagitta_calibrated_ptms = 0;
            double g0;
            double g1;
            double g2;
            double g3;
            double g4;
            double extra_g;
            int charge = 1;
            int detRegion = 0;
            std::vector<float> cbParsA;
            std::vector<float> cbCovMat;
            double weightMS = 0;
            double weightID = 0;
            double smearDeltaMS = 0;
            double smearDeltaID = 0;
            double smearDeltaCB = 0;
            double smearDeltaCBOnly = 0;
            double smearDeltaCBDirect = 0;
            int sel_category = -1;
            double uncorrected_ptcb = 0;
            double uncorrected_ptid = 0;
            double uncorrected_ptms = 0;
        };

    public:
        // Constructor
        MuonCalibrationAndSmearingTool(const std::string& name);

        // Copy constructor
        MuonCalibrationAndSmearingTool(const MuonCalibrationAndSmearingTool& tool);

        // Destructor
        virtual ~MuonCalibrationAndSmearingTool();

        virtual StatusCode initialize() override;

        double ExpectedResolution(const std::string& DetType, const xAOD::Muon& mu, const bool mc) const;
        double ExpectedResolution(const int DetType, const xAOD::Muon& mu, const bool mc) const;

        // Expert method to apply the MC correction on a modifyable trackParticle for ID- or MS-only corrections
        virtual CorrectionCode applyCorrectionTrkOnly(xAOD::TrackParticle& inTrk, const int DetType) const override;

        virtual CorrectionCode applyStatCombination(AmgVector(5) parsID, AmgSymMatrix(5) covID, AmgVector(5) parsMS, AmgSymMatrix(5) covMS,
                                                    int charge, AmgVector(5) & parsCB, AmgSymMatrix(5) & covCB, double& chi2) const;
        virtual CorrectionCode applyStatCombination(xAOD::Muon& mu, InfoHelper& muonInfo) const;
        virtual CorrectionCode applySagittaBiasCorrectionAuto(const int DetType, xAOD::Muon& mu, bool isMC, const unsigned int SystCase,
                                                              InfoHelper& muonInfo) const;
        virtual CorrectionCode CorrectForCharge(double p2, double& pt, int q, bool isMC, double p2Kin = 0) const;
        virtual CorrectionCode applySagittaBiasCorrection(const unsigned int SgCorrType, xAOD::Muon& mu, unsigned int iter, bool stop,
                                                          bool isMC, InfoHelper& muonInfo, const unsigned int SystCase = 0) const;

    protected:
        // Regions helpers
        StatusCode Regions(std::string inRegionFile, int doMacroRegionsFlag = 0);
        void PrintRegions() const;
        unsigned int GetNRegions() const;
        int GetRegion(const double eta, const double phi) const;
        float GetRegionInnerEta(const int r_i) const;  // Return Eta closer to the origin
        std::string GetRegionName(const int r_i) const;
        std::string GetRegionName(const double eta, const double phi) const;
        double GetSmearing(int DetType, InfoHelper& muonInfo, bool doDirectCB) const;
        double GetSystVariation(int DetType, double var, InfoHelper& muonInfo, bool doDirectCB) const;
        StatusCode SetInfoHelperCorConsts(InfoHelper& inMuonInfo) const;
        void CalcCBWeights(xAOD::Muon&, InfoHelper& muonInfo) const;
        double CalculatePt(const int DetType, const double inSmearID, const double inSmearMS, const double scaleVarID,
                           const double scaleMS_scale, const double scaleMS_egLoss, const double scaleCB_scale, const double scaleCB_egLoss,
                           InfoHelper& muonInfo) const;
        StatusCode FillValues();
        void Clean();
        double ScaleApply(const double pt, double S, const double S_EnLoss, InfoHelper& muonInfo) const;
        void CleanScales();
        void CollectMacroRegionsSL();           // Small and large regions are collected together
        void CollectMacroRegionsSL_UpDn();      // Small,Large,Up,Down regions are collected together
        void CollectMacroRegionsSL_SplitBAR();  // Large,Small sectors split plus Feet(12+14) and 11+15 sector split in Barrel
        void CollectSectors();

        StatusCode SetData(std::string);
        StatusCode SetAlgorithm(std::string);
        StatusCode SetRelease(std::string);
        StatusCode SetType(std::string);

        virtual unsigned int setSagittaHistogramsSingle(TProfile2D* pCB = nullptr, unsigned int track = 0);
        virtual double sagitta(TProfile2D* corrM, TLorentzVector& lv) const;

        virtual void ConvertToSagittaBias(TH2F* h, float mean = 1);
        virtual TProfile2D* GetHist(std::string fname = "", std::string hname = "inclusive", double GlobalScale = MZPDG);
        virtual TProfile2D* GetHistSingleMethod(std::string fname = "", std::string hname = "");
        virtual bool isBadMuon(const xAOD::Muon& mu, InfoHelper& muonInfo) const;
        int ConvertToMacroCategory(const int raw_mst_category) const;
        // private:
        // fake assignment operator missing actual implementation
        MuonCalibrationAndSmearingTool& operator=(const MuonCalibrationAndSmearingTool&);

        struct ParameterSet {
            double SmearTypeID;
            double SmearTypeMS;
            double SmearTypeCB;
            double ScaleID;
            double ScaleMS_scale;
            double ScaleMS_egLoss;
            double ScaleCB_scale;
            double ScaleCB_egLoss;
            double SagittaRho;
            double SagittaBias;
            double SagittaDataStat;
            double SagittaEtaSlice;
        };

        // calculate the parameter set for the given systematic
        StatusCode calcSystematicVariation(const SystematicSet& systConfig, ParameterSet& param) const;

        SG::ReadHandleKey<xAOD::EventInfo> m_eventInfo{this, "EventInfoContName", "EventInfo", "event info key"};

        bool m_expertMode = false;
        bool m_expertMode_isData = false;

        bool m_useExternalSeed;
        int m_externalSeed;

        std::string m_year, m_algo, m_type, m_release;
        std::string m_FilesPath, m_sysScheme;
        bool m_extra_highpt_smearing;
        bool m_2stations_highpt_smearing;
        bool m_extra_decorations;
        bool m_toroidOff;
        int m_Tsmear;
        int m_Tdata;
        int m_Trel;
        int m_Talgo;
        double m_extraRebiasSys = 0.0;
        double m_useNsigmaForICombine;
        bool m_doDirectCBCalib = false;
        bool m_doEtaSagittaSys = false;
        std::vector<double> m_scale_ID, m_enLoss_MS, m_scale_MS, m_scale_CB;

        // sys variations (stat error added in quadrature), one if it's simmetrized, 2 if Up != Dw.
        std::vector<double> m_scaleSyst_ID, m_enLossSyst_MS, m_scaleSyst_MS, m_scaleSyst_CB;
        std::vector<double> m_scaleSystUp_ID, m_enLossSystUp_MS, m_scaleSystUp_MS;
        std::vector<double> m_scaleSystDw_ID, m_enLossSystDw_MS, m_scaleSystDw_MS;

        std::vector<double> m_scaleBins;
        std::vector<double> m_p1_ID, m_p2_ID, m_p2_ID_TAN, m_p0_MS, m_p1_MS, m_p2_MS;
        std::vector<double> m_E_p1_ID, m_E_p2_ID, m_E_p2_ID_TAN, m_E_p0_MS, m_E_p1_MS, m_E_p2_MS;
        // syst. errors on resolution parameters corrections:
        // one if it's simmetrized, then Stat and Sys err are separate in cfg file.
        std::vector<double> m_S_p1_ID, m_S_p2_ID, m_S_p2_ID_TAN, m_S_p0_MS, m_S_p1_MS, m_S_p2_MS;
        // Two if Up != Dw, Stat and Sys err added in quadrature in cfg file.
        std::vector<double> m_SUp_p1_ID, m_SUp_p2_ID, m_SUp_p2_ID_TAN, m_SUp_p0_MS, m_SUp_p1_MS, m_SUp_p2_MS;
        std::vector<double> m_SDw_p1_ID, m_SDw_p2_ID, m_SDw_p2_ID_TAN, m_SDw_p0_MS, m_SDw_p1_MS, m_SDw_p2_MS;
        std::vector<double> m_MC_p1_ID, m_MC_p2_ID, m_MC_p2_ID_TAN, m_MC_p0_MS, m_MC_p1_MS, m_MC_p2_MS;

        std::vector<double> m_S_0_CB, m_SUp_0_CB, m_SDw_0_CB;
        std::vector<double> m_S_1_CB, m_SUp_1_CB, m_SDw_1_CB;
        std::vector<double> m_R_0_CB, m_RUp_0_CB, m_RDw_0_CB;
        std::vector<double> m_R_1_CB, m_RUp_1_CB, m_RDw_1_CB;
        std::vector<double> m_R_2_CB, m_RUp_2_CB, m_RDw_2_CB;
        
        // Special "p2" systematics and corrections for non-three-station muons
        // Maps have two keys: detector region and category
        std::map<std::pair<int, int>, std::pair<double, double> > m_extra_p1_p2_MS_AlignedOnly, m_extra_p1_p2_MS_AlignedAndCorrected,
            m_extra_p1_p2_MS_Misaligned;

        std::vector<std::string> m_names;
        bool m_loadNames;
        int m_nb_regions;
        std::vector<float> m_eta_min, m_eta_max, m_phi_min, m_phi_max;

        bool m_doMacroRegions;
        std::map<int, int> m_MacroRegionIdxMap;
        std::vector<std::string> m_MacroRegionName;
        std::vector<double> m_MacroRegionInnerEta;

        SystematicsCache<ParameterSet> m_Parameters{this};
        const ParameterSet* m_currentParameters{nullptr};

        double m_StatCombPtThreshold;
        double m_HighPtSystThreshold;
        bool m_useStatComb;

        unsigned int m_sgItersID = 0U;
        unsigned int m_sgItersCB = 0U;
        unsigned int m_sgItersME = 0U;
        bool m_sgIetrsManual = false;
        double m_fixedRho = 0.0;
        bool m_useFixedRho = false;
        double m_sagittaMapUnitConversion;

        std::vector<std::unique_ptr<TProfile2D> > m_sagittasCB;
        std::vector<std::unique_ptr<TProfile2D> > m_sagittasID;
        std::vector<std::unique_ptr<TProfile2D> > m_sagittasME;

        bool m_SagittaCorrPhaseSpace;
        bool m_doSagittaCorrection;
        bool m_doSagittaMCDistortion = false;
        bool m_doNotUseAMGMATRIXDECOR;
        float m_IterWeight;

        std::unique_ptr<TProfile2D> m_sagittaPhaseSpaceCB;
        std::unique_ptr<TProfile2D> m_sagittaPhaseSpaceID;
        std::unique_ptr<TProfile2D> m_sagittaPhaseSpaceME;

        std::string m_SagittaRelease;
        std::vector<unsigned int> m_SagittaIterations;
        std::vector<double> m_GlobalZScales;
        unsigned int m_saggitaMapsInputType;

        asg::AnaToolHandle<CP::IMuonSelectionTool> m_MuonSelectionTool;

    };  // class MuonCalibrationAndSmearingTool

}  // namespace CP

#endif
