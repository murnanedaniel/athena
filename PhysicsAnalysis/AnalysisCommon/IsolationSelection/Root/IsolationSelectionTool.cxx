/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
 */

#include <xAODPrimitives/IsolationType.h>
#include <IsolationSelection/IsolationSelectionTool.h>
#include <IsolationSelection/IsolationConditionHist.h>
#include <IsolationSelection/IsolationConditionFormula.h>
#include <IsolationSelection/IsolationConditionCombined.h>
#include <IsolationSelection/Interp3D.h>
#include "PathResolver/PathResolver.h"
#include <TROOT.h>
#include <TFile.h>
#include <TObjString.h>
#include <TH3.h>
#include <TF2.h>

namespace CP {
    IsolationSelectionTool::IsolationSelectionTool(const std::string& name) :
                asg::AsgTool(name),
                /// input file
                m_calibFile(nullptr),

                /// AcceptInfo's
                m_photonAccept("IsolationSelectionToolPhotonAcceptInfo"),
                m_electronAccept("IsolationSelectionToolElectronAcceptInfo"),
                m_muonAccept("IsolationSelectionToolMuonAcceptInfo"),
                m_objAccept("IsolationSelectionToolObjAcceptInfo"),
                m_iparWPs(0),
                m_Interp(nullptr),
                m_TwikiLoc("https://twiki.cern.ch/twiki/bin/view/AtlasProtected/IsolationSelectionTool#List_of_current_official_working") {
        declareProperty("CalibFileName", m_calibFileName = "", "The config file to use");
        declareProperty("MuonWP",       m_muWPname = "Undefined" , "Working point for muon");
 	declareProperty("ElectronWP",   m_elWPname = "Undefined" , "Working point for electron");
 	declareProperty("PhotonWP",     m_phWPname = "Undefined" , "Working point for photon");
	declareProperty("MuonWPVec",    m_muWPvec                , "Vector of working points for muon");
 	declareProperty("ElectronWPVec",m_elWPvec                , "Vector of working points for electron");
 	declareProperty("PhotonWPVec",  m_phWPvec                , "Vector of working points for photon");
	declareProperty("MuonKey",      m_muWPKey  = "/Muons/DFCommonGoodMuon/mu_cutValues_", "path of the cut map for muon");
	declareProperty("ElectronKey",  m_elWPKey  = "/ElectronPhoton/LHTight/el_cutValues_", "path of the cut map for electron");
	declareProperty("PhotonKey",    m_phWPKey  = "/ElectronPhoton/LHTight/el_cutValues_", "path of the cut map for photon");

        declareProperty("doCutInterpolationMuon", m_doInterpM = false, "flag to perform cut interpolation, muon");
        declareProperty("doCutInterpolationElec", m_doInterpE = true, "flag to perform cut interpolation, electron");
    }
    void IsolationSelectionTool::clearWPs(std::vector<IsolationWP*>& WP) {
        for (auto& c : WP) {
            if (c) delete c;
        }
        WP.clear();
    }

    void IsolationSelectionTool::clearPhotonWPs() {
        clearWPs(m_phWPs);
    }
    void IsolationSelectionTool::clearElectronWPs() {
        clearWPs(m_elWPs);
    }
    void IsolationSelectionTool::clearMuonWPs() {
        clearWPs(m_muWPs);
    }
    void IsolationSelectionTool::clearObjWPs() {
        clearWPs(m_objWPs);
    }
    const std::vector<IsolationWP*>& IsolationSelectionTool::getMuonWPs() const {
        return m_muWPs;
    }
    const std::vector<IsolationWP*>& IsolationSelectionTool::getElectronWPs() const {
        return m_elWPs;
    }
    const std::vector<IsolationWP*>& IsolationSelectionTool::getPhotonWPs() const {
        return m_phWPs;
    }
    const std::vector<IsolationWP*>& IsolationSelectionTool::getObjWPs() const {
        return m_objWPs;
    }
    IsolationSelectionTool::~IsolationSelectionTool() {
        clearMuonWPs();
        clearPhotonWPs();
        clearElectronWPs();
        clearObjWPs();
    }

    StatusCode IsolationSelectionTool::initialize() {
        /// Greet the user:
        ATH_MSG_INFO("Initialising...");
	if(m_calibFileName!="") {
	  std::string filename = PathResolverFindCalibFile(m_calibFileName);

	  ATH_MSG_INFO("Reading input file " << m_calibFileName << " from " << filename);
	  m_calibFile = std::shared_ptr < TFile > (new TFile(filename.c_str(), "READ"));

	  TObjString* versionInfo = (TObjString*) m_calibFile->Get("VersionInfo");
	  if (versionInfo) ATH_MSG_INFO("VersionInfo:" << versionInfo->String());
	  else ATH_MSG_WARNING("VersionInfo of input file (" << filename << ") is missing.");
	}

        if (m_doInterpE || m_doInterpM) {
            // special setting for electrons
            // do not apply interpolation in crack vicinity for topoetcone
            std::vector<std::pair<double, double>> rangeEtaNoInt;
            std::pair<double, double> apair(1.26, 1.665);
            rangeEtaNoInt.push_back(apair);
            // do not apply interpolation between Z defined and J/Psi defined cuts (pT < > 15 GeV/c) for both calo and track iso
            std::vector<std::pair<double, double>> rangePtNoInt;
            apair.first = 12.5;
            apair.second = 17.5;
            rangePtNoInt.push_back(apair);
            std::map<std::string, Interp3D::VetoInterp> amap;
            Interp3D::VetoInterp veto;
            veto.xRange = rangePtNoInt;
            veto.yRange = std::vector<std::pair<double, double>>();
            amap.insert(std::make_pair(std::string("el_cutValues_ptvarcone20"), veto));
            veto.yRange = rangeEtaNoInt;
            amap.insert(std::make_pair(std::string("el_cutValues_topoetcone20"), veto));
            m_Interp = new Interp3D(amap);
            m_Interp->debug(false);
            //m_Interp->debug();
        }

        /// setup working points
        if (m_phWPname != "Undefined") ATH_CHECK(addPhotonWP(m_phWPname));
        if (m_elWPname != "Undefined") ATH_CHECK(addElectronWP(m_elWPname));
        if (m_muWPname != "Undefined") ATH_CHECK(addMuonWP(m_muWPname));
	for (auto c: m_muWPvec) ATH_CHECK(addMuonWP(c));
	for (auto c: m_elWPvec) ATH_CHECK(addElectronWP(c));
	for (auto c: m_phWPvec) ATH_CHECK(addPhotonWP(c));

        /// Return gracefully:
        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::setIParticleCutsFrom(xAOD::Type::ObjectType ObjType) {
        if (ObjType == xAOD::Type::Electron) {
            m_iparAcceptInfo = &m_electronAccept;
            m_iparWPs = &m_elWPs;
        } else if (ObjType == xAOD::Type::Muon) {
            m_iparAcceptInfo = &m_muonAccept;
            m_iparWPs = &m_muWPs;
        } else if (ObjType == xAOD::Type::Photon) {
            m_iparAcceptInfo = &m_photonAccept;
            m_iparWPs = &m_phWPs;
        } else {
            return StatusCode::FAILURE;
        }

        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::finalize() {
        /// Return gracefully:
        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::addCutToWP(IsolationWP* wp, std::string key, const xAOD::Iso::IsolationType t, const std::string expression, const xAOD::Iso::IsolationType isoCutRemap) {
        if(m_calibFile==nullptr) {
	  ATH_MSG_ERROR("Calibration File (" << m_calibFileName << ") is missing.");
	  return StatusCode::FAILURE;
	}

        std::string varname(xAOD::Iso::toCString(isoCutRemap));
        key += varname;

        TDirectory* tempDir = getTemporaryDirectory();
        tempDir->cd();

        std::shared_ptr<TH3F> histogram((TH3F*) m_calibFile->Get(key.c_str())->Clone());
        IsolationConditionHist *ich = new IsolationConditionHist(varname, t, expression, histogram);
        if ((m_doInterpM && key.find("Muon") != std::string::npos) || (m_doInterpE && key.find("Electron") != std::string::npos)) ich->setInterp(m_Interp);
        wp->addCut(ich);

        return StatusCode::SUCCESS;
   }

    StatusCode IsolationSelectionTool::addCutToWP(IsolationWP* wp, std::string key, const xAOD::Iso::IsolationType t, const std::string expression) {
        return addCutToWP(wp,key,t,expression,t);
    }

    StatusCode IsolationSelectionTool::addMuonWP(std::string muWPname) {
        IsolationWP* wp = new IsolationWP(muWPname);
        if (muWPname == "HighPtTrackOnly") {
            wp->addCut(new IsolationConditionFormula("ptcone20_Tight_1p25", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "1.25E03"));  //units are MeV!
        } else if (muWPname == "TightTrackOnly_FixedRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTrackOnly_lowPt", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.06*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTrackOnly_highPt", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "0.06*(x>50e3?x:1e9)"));
        } else if (muWPname == "Tight_FixedRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTight_track_lowPt", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.04*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTight_track_highPt", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "0.04*(x>50e3?x:1e9)"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTight_calo", xAOD::Iso::topoetcone20, "0.15*x"));
        } else if (muWPname == "Loose_FixedRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuLoose_track_lowPt", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.15*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuLoose_track_highPt", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "0.15*(x>50e3?x:1e9)"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuLoose_calo", xAOD::Iso::topoetcone20, "0.30*x"));
        } else if (muWPname == "TightTrackOnly_VarRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTrackOnly", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.06*x"));
        } else if (muWPname == "Tight_VarRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTight_track", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.04*x"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuTight_calo", xAOD::Iso::topoetcone20, "0.15*x"));
        } else if (muWPname == "Loose_VarRad") {
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuLoose_track", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000, "0.15*x"));
            wp->addCut(new IsolationConditionFormula("MuonFixedCutHighMuLoose_calo", xAOD::Iso::topoetcone20, "0.30*x"));
        } else if (muWPname == "PflowTight_FixedRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypesHighPt;
            std::vector<xAOD::Iso::IsolationType> isoTypesLowPt;
            isoTypesHighPt.push_back(xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypesHighPt.push_back(xAOD::Iso::neflowisol20);
            isoTypesLowPt.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypesLowPt.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("MuonPFlowTightLowPt", isoTypesLowPt, TF2("pflowTFunctionLowPt","fabs(x)+0.4*(y>0?y:0)"), "0.045*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionCombined("MuonPFlowTightHighPt", isoTypesHighPt, TF2("pflowTFunctionHighPt","fabs(x)+0.4*(y>0?y:0)"), "0.045*(x>50e3?x:1e9)"));
        } else if (muWPname == "PflowTight_VarRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypes;
            isoTypes.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypes.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("MuonPFlowTight", isoTypes, TF2("pflowTFunction","fabs(x)+0.4*(y>0?y:0)"), "0.045*x"));
        } else if (muWPname == "PflowLoose_FixedRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypesHighPt;
            std::vector<xAOD::Iso::IsolationType> isoTypesLowPt;
            isoTypesHighPt.push_back(xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypesHighPt.push_back(xAOD::Iso::neflowisol20);
            isoTypesLowPt.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypesLowPt.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("MuonPFlowLooseLowPt", isoTypesLowPt, TF2("pflowLFunctionLowPt","fabs(x)+0.4*(y>0?y:0)"), "0.16*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionCombined("MuonPFlowLooseHighPt", isoTypesHighPt, TF2("pflowLFunctionHighPt","fabs(x)+0.4*(y>0?y:0)"), "0.16*(x>50e3?x:1e9)"));
        } else if (muWPname == "PflowLoose_VarRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypes;
            isoTypes.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt500);
            isoTypes.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("MuonPFlowLoose", isoTypes, TF2("pflowTFunction","fabs(x)+0.4*(y>0?y:0)"), "0.16*x"));
        } else {
            ATH_MSG_ERROR("Unknown muon isolation WP: " << muWPname);
            delete wp;
            return StatusCode::FAILURE;
        }
        m_muWPs.push_back(wp);
        m_muonAccept.addCut(wp->name(), wp->name());
        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::addPhotonWP(std::string phWPname) {
        auto wp = new IsolationWP(phWPname);
        if (phWPname == "TightCaloOnly") {
            wp->addCut(new IsolationConditionFormula("PhFixedCut_calo40", xAOD::Iso::topoetcone40, "0.022*x+2450"));
        } else if (phWPname == "FixedCutTight") {
            wp->addCut(new IsolationConditionFormula("PhFixedCut_calo40", xAOD::Iso::topoetcone40, "0.022*x+2450"));
            wp->addCut(new IsolationConditionFormula("PhFixedCut_track20", xAOD::Iso::ptcone20, "0.05*x"));
        } else if (phWPname == "FixedCutLoose") {
            wp->addCut(new IsolationConditionFormula("PhFixedCut_calo20", xAOD::Iso::topoetcone20, "0.065*x"));
            wp->addCut(new IsolationConditionFormula("PhFixedCut_track20", xAOD::Iso::ptcone20, "0.05*x"));
        } else if (phWPname == "Tight") {
            wp->addCut(new IsolationConditionFormula("PhFixedCut_calo40", xAOD::Iso::topoetcone40, "0.022*x+2450"));
            wp->addCut(new IsolationConditionFormula("PhFixedCut_Tighttrack20", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "0.05*x"));
        } else if (phWPname == "Loose") {
            wp->addCut(new IsolationConditionFormula("PhFixedCut_calo20", xAOD::Iso::topoetcone20, "0.065*x"));
            wp->addCut(new IsolationConditionFormula("PhFixedCut_Tighttrack20", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVA_pt1000, "0.05*x"));
        } else {
            ATH_MSG_ERROR("Unknown photon isolation WP: " << phWPname);
            delete wp;
            return StatusCode::FAILURE;
        }

        m_phWPs.push_back(wp);
        m_photonAccept.addCut(wp->name(), wp->name());

        // Return gracefully:
        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::addElectronWP(std::string elWPname) {
        auto wp = new IsolationWP(elWPname);

        if(elWPname == "HighPtCaloOnly"){
            wp->addCut(new IsolationConditionFormula("FCHighPtCaloOnly_calo",    xAOD::Iso::topoetcone20, "std::max(0.015*x,3.5E3)")); //units are MeV!
        } else if (elWPname == "Tight_VarRad") {
            wp->addCut(new IsolationConditionFormula("ElecTight_track", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000, "0.06*x"));
            wp->addCut(new IsolationConditionFormula("ElecTight_calo", xAOD::Iso::topoetcone20, "0.06*x"));
        } else if (elWPname == "Loose_VarRad") {
            wp->addCut(new IsolationConditionFormula("ElecLoose_track", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000, "0.15*x"));
            wp->addCut(new IsolationConditionFormula("ElecLoose_calo", xAOD::Iso::topoetcone20, "0.20*x"));
        } else if (elWPname == "TightTrackOnly_VarRad") {
            wp->addCut(new IsolationConditionFormula("ElecTightTrackOnly", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000, "0.06*x"));
        } else if (elWPname == "TightTrackOnly_FixedRad") {
            wp->addCut(new IsolationConditionFormula("ElecTightTrackOnly_lowPt", xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000, "0.06*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionFormula("ElecTightTrackOnly_highPt", xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVALooseCone_pt1000, "0.06*(x>50e3?x:1e9)"));
        } else if (elWPname == "PflowTight_FixedRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypesHighPt;
            std::vector<xAOD::Iso::IsolationType> isoTypesLowPt;
            isoTypesHighPt.push_back(xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypesHighPt.push_back(xAOD::Iso::neflowisol20);
            isoTypesLowPt.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypesLowPt.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("ElecPFlowTightLowPt", isoTypesLowPt, TF2("pflowTFunctionLowPt","fabs(x)+0.4*(y>0?y:0)"), "0.045*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionCombined("ElecPFlowTightHighPt", isoTypesHighPt, TF2("pflowTFunctionHighPt","fabs(x)+0.4*(y>0?y:0)"), "0.045*(x>50e3?x:1e9)"));
        } else if (elWPname == "PflowTight") {
            std::vector<xAOD::Iso::IsolationType> isoTypes;
            isoTypes.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypes.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("ElecPFlowTight", isoTypes, TF2("pflowLFunction","fabs(x)+0.4*(y>0?y:0)"), "0.045*x"));
        } else if (elWPname == "PflowLoose_FixedRad") {
            std::vector<xAOD::Iso::IsolationType> isoTypesHighPt;
            std::vector<xAOD::Iso::IsolationType> isoTypesLowPt;
            isoTypesHighPt.push_back(xAOD::Iso::ptcone20_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypesHighPt.push_back(xAOD::Iso::neflowisol20);
            isoTypesLowPt.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypesLowPt.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("ElecPFlowLooseLowPt", isoTypesLowPt, TF2("pflowLFunctionLowPt","fabs(x)+0.4*(y>0?y:0)"), "0.16*(x>50e3?1e9:x)"));
            wp->addCut(new IsolationConditionCombined("ElecPFlowLooseHighPt", isoTypesHighPt, TF2("pflowLFunctionHighPt","fabs(x)+0.4*(y>0?y:0)"), "0.16*(x>50e3?x:1e9)"));
        } else if (elWPname == "PflowLoose") {
            std::vector<xAOD::Iso::IsolationType> isoTypes;
            isoTypes.push_back(xAOD::Iso::ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt500);
            isoTypes.push_back(xAOD::Iso::neflowisol20);
            wp->addCut(new IsolationConditionCombined("ElecPFlowLoose", isoTypes, TF2("pflowLFunction","fabs(x)+0.4*(y>0?y:0)"), "0.16*x"));
	} else {
	  ATH_MSG_ERROR("Unknown electron isolation WP: " << elWPname);
            delete wp;
            return StatusCode::FAILURE;
        }

        m_elWPs.push_back(wp);
        m_electronAccept.addCut(wp->name(), wp->name());

        // Return gracefully:
        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::addUserDefinedWP(std::string WPname, xAOD::Type::ObjectType ObjType, std::vector<std::pair<xAOD::Iso::IsolationType, std::string> >& cuts, std::string key, IsoWPType type) {
        std::vector<IsolationWP*>* wps(nullptr);
        asg::AcceptInfo *ac = nullptr;
        if (ObjType == xAOD::Type::Electron) {
            if (key == "") key = m_elWPKey;
            wps = &m_elWPs;
            ac = &m_electronAccept;
        } else if (ObjType == xAOD::Type::Muon) {
            if (key == "") key = m_muWPKey;
            wps = &m_muWPs;
            ac = &m_muonAccept;
        } else if (ObjType == xAOD::Type::Photon) {
            if (key == "") key = m_phWPKey;
            wps = &m_phWPs;
            ac = &m_photonAccept;
        } else if (ObjType == xAOD::Type::Other) {
            if (key == "") return StatusCode::FAILURE;
            wps = &m_objWPs;
            ac = &m_objAccept;
        } else {
            return StatusCode::FAILURE;
        }

        auto wp = new IsolationWP(WPname);
        if (type == Efficiency) {
            for (auto& c : cuts)
	      ATH_CHECK(addCutToWP(wp, key, c.first, c.second));
        } else if (type == Cut) {
            for (auto& c : cuts)
                wp->addCut(new IsolationConditionFormula(xAOD::Iso::toCString(c.first), c.first, c.second));
        } else {
            ATH_MSG_ERROR("Unknown isolation WP type -- should not happen.");
            delete wp;
            return StatusCode::FAILURE;
        }

        wps->push_back(wp);
        ac->addCut(wp->name(), wp->name());

        return StatusCode::SUCCESS;
    }

    StatusCode IsolationSelectionTool::addWP(std::string WP, xAOD::Type::ObjectType ObjType) {
        if (ObjType == xAOD::Type::Electron) {
            return addElectronWP(WP);
        } else if (ObjType == xAOD::Type::Muon) {
            return addMuonWP(WP);
        } else if (ObjType == xAOD::Type::Photon) {
            return addPhotonWP(WP);
        }

        return StatusCode::FAILURE;
    }
    StatusCode IsolationSelectionTool::addWP(IsolationWP* wp, xAOD::Type::ObjectType ObjType) {
        if (ObjType == xAOD::Type::Electron) {
            m_elWPs.push_back(wp);
            m_electronAccept.addCut(wp->name(), wp->name());
        } else if (ObjType == xAOD::Type::Muon) {
            m_muWPs.push_back(wp);
            m_muonAccept.addCut(wp->name(), wp->name());
        } else if (ObjType == xAOD::Type::Photon) {
            m_phWPs.push_back(wp);
            m_photonAccept.addCut(wp->name(), wp->name());
        } else if (ObjType == xAOD::Type::Other) {
            m_objWPs.push_back(wp);
            m_objAccept.addCut(wp->name(), wp->name());
        } else {
            return StatusCode::FAILURE;
        }

        return StatusCode::SUCCESS;
    }
    template<typename T> void IsolationSelectionTool::evaluateWP(const T& x, const std::vector<IsolationWP*>& WP, asg::AcceptData& accept) const {
        accept.clear();
        for (auto& i : WP) {
            if (i->accept(x)) accept.setCutResult(i->name(), true);
        }
    }
    asg::AcceptData IsolationSelectionTool::accept(const xAOD::Photon& x) const {
        asg::AcceptData accept (&m_photonAccept);
        evaluateWP(x, m_phWPs, accept);
        return accept;
    }

    asg::AcceptData IsolationSelectionTool::accept(const xAOD::Electron& x) const {
        asg::AcceptData accept (&m_electronAccept);
        evaluateWP(x, m_elWPs, accept);
        return accept;
    }

    asg::AcceptData IsolationSelectionTool::accept(const xAOD::Muon& x) const {
        asg::AcceptData accept (&m_muonAccept);
        evaluateWP(x, m_muWPs, accept);
        return accept;
    }

    asg::AcceptData IsolationSelectionTool::accept(const xAOD::IParticle& x) const {

        if (x.type() == xAOD::Type::Electron) {
            asg::AcceptData accept (&m_electronAccept);
            evaluateWP(x, m_elWPs, accept);
            return accept;
        } else if (x.type() == xAOD::Type::Muon) {
            asg::AcceptData accept (&m_muonAccept);
            evaluateWP(x, m_muWPs, accept);
            return accept;
        } else if (x.type() == xAOD::Type::Photon) {
            asg::AcceptData accept (&m_photonAccept);
            evaluateWP(x, m_phWPs, accept);
            return accept;
        }

        else if (m_iparAcceptInfo && m_iparWPs) {
          asg::AcceptData accept (m_iparAcceptInfo);
            evaluateWP(x, *m_iparWPs, accept);
            return accept;
        }
        ATH_MSG_ERROR("Someting here makes really no  sense");
        return asg::AcceptData (&m_objAccept);
    }

    asg::AcceptData IsolationSelectionTool::accept(const strObj& x) const {
        if (x.type == xAOD::Type::Electron) {
            asg::AcceptData accept (&m_electronAccept);
            evaluateWP(x, m_elWPs, accept);
            return accept;
        } else if (x.type == xAOD::Type::Muon) {
            asg::AcceptData accept (&m_muonAccept);
            evaluateWP(x, m_muWPs, accept);
            return accept;
        } else if (x.type == xAOD::Type::Photon) {
            asg::AcceptData accept (&m_photonAccept);
            evaluateWP(x, m_phWPs, accept);
            return accept;
        } else {
            asg::AcceptData accept (&m_objAccept);
            evaluateWP(x, m_objWPs, accept);
            return accept;
        }
        return asg::AcceptData (&m_objAccept);
    }

    const asg::AcceptInfo& IsolationSelectionTool::getPhotonAcceptInfo() const {
        return m_photonAccept;
    }

    const asg::AcceptInfo& IsolationSelectionTool::getElectronAcceptInfo() const {
        return m_electronAccept;
    }

    const asg::AcceptInfo& IsolationSelectionTool::getMuonAcceptInfo() const {
        return m_muonAccept;
    }
    const asg::AcceptInfo& IsolationSelectionTool::getObjAcceptInfo() const {
        return m_objAccept;
    }

    TDirectory* IsolationSelectionTool::getTemporaryDirectory(void) const {
        //
        // Create a unique directory in memory to hold the histograms:
        //
        gROOT->cd();
        TDirectory* tempDir = 0;
        int counter = 0;
        while (not tempDir) {
            // First, let's find a directory name that doesn't exist yet:
            std::stringstream dirname;
            dirname << "IsolationSelectionTempDir_%i" << counter;
            if (gROOT->GetDirectory((dirname.str()).c_str())) {
                ++counter;
                continue;
            }
            // Let's try to make this directory:
            tempDir = gROOT->mkdir((dirname.str()).c_str());
            if (not tempDir) {
                //      ::Fatal("IsolationSelectionTool::getTemporaryDirectory",
                //      ERROR_MESSAGE("Temporary directory could not be created"));
                throw std::runtime_error("Temporary directory could not be created");
            }
        }

        return tempDir;
    }
}
