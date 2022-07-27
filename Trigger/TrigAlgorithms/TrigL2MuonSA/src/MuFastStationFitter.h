/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef  TRIGL2MUONSA_MUFASTSTATIONFITTER_H
#define  TRIGL2MUONSA_MUFASTSTATIONFITTER_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "MdtData.h"
#include "TgcFitResult.h"
#include "RpcFitResult.h"
#include "SuperPointData.h"
#include "TrackData.h"
#include "MuonRoad.h"
#include "PtEndcapLUT.h"
#include "TrigMuonToolInterfaces/ITrigMuonBackExtrapolator.h"
#include "AlphaBetaEstimate.h"
#include "PtFromAlphaBeta.h"
#include "NswStationFitter.h"
#include "StgcData.h"
#include "MmData.h"

namespace TrigL2MuonSA {

class MuFastStationFitter: public AthAlgTool
{
   public:

      MuFastStationFitter(const std::string& type,
                          const std::string& name,
                          const IInterface*  parent);
    
      virtual StatusCode initialize() override;
    
   public:
      StatusCode findSuperPoints(const TrigRoiDescriptor* p_roids,
				 const TrigL2MuonSA::MuonRoad& muonRoad,
				 TrigL2MuonSA::RpcFitResult& rpcFitResult,
				 std::vector<TrigL2MuonSA::TrackPattern>& v_trackPatterns) const;
      StatusCode findSuperPointsSimple(const TrigRoiDescriptor* p_roids,
				       const TrigL2MuonSA::MuonRoad& muonRoad,
				       TrigL2MuonSA::TgcFitResult& tgcFitResult,
				       std::vector<TrigL2MuonSA::TrackPattern>& v_trackPatterns,
				       TrigL2MuonSA::StgcHits& stgcHits,
				       TrigL2MuonSA::MmHits& mmHits) const;

      StatusCode findSuperPoints(const TrigRoiDescriptor* p_roids,
                                 const TrigL2MuonSA::MuonRoad& muonRoad,
                                 TrigL2MuonSA::TgcFitResult& tgcFitResult,
                                 std::vector<TrigL2MuonSA::TrackPattern>& v_trackPatterns,
                                 TrigL2MuonSA::StgcHits& stgcHits,
                                 TrigL2MuonSA::MmHits& mmHits) const;

      StatusCode superPointFitter(TrigL2MuonSA::TrackPattern& trackPattern) const;

      StatusCode superPointFitter(TrigL2MuonSA::TrackPattern& trackPattern,
                                  const TrigL2MuonSA::MuonRoad&  muonRoad) const;

      StatusCode setMCFlag(const BooleanProperty& use_mcLUT);

   private:

      BooleanProperty m_use_mcLUT {true};
		
      Gaudi::Property< double > m_endcapinn_mdt_chi2_limit {
	this, "ENDCAPINN_MDT_CHI2_LIMIT", 20., ""};
      Gaudi::Property< double > m_endcapmid_mdt_chi2_limit {
	this, "ENDCAPMID_MDT_CHI2_LIMIT", 20., ""};
      Gaudi::Property< double > m_endcapout_mdt_chi2_limit {
	this, "ENDCAPOUT_MDT_CHI2_LIMIT", 20., ""};
      Gaudi::Property< double > m_endcapee_mdt_chi2_limit {
	this, "ENDCAPEE_MDT_CHI2_LIMIT",  20., ""};

      Gaudi::Property< double > m_rwidth_Endcapinn_first {
	this, "RWIDTH_EndcapINN_FIRST",  150., ""};
      Gaudi::Property< double > m_rwidth_Endcapinn_second {
	this, "RWIDTH_EndcapINN_SECOND", 80., ""};
      Gaudi::Property< double > m_rwidth_Endcapmid_first {
	this, "RWIDTH_EndcapMID_FIRST", 150., ""};
      Gaudi::Property< double > m_rwidth_Endcapmid_second {
	this, "RWIDTH_EndcapMID_SECOND", 100., ""};
      Gaudi::Property< double > m_rwidth_Endcapout_first {
	this, "RWIDTH_EndcapOUT_FIRST", 120., ""};
      Gaudi::Property< double > m_rwidth_Endcapout_second {
	this, "RWIDTH_EndcapOUT_SECOND", 60., ""};
      Gaudi::Property< double > m_rwidth_Endcapee_first {
	this, "RWIDTH_EndcapEE_FIRST", 150., ""};
      Gaudi::Property< double > m_rwidth_Endcapee_second {
	this, "RWIDTH_EndcapEE_SECOND", 100., ""};

      Gaudi::Property< double > m_mdt_driftspace_uplimit {
	this, "MDT_DRFITSPACE_UPLIMIT", 14.8, ""};
      Gaudi::Property< double > m_mdt_driftspace_downlimit {
	this, "MDT_DRFITSPACE_DOWNLIMIT", 0.1, ""};
      Gaudi::Property< double > m_mdt_drifttime_limit {
	this, "MDT_DRFITTIME_LIMIT", 1700., ""};

      ToolHandle<ITrigMuonBackExtrapolator> m_backExtrapolator {
	this, "BackExtrapolator", "TrigMuonBackExtrapolator", "public tool for back extrapolating the muon tracks to the IV"};

   private:
      float SetDriftSpace(float tdr, float rad, float zeta, float phim, float phiDir) const;
      void  Xline(float *, float *, float *, int *, int ,
                  float *, float *, float *, float *, float *, float *) const;
      void  Circfit (int, float *, float *, float *, float *, int *,
                     float *, float *, float DAB[2][2], float *)  const;
      void  Circles (int, float *, float *, float *, float *, int *,
                     float *, float *, float DAB[2][2], float *, float *) const;
      int   Evlfit (int, TrigL2MuonSA::PBFitResult& fitres) const;

      ToolHandle<AlphaBetaEstimate>          m_alphaBetaEstimate {
	this, "AlphaBetaEstimate", "TrigL2MuonSA::AlphaBetaEstimate"};
      ToolHandle<PtFromAlphaBeta>            m_ptFromAlphaBeta {
	this, "PtFromAlphaBeta", "TrigL2MuonSA::PtFromAlphaBeta", ""};
      ToolHandle<NswStationFitter> m_nswStationFitter {this, "NswStationFitter", "TrigL2MuonSA::NswStationFitter"};

      void findLayerCombination(std::vector<unsigned int> &a, int n, int r,std::vector<std::vector<unsigned int> > &c, int &nr) const;
      void findSubLayerCombination(std::vector<unsigned int> &a, int n,int r, std::vector<unsigned int> &b, int index ,int num,
                                   std::vector<std::vector<unsigned int> > &c, int &nr) const;
      void makeReferenceLine(TrigL2MuonSA::TrackPattern& trackPattern,const TrigL2MuonSA::MuonRoad&    muonRoad) const;
      void Circles (int, float *, float *, float *, float *, int *,
                     float *, float *, float DAB[2][2], float *, float *, float *, float *, float *) const;

      double fromAlphaPtToInn(TrigL2MuonSA::TgcFitResult& tgcFitResult,TrigL2MuonSA::TrackPattern& trackPattern) const;
      void updateInnSP(TrigL2MuonSA::TrackPattern& trackPattern, double &aw,double &tgc_aw, double &bw) const;
      void stationSPFit(TrigL2MuonSA::MdtHits*    mdtSegment, TrigL2MuonSA::SuperPoint* superPoint,
                        TrigL2MuonSA::PBFitResult& pbFitResult, int s_address, int i_station,double aw, float phiDir) const;


};

} // namespace TrigL2MuonSA

#endif  // MUFASTSTATIONFITTER_H
