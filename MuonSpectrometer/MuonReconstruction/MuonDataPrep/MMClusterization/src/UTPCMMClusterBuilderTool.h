/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/
#ifndef UTPCMMClusterBuilderTool_h
#define UTPCMMClusterBuilderTool_h

#include "GaudiKernel/ToolHandle.h"
#include "MMClusterization/IMMClusterBuilderTool.h"
#include "MuonPrepRawData/MMPrepData.h"
#include "AthenaBaseComps/AthAlgTool.h"

#include <numeric>

#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLinearFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"

class MmIdHelper;
namespace MuonGM
{
  class MuonDetectorManager;
}

// 
// Simple clusterization tool for MicroMegas
//
namespace Muon
{
  
  class UTPCMMClusterBuilderTool : virtual public IMMClusterBuilderTool, public AthAlgTool {

  public:
    /** Default constructor */
    UTPCMMClusterBuilderTool(const std::string&, const std::string&, const IInterface*);
     
    /** Default destructor */
    virtual ~UTPCMMClusterBuilderTool() = default;

    /** standard initialize method */
    virtual StatusCode initialize() override;

    virtual
    StatusCode getClusters(std::vector<Muon::MMPrepData>& MMprds, 
			   std::vector<Muon::MMPrepData*>& clustersVec) override;

  private: 

    /// Muon Detector Descriptor
    const MuonGM::MuonDetectorManager* m_muonMgr;
    const MmIdHelper* m_mmIdHelper;


    // params for the hough trafo
    double m_alphaMin,m_alphaMax,m_alphaResolution,m_selectionCut;
    double m_dMin,m_dMax,m_dResolution;
    int m_houghMinCounts;

    double m_timeOffset,m_dHalf,m_vDrift;

    double m_toRad=TMath::Pi()/180.;


    StatusCode runHoughTrafo(std::vector<int>& flag,std::vector<double>& xpos, std::vector<double>& time,std::vector<int>& idx_selected);
    StatusCode fillHoughTrafo(std::unique_ptr<TH2D>& cummulator,std::vector<int>& flag, std::vector<double>& xpos, std::vector<double>& time);
    StatusCode houghInitCummulator(std::unique_ptr<TH2D>& cummulator,double xmax,double xmin);

    StatusCode findAlphaMax(std::unique_ptr<TH2D>& h_hough, std::vector<std::tuple<double,double>> &maxPos);
    StatusCode selectTrack(std::vector<std::tuple<double,double>> &tracks,std::vector<double>& xpos, std::vector<double>& time,std::vector<int>& flag,std::vector<int>& idx_selected);

    StatusCode transformParameters(double alpha, double d, double dRMS, double& slope,double& intercept, double& interceptRMS);
    StatusCode finalFit(std::vector<double>& xpos, std::vector<double>& time, std::vector<int>& idxSelected,double& x0, double &sigmaX0, double &fitAngle, double &chiSqProb);
};


}
#endif
