/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


#include "HIEventUtils/ZdcRecTool.h"
#include "TGraph.h"
#include "TEnv.h"
#include "TSystem.h"

namespace ZDC
{

  ZdcRecTool::ZdcRecTool(const std::string& name) : asg::AsgTool(name),
              m_name(name),
              m_init(false)
  {
    declareProperty("ZdcModuleContainerName",m_zdcModuleContainerName="ZdcModules","Location of ZDC processed data");
    declareProperty("ZdcRecConfigPath",m_zdcRecConfigPath="$ROOTCOREDIR/data/HIEventUtils/","ZDC Rec config file path");

    ATH_MSG_DEBUG("Creating ZdcRecoTool named " << m_name);
    ATH_MSG_INFO("ZDC config file path " << m_zdcRecConfigPath);

  }

  ZdcRecTool::~ZdcRecTool()
  {
    ATH_MSG_DEBUG("Deleting ZdcRecoTool named " << m_name);
  }

  StatusCode ZdcRecTool::initializeTool()
  {
    m_init = true;
    m_tf1SincInterp = new TF1("SincInterp",SincInterp,-5.,160.,8);
    m_tf1SincInterp->SetNpx(300);
    m_tf1FermiExpFit = new TF1("FermiExpFit",FermiExpFit,-5.,160.,7);
    m_tf1FermiExpFit->SetNpx(300);

    char* path = gSystem->ExpandPathName(m_zdcRecConfigPath.c_str());
    ATH_MSG_INFO("Resolved file path " << path);

    TString zdcConf(path);
    zdcConf += "/ZdcRecConfig.conf";

    TEnv env(zdcConf);

    ATH_MSG_INFO("ZdcC calibration LG " << env.GetValue("ZdcCCalibrationLG","1.;1.;1.;1."));

    return StatusCode::SUCCESS;

  }

  StatusCode ZdcRecTool::recoZdcModule(const xAOD::ZdcModule& module)
  {

    ATH_MSG_DEBUG("Not processing ZDC module S/T/M/C = "
    << module.zdcSide() << " "
    << module.zdcType() << " "
    << module.zdcModule() << " "
    << module.zdcChannel()
  );


  return StatusCode::SUCCESS;
}

StatusCode ZdcRecTool::recoZdcModules(const xAOD::ZdcModuleContainer& moduleContainer)
{
  if (!m_eventReady)
  {
    ATH_MSG_INFO("Event not ready for ZDC reco!");
    return StatusCode::FAILURE;
  }

  for (const auto zdcModule : moduleContainer)
  {
    ATH_CHECK(recoZdcModule(*zdcModule));
  }
  return StatusCode::SUCCESS;
}

StatusCode ZdcRecTool::reprocessZdc()
{
  if (!m_init)
  {
    ATH_MSG_INFO("Tool not initialized!");
    return StatusCode::FAILURE;
  }
  m_eventReady = false;
  ATH_CHECK(evtStore()->retrieve(m_zdcModules,m_zdcModuleContainerName));
  m_eventReady = true;

  ATH_CHECK(recoZdcModules(*m_zdcModules));

  return StatusCode::SUCCESS;
}

bool ZdcRecTool::sigprocMaxFinder(const std::vector<unsigned short>& adc, float deltaT, float& amp, float& time, float& qual)
{
  size_t nsamp = adc.size();
  float presamp = adc.at(0);
  unsigned short max_adc = 0;
  int max_index = -1;
  for (size_t i = 0;i<nsamp;i++)
  {
    if (adc[i]>max_adc)
    {
      max_adc = adc[i];
      max_index = i;
    }
  }
  amp = max_adc - presamp;
  time = max_index*deltaT;
  qual = 1.;

  if(max_index==-1)
  {
    qual=0.;
    return false;
  }

  return true;
}

bool ZdcRecTool::sigprocPeakFitter(const std::vector<unsigned short>& adc, float deltaT, float& amp, float& time, float& qual)
{
  TGraph g;
  size_t nsamp = adc.size();
  float presamp = adc.at(0);
  unsigned short max_adc = 0;
  for (size_t i = 0;i<nsamp;i++)
  {
    if (adc.at(i)>max_adc)
    {
      max_adc = adc.at(i);
    }
    g.SetPoint(i,i*deltaT,adc.at(i)-presamp);
  }
  double timeScale = deltaT/12.5;
  m_tf1FermiExpFit->SetParameters(max_adc-presamp,30*timeScale, 2.5*timeScale, 25*timeScale);
  m_tf1FermiExpFit->SetParLimits(0,0,1024);
  m_tf1FermiExpFit->SetParLimits(1,10*timeScale,50*timeScale);
  m_tf1FermiExpFit->FixParameter(2,2.5*timeScale);
  m_tf1FermiExpFit->SetParLimits(3,10*timeScale,40*timeScale);
  g.Fit("FermiExpFit","QWR","");
  amp = m_tf1FermiExpFit->GetMaximum();
  time = m_tf1FermiExpFit->GetMaximumX();
  qual = 1.;
  return true;
}

bool ZdcRecTool::sigprocSincInterp(const std::vector<unsigned short>& adc, float deltaT, float& amp, float& time, float& qual)
{
  size_t nsamp = adc.size();
  float presamp = adc.at(0);
  m_tf1SincInterp->SetParameter(0,deltaT);
  for (size_t i = 0;i<nsamp;i++)
  {
    m_tf1SincInterp->SetParameter(i+1,adc.at(i)-presamp);
  }
  amp = m_tf1SincInterp->GetMaximum();
  time = m_tf1SincInterp->GetMaximumX();
  qual = 1.;
  return true;
}

double SincInterp(double* xvec, double* pvec)
{
  // pvec are the sample values
  double ret = 0;
  double T = pvec[0]; // deltaT
  double t = xvec[0];
  for (int isamp = 0;isamp<7;isamp++)
  {
    double arg = (t - isamp*T)/T;
    if (arg!=0.0)
    {
      ret += pvec[isamp+1] * std::sin(TMath::Pi()*arg)/(TMath::Pi()*arg);
    }
  }
  return ret;
}

double FermiExpFit(double* xvec, double* pvec)
{
  const float offsetScale = 3;

  double t = xvec[0];

  double amp = pvec[0];
  double t0 = pvec[1];
  double tau1 = pvec[2];
  double tau2 = pvec[3];

  double tauRatio = tau2/tau1;
  double tauRatioMinunsOne = tauRatio - 1;

  double norm = ( std::exp(-offsetScale/tauRatio)*pow(1./tauRatioMinunsOne, 1./(1 + tauRatio))/
  ( 1 + pow(1./tauRatioMinunsOne, 1./(1 + 1/tauRatio)))   );

  double deltaT = t - (t0 - offsetScale*tau1);
  if (deltaT < 0) deltaT = 0;

  double expTerm = std::exp(-deltaT/tau2);
  double fermiTerm = 1./(1. + std::exp(-(t - t0)/tau1));

  return amp*expTerm*fermiTerm/norm;
}

} // namespace ZDC

