/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

//*************************************************************************
//                          eFEXtauAlgo  -  description
//                             --------------------
//    begin                 : 06 05 2020
//    email                 : nicholas.andrew.luongo@cern.ch
//*************************************************************************

#include <iostream>

#include "L1CaloFEXSim/eFEXtauAlgo.h"
#include "L1CaloFEXSim/eFEXtauTOB.h"
#include "L1CaloFEXSim/eTower.h"

  // default constructor for persistency
LVL1::eFEXtauAlgo::eFEXtauAlgo(const std::string& type, const std::string& name, const IInterface* parent):
    AthAlgTool(type, name, parent)
  {
    declareInterface<IeFEXtauAlgo>(this);
  }

  /** Destructor */
LVL1::eFEXtauAlgo::~eFEXtauAlgo()
{
}

StatusCode LVL1::eFEXtauAlgo::initialize()
{
  ATH_CHECK(m_eTowerContainerKey.initialize());
  return StatusCode::SUCCESS;
}

StatusCode LVL1::eFEXtauAlgo::safetyTest(){

  SG::ReadHandle<eTowerContainer> eTowerContainer(m_eTowerContainerKey/*,ctx*/);
  if(!eTowerContainer.isValid()){
    ATH_MSG_FATAL("Could not retrieve eTowerContainer " << m_eTowerContainerKey.key() );
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;

}

void LVL1::eFEXtauAlgo::setup(int inputTable[3][3], int efex_id, int fpga_id, int central_eta){

  std::copy(&inputTable[0][0], &inputTable[0][0] + 9, &m_eFexalgoTowerID[0][0]);

  buildLayers(efex_id, fpga_id, central_eta);
  setSupercellSeed();
  setUnDAndOffPhi();

}

std::unique_ptr<LVL1::eFEXtauTOB> LVL1::eFEXtauAlgo::getTauTOB()
{
  std::unique_ptr<eFEXtauTOB> tob = std::make_unique<eFEXtauTOB>();
  unsigned int et = getEt();
  tob->setEt(et);
  tob->setRcoreCore(rCoreCore());
  tob->setRcoreEnv(rCoreEnv());
  tob->setRhadCore(rHadCore());
  tob->setRhadEnv(rHadEnv());
  tob->setBitwiseEt(getBitwiseEt());
  tob->setIso(getRealRCore());
  tob->setSeedUnD(getUnD());
  return tob;
}

// Build arrays holding cell ETs for each layer plus entire tower
void LVL1::eFEXtauAlgo::buildLayers(int efex_id, int fpga_id, int central_eta)
{

  SG::ReadHandle<eTowerContainer> eTowerContainer(m_eTowerContainerKey/*,ctx*/);

  for(unsigned int ieta = 0; ieta < 3; ieta++)
  {
    for(unsigned int iphi = 0; iphi < 3; iphi++)
    {
      if (((efex_id%3 == 0) && (fpga_id == 0) && (central_eta == 0) && (ieta == 0)) || ((efex_id%3 == 2) && (fpga_id == 3) && (central_eta == 5) && (ieta == 2))) 
      {
	      m_twrcells[ieta][iphi] = 0;
	      m_em0cells[ieta][iphi] = 0;
	      m_em3cells[ieta][iphi] = 0;
	      m_hadcells[ieta][iphi] = 0;
	      for(unsigned int i = 0; i < 4; i++)
        {  
	        m_em1cells[4 * ieta + i][iphi] = 0;
	        m_em2cells[4 * ieta + i][iphi] = 0;
        }
      } 
      else 
      {  
	      const LVL1::eTower * tmpTower = eTowerContainer->findTower(m_eFexalgoTowerID[iphi][ieta]);
	      m_twrcells[ieta][iphi] = tmpTower->getTotalET();
	      m_em0cells[ieta][iphi] = tmpTower->getLayerTotalET(0);
	      m_em3cells[ieta][iphi] = tmpTower->getLayerTotalET(3);
	      m_hadcells[ieta][iphi] = tmpTower->getLayerTotalET(4);
	      for(unsigned int i = 0; i < 4; i++)
        {  
	        m_em1cells[4 * ieta + i][iphi] = tmpTower->getET(1, i);
	        m_em2cells[4 * ieta + i][iphi] = tmpTower->getET(2, i);
        }
      }
    }
  }
  m_cellsSet = true;
}

// Check if central tower qualifies as a seed tower for the tau algorithm
bool LVL1::eFEXtauAlgo::isCentralTowerSeed()
{
  // Need layer cell ET arrays to be built
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot accurately determine if a seed tower.");
  }
  
  bool out = true;
  
  // Get central tower ET
  unsigned int centralET = m_twrcells[1][1];
  
  // Loop over all cells and check that the central tower is a local maximum
  for (unsigned int beta = 0; beta < 3; beta++)
  {
    for (unsigned int bphi = 0; bphi < 3; bphi++)
    {
      // Don't need to compare central cell with itself
      if ((beta == 1) && (bphi == 1)){
        continue;
      }
      
      // Cells to the up and right must have strictly lesser ET
      if (((beta == 2) && (bphi == 0)) || ((beta == 2) && (bphi == 1)) || ((beta == 2) && (bphi == 2)) || ((beta == 1) && (bphi == 2)))
      {
        if (centralET <= m_twrcells[beta][bphi])
        {
          out = false;
        }
      }
      // Cells down and to the left must have lesser or equal ET. If strictly lesser would create zero TOB if two adjacent cells had equal energy
      else if (((beta == 0) && (bphi == 0)) || ((beta == 0) && (bphi == 1)) || ((beta == 0) && (bphi == 2)) || ((beta == 1) && (bphi == 0)))
      { 
        if (centralET < m_twrcells[beta][bphi])
        {
          out = false;
        }
      }
    }
  }
  
  return out;
}

// Calculate reconstructed ET value
unsigned int LVL1::eFEXtauAlgo::getEt()
{
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot accurately calculate Et.");
  }

  unsigned int out = 0;

  out += m_em0cells[0][1];
  out += m_em0cells[1][1];
  out += m_em0cells[2][1];
  out += m_em0cells[0][m_offPhi];
  out += m_em0cells[1][m_offPhi];
  out += m_em0cells[2][m_offPhi];

  out += m_em1cells[m_seed][1];
  out += m_em1cells[m_seed + 1][1];
  out += m_em1cells[m_seed + 2][1];
  out += m_em1cells[m_seed - 1][1];
  out += m_em1cells[m_seed - 2][1];
  out += m_em1cells[m_seed][m_offPhi];
  out += m_em1cells[m_seed + 1][m_offPhi];
  out += m_em1cells[m_seed + 2][m_offPhi];
  out += m_em1cells[m_seed - 1][m_offPhi];
  out += m_em1cells[m_seed - 2][m_offPhi];

  out += m_em2cells[m_seed][1];
  out += m_em2cells[m_seed + 1][1];
  out += m_em2cells[m_seed + 2][1];
  out += m_em2cells[m_seed - 1][1];
  out += m_em2cells[m_seed - 2][1];
  out += m_em2cells[m_seed][m_offPhi];
  out += m_em2cells[m_seed + 1][m_offPhi];
  out += m_em2cells[m_seed + 2][m_offPhi];
  out += m_em2cells[m_seed - 1][m_offPhi];
  out += m_em2cells[m_seed - 2][m_offPhi];

  out += m_em3cells[0][1];
  out += m_em3cells[1][1];
  out += m_em3cells[2][1];
  out += m_em3cells[0][m_offPhi];
  out += m_em3cells[1][m_offPhi];
  out += m_em3cells[2][m_offPhi];

  out += m_hadcells[0][1];
  out += m_hadcells[1][1];
  out += m_hadcells[2][1];
  out += m_hadcells[0][m_offPhi];
  out += m_hadcells[1][m_offPhi];
  out += m_hadcells[2][m_offPhi];

  return out;
}

unsigned int LVL1::eFEXtauAlgo::rCoreCore()
{
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot calculate rCore core value");
  }
     
  unsigned int out = 0;

  out += m_em2cells[m_seed][1];
  out += m_em2cells[m_seed + 1][1];
  out += m_em2cells[m_seed - 1][1];
  out += m_em2cells[m_seed][m_offPhi];
  out += m_em2cells[m_seed + 1][m_offPhi];
  out += m_em2cells[m_seed - 1][m_offPhi];

  return out;

}

unsigned int LVL1::eFEXtauAlgo::rCoreEnv()
{
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot calculate rCore environment value");
  }
     
  unsigned int out = 0;

  out += m_em2cells[m_seed + 2][1];
  out += m_em2cells[m_seed - 2][1];
  out += m_em2cells[m_seed + 3][1];
  out += m_em2cells[m_seed - 3][1];
  out += m_em2cells[m_seed + 4][1];
  out += m_em2cells[m_seed - 4][1];
  out += m_em2cells[m_seed + 2][m_offPhi];
  out += m_em2cells[m_seed - 2][m_offPhi];
  out += m_em2cells[m_seed + 3][m_offPhi];
  out += m_em2cells[m_seed - 3][m_offPhi];
  out += m_em2cells[m_seed + 4][m_offPhi];
  out += m_em2cells[m_seed - 4][m_offPhi];

  return out;

}

void LVL1::eFEXtauAlgo::getRCore(std::vector<unsigned int> & rCoreVec)
{
  unsigned int core = rCoreCore();
  unsigned int env = rCoreEnv();

  rCoreVec.push_back(core);
  rCoreVec.push_back(env);

}

// Calculate float isolation variable
float LVL1::eFEXtauAlgo::getRealRCore()
{
  unsigned int core = rCoreCore();
  unsigned int env = rCoreEnv();

  unsigned int num = core;
  unsigned int denom = core + env;

  float out = denom ? (float)num / (float)denom : 0;

  return out;
}

unsigned int LVL1::eFEXtauAlgo::rHadCore()
{
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot calculate rHad core value");
  }
     
  unsigned int out = 0;

  out += m_hadcells[0][1];
  out += m_hadcells[1][1];
  out += m_hadcells[2][1];
  out += m_hadcells[0][m_offPhi];
  out += m_hadcells[1][m_offPhi];
  out += m_hadcells[2][m_offPhi];

  return out;

}

unsigned int LVL1::eFEXtauAlgo::rHadEnv()
{
  if (m_cellsSet == false){
    ATH_MSG_DEBUG("Layers not built, cannot calculate rHad environment value");
  }
     
  unsigned int out = 0;

  out += m_em2cells[m_seed][1];
  out += m_em2cells[m_seed - 1][1];
  out += m_em2cells[m_seed + 1][1];
  out += m_em2cells[m_seed - 2][1];
  out += m_em2cells[m_seed + 2][1];
  out += m_em2cells[m_seed][m_offPhi];
  out += m_em2cells[m_seed - 1][m_offPhi];
  out += m_em2cells[m_seed + 1][m_offPhi];
  out += m_em2cells[m_seed - 2][m_offPhi];
  out += m_em2cells[m_seed + 2][m_offPhi];
  out += m_em1cells[m_seed][1];
  out += m_em1cells[m_seed - 1][1];
  out += m_em1cells[m_seed + 1][1];
  out += m_em1cells[m_seed][m_offPhi];
  out += m_em1cells[m_seed - 1][m_offPhi];
  out += m_em1cells[m_seed + 1][m_offPhi];

  return out;

}

// Calculate the hadronic fraction isolation variable
void LVL1::eFEXtauAlgo::getRHad(std::vector<unsigned int> & rHadVec)
{
  unsigned int core = rHadCore();
  unsigned int env = rHadEnv();

  rHadVec.push_back(core);
  rHadVec.push_back(env);

}

float LVL1::eFEXtauAlgo::getRealRHad()
{
  unsigned int core = rHadCore();
  unsigned int env = rHadEnv();

  unsigned int num = core;
  unsigned int denom = core + env;

  float out = denom ? (float)num / (float)denom : 0;

  return out;

}
// Set the off phi value used to calculate ET and isolation
void LVL1::eFEXtauAlgo::setUnDAndOffPhi()
{
  if (m_cellsSet == false){ 
    ATH_MSG_DEBUG("Layers not built, cannot accurately set phi direction.");
  }  

  unsigned int upwardEt = 0;
  upwardEt += m_em2cells[m_seed][2];
  upwardEt += m_em1cells[m_seed][2];
  
  unsigned int downwardEt = 0;
  downwardEt += m_em2cells[m_seed][0];
  downwardEt += m_em1cells[m_seed][0];

  if (downwardEt > upwardEt) {
    m_offPhi = 0;
    m_und = false;
  }
  else {
    m_offPhi = 2;
    m_und = true;
  }
}

// Find the supercell seed eta value, must be in central cell so in the range 4-7 inclusive
void LVL1::eFEXtauAlgo::setSupercellSeed()
{
  unsigned int seed = 4;
  int max_et = 0;
  int cell_et = 0;
  for(unsigned int i = 4; i < 8; i++)
  {
    cell_et = m_em2cells[i][1];
    if (cell_et > max_et)
    {
      seed = i;
      max_et = cell_et;
    }
  }
  m_seed = seed;
}

// Return the bitwise value of the given Et
// See eFEXtauBaseAlgo for a first attempt at this
unsigned int LVL1::eFEXtauAlgo::getBitwiseEt()
{
    unsigned int out = 0;
    return out;
}

bool LVL1::eFEXtauAlgo::getUnD()
{
    return m_und;
}

unsigned int LVL1::eFEXtauAlgo::getSeed()
{
    return m_seed;
}
