/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

//***************************************************************************
//                                  eFEXFormTOBs
//                              -------------------
//     begin                : 30 04 2021
//     email                : nicholas.andrew.luongo@cern.ch
//  ***************************************************************************/

#include "L1CaloFEXSim/eFEXFormTOBs.h"

namespace LVL1 {

  // default constructor for persistency

eFEXFormTOBs::eFEXFormTOBs(const std::string& type, const std::string& name, const IInterface* parent):
  AthAlgTool(type, name, parent)
  {
    declareInterface<IeFEXFormTOBs>(this);
  }

/** Desctructor */
eFEXFormTOBs::~eFEXFormTOBs()
{
}

StatusCode eFEXFormTOBs::initialize()
{
  return StatusCode::SUCCESS;
}

uint32_t eFEXFormTOBs::formTauTOBWord(int & fpga, int & eta, int & phi, unsigned int & et, unsigned int & iso, unsigned int & seed, unsigned int & und, unsigned int & ptMinTopo)
{

  uint32_t tobWord = 0;

  //rescale from 25 MeV eFEX steps to 100 MeV for the TOB
  unsigned int etTob = 0;
  etTob = et*m_eFexStep/m_eFexTobStep; 

  // Truncate at 12 bits, set to max value of 4095, 0xfff, or 111111111111
  if (etTob > 0xfff) etTob = 0xfff;

  // Create bare minimum tob word with et, eta, phi, and fpga index, bitshifted to the appropriate locations
  tobWord = tobWord + (fpga << 30) + (eta << 27) + (phi << 24) + (iso << 18) + (seed << 16) + (und << 15) + etTob;

  ATH_MSG_DEBUG("Tau tobword: " << std::bitset<32>(tobWord) );

  unsigned int minEtThreshold = ptMinTopo;
  if (etTob < minEtThreshold) return 0;
  else return tobWord;
}

uint32_t eFEXFormTOBs::formEmTOBWord(int & fpga, int & eta, int & phi, unsigned int & rhad, unsigned int & wstot, unsigned int & reta, unsigned int & seed, unsigned int & et, unsigned int & ptMinTopo)
{
  uint32_t tobWord = 0;

  unsigned int etTob = 0;
  etTob = et*m_eFexStep/m_eFexTobStep; //rescale from 25 MeV eFEX steps to 100 MeV for the TOB

  // Truncate at 12 bits, set to max value of 4095, 0xfff, or 111111111111
  if (etTob > 0xfff) etTob = 0xfff;

  // Create bare minimum tob word with et, eta, phi, and fpga index, bitshifted to the appropriate locations
  tobWord = tobWord + (fpga << 30) + (eta << 27) + (phi << 24) + (rhad << 22) + (wstot << 20) + (reta << 18) + (seed << 16) + etTob;

  ATH_MSG_DEBUG("EM tobword: " << std::bitset<32>(tobWord) );

  // Some arbitrary cut so that we're not flooded with tobs, to be taken from the Trigger menu in the future!
  unsigned int minEtThreshold = ptMinTopo;
  if (etTob < minEtThreshold) return 0;
  else return tobWord;
}

} // end of namespace bracket
