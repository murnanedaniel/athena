/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TRT_G4Utilities/TRTParametersOfBarrelStraws.hh"
#include "TRT_G4Utilities/TRTParameters.hh"
#include "TRT_G4Utilities/TRTUtilities.hh"
#include "TRT_G4Utilities/TRTOutputFile.hh"
#include <fstream>



  // Called by TRTConstructionOfBarrelStraws::TRTConstructionOfBarrelStraws

TRTParametersOfBarrelStraws::TRTParametersOfBarrelStraws()
{
  m_pParameters = TRTParameters::GetPointer();

  DefineParameters();

  if (m_pParameters->GetInteger("PrintParametersOfBarrelStraws"))
    PrintParameters();
}


// Called by TRTConstructionOfBarrelStraws::~TRTConstructionOfBarrelStraws

TRTParametersOfBarrelStraws::~TRTParametersOfBarrelStraws()
{
}


// Called by TRTParametersOfBarrelStraws

void TRTParametersOfBarrelStraws::DefineParameters()
{
  // Parameters of straw hole:
  m_outerRadiusOfStrawHole = m_pParameters->GetDouble("OuterRadiusOfStrawHole");
  m_lengthOfStrawHole = m_pParameters->GetDouble("LengthOfBarrelVolume");

  // Parameters of straw:
  m_outerRadiusOfStraw = m_pParameters->GetDouble("OuterRadiusOfStraw");
  m_lengthOfStraw = m_lengthOfStrawHole;

  // Parameters of wire:
  m_outerRadiusOfWire = m_pParameters->GetDouble("OuterRadiusOfWire");
  m_lengthOfWire = m_lengthOfStraw;

  // Parameters of twister:
  m_innerRadiusOfTwister = m_outerRadiusOfWire;
  m_outerRadiusOfTwister = m_pParameters->GetDouble("InnerRadiusOfStraw");
  m_lengthOfTwister = m_pParameters->GetDouble("LengthOfTwister");

  // Parameters of dead regions:
  m_innerRadiusOfDeadRegion = m_outerRadiusOfWire;
  m_outerRadiusOfDeadRegion = m_outerRadiusOfTwister;
  m_lengthOfDeadRegion = m_pParameters->GetDouble("LengthOfDeadRegion");
  m_lengthOfLongDeadRegion =
    m_pParameters->GetDouble("LengthOfInsensitiveRegionA") - m_lengthOfTwister / 2.;

  m_positionOfDeadRegionLA = (m_lengthOfStraw - m_lengthOfDeadRegion) / 2.;
  m_positionOfDeadRegionLB = (m_lengthOfTwister + m_lengthOfDeadRegion) / 2.;
  m_positionOfDeadRegionSA = m_positionOfDeadRegionLA;
  m_positionOfLongDeadRegionSB = (m_lengthOfTwister + m_lengthOfLongDeadRegion) / 2.;

  // Parameters of gas:
  m_innerRadiusOfGas = m_innerRadiusOfDeadRegion;
  m_outerRadiusOfGas = m_outerRadiusOfDeadRegion;
  m_lengthOfGasL = (m_lengthOfStraw - m_lengthOfTwister) / 2. -
    m_lengthOfDeadRegion * 2.;
  m_lengthOfGasS = (m_lengthOfStraw - m_lengthOfTwister) / 2. - m_lengthOfDeadRegion -
    m_lengthOfLongDeadRegion;

  m_positionOfGasL = (m_lengthOfGasL + m_lengthOfTwister) / 2. + m_lengthOfDeadRegion;
  m_positionOfGasS = (m_lengthOfGasS + m_lengthOfTwister) / 2. +
    m_lengthOfLongDeadRegion;
}


// Called by TRTParametersOfBarrelStraws

void TRTParametersOfBarrelStraws::PrintParameters() const
{
  TRTOutputFile* pOutputFile = TRTOutputFile::GetPointer();

  std::ofstream& output = pOutputFile->GetReference();

  output << std::endl;
  output << "***** TRTParametersOfBarrelStraws::PrintParameters *****"
         << std::endl;

  output << "Parameters of barrel straw hole:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_outerRadiusOfStrawHole,
                                    m_lengthOfStrawHole);

  output << std::endl << "Parameters of barrel straw:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_outerRadiusOfStraw, m_lengthOfStraw);

  output << std::endl << "Parameters of barrel straw wire:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_outerRadiusOfWire, m_lengthOfWire);

  output << std::endl << "Parameters of barrel straw twister:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfTwister,
                                    m_outerRadiusOfTwister, m_lengthOfTwister);

  output << std::endl << "Parameters of gas S for short straw:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfGas, m_outerRadiusOfGas,
    m_lengthOfGasS, m_positionOfGasS);

  output << std::endl << "Parameters of dead region SA for short straw:"
         << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfDeadRegion,
    m_outerRadiusOfDeadRegion, m_lengthOfDeadRegion, m_positionOfDeadRegionSA);

  output << std::endl << "Parameters of long dead region SB for short straw:"
         << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfDeadRegion,
    m_outerRadiusOfDeadRegion, m_lengthOfLongDeadRegion,
    m_positionOfLongDeadRegionSB);

  output << std::endl << "Parameters of gas L for long straw:" << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfGas, m_outerRadiusOfGas,
    m_lengthOfGasL, m_positionOfGasL);

  output << std::endl << "Parameters of dead region LA for long straw:"
         << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfDeadRegion,
    m_outerRadiusOfDeadRegion, m_lengthOfDeadRegion, m_positionOfDeadRegionLA);

  output << std::endl << "Parameters of dead region LB for long straw:"
         << std::endl;
  TRTUtilities::PrintTubeParameters(*pOutputFile,
                                    m_innerRadiusOfDeadRegion,
    m_outerRadiusOfDeadRegion, m_lengthOfDeadRegion, m_positionOfDeadRegionLB);

  output << std::endl;
}
