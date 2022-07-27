/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


#ifndef TRTParametersOfBarrelStraws_hh
#define TRTParametersOfBarrelStraws_hh

#include "globals.hh"

class TRTParameters;


class TRTParametersOfBarrelStraws
{
  friend class TRTConstructionOfBarrelStraws;

  public:
    TRTParametersOfBarrelStraws();
    ~TRTParametersOfBarrelStraws();

  private:
    TRTParametersOfBarrelStraws (const TRTParametersOfBarrelStraws&); 
    TRTParametersOfBarrelStraws& operator= (const TRTParametersOfBarrelStraws&); 
    void DefineParameters();
    void PrintParameters() const;

    double m_outerRadiusOfStrawHole = 0.0;
    double m_lengthOfStrawHole = 0.0;

    double m_outerRadiusOfStraw = 0.0;
    double m_lengthOfStraw = 0.0;

    double m_innerRadiusOfGas = 0.0;
    double m_outerRadiusOfGas = 0.0;
    double m_lengthOfGasL = 0.0;
    double m_lengthOfGasS = 0.0;

    double m_positionOfGasL = 0.0;
    double m_positionOfGasS = 0.0;

    double m_innerRadiusOfDeadRegion = 0.0;
    double m_outerRadiusOfDeadRegion = 0.0;
    double m_lengthOfDeadRegion = 0.0;
    double m_lengthOfLongDeadRegion = 0.0;

    double m_positionOfDeadRegionLA = 0.0;
    double m_positionOfDeadRegionLB = 0.0;
    double m_positionOfDeadRegionSA = 0.0;
    double m_positionOfLongDeadRegionSB = 0.0;

    double m_innerRadiusOfTwister = 0.0;
    double m_outerRadiusOfTwister = 0.0;
    double m_lengthOfTwister = 0.0;

    double m_outerRadiusOfWire = 0.0;
    double m_lengthOfWire = 0.0;

    const TRTParameters* m_pParameters;
};

#endif
