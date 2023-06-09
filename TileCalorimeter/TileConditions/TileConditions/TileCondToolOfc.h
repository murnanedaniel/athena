/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TILECONDITIONS_TILECONDTOOLOFC_H
#define TILECONDITIONS_TILECONDTOOLOFC_H

// Gaudi includes
#include "GaudiKernel/ToolHandle.h"

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"

// Tile includes
#include "TileConditions/ITileCondToolOfc.h"
#include "TileConditions/TileInfo.h"
// #include "TileConditions/TileCondToolNoiseSample.h"
#include "TileConditions/TileCondToolAutoCr.h"
#include "TileConditions/TileCondToolPulseShape.h"

#include "CLHEP/Matrix/Matrix.h"
#include <CLHEP/Matrix/Vector.h>
#include <map>
#include <vector>
#include <string>
#include <TMatrixD.h>
#include <iostream>
#include <iomanip>

#define MAX_ADCS 13248

class CaloDetDescrManager;
class IdContext;

/**
 *
 * @class TileCondToolOfc
 * @brief Calculates OFCs on the fly using pulse shapes and A/C matrix from database
 *
 * Optionally, it can create cache table of OFCs with 1-ns step to minimize
 * CPU time. Also, by request, unity A/C matrix can be used.
 */
class TileCondToolOfc: public AthAlgTool, public ITileCondToolOfc {
  public:

    TileCondToolOfc(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TileCondToolOfc();

    StatusCode initialize();
    StatusCode finalize();
    static const InterfaceID& interfaceID();

    void CalcWeights(unsigned int drawerIdx, unsigned int channel, int gain, float phase, bool of2);

    //===============================================================
    //=== ITileCondTollOfc methods
    //===============================================================

    const TileOfcWeightsStruct* getOfcWeights(unsigned int drawerIdx, unsigned int channel, unsigned int adc, float& phase, bool of2);

    int getNSamples(void);

  private:

    typedef std::map<int, TileOfcWeightsStruct*> ADCPhaseCache;
    ADCPhaseCache* m_ofc_phase_cache[MAX_ADCS * 2];

    ToolHandle<TileCondToolPulseShape> m_tileToolPulseShape;
    ToolHandle<TileCondToolAutoCr> m_tileToolAutoCr;
    //  ToolHandle<TileCondToolNoiseSample> m_tileToolNoiseSample;
    const TileInfo* m_tileInfo;
    TileOfcWeightsStruct* m_weights;

    //=== properties
    int m_nSamples;
    int m_t0Sample;
    bool m_deltaCorrelation;
    unsigned int m_cache;
};

#endif
