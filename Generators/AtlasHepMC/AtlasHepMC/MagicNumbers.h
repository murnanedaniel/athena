/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/* Author: Andrii Verbytskyi andrii.verbytskyi@mpp.mpg.de */

#ifndef ATLASHEPMC_MAGICNUMBERS_H
#define ATLASHEPMC_MAGICNUMBERS_H

#include <limits>

namespace HepMC {
/// @brief Constant defining the barcode threshold distinguishing generator record entries from detector sim ones
/// @todo The sim barcodes start at 1M in MC15, so we should update the 200k threshold,
///   but >= 200k is still a valid test for b = 1M so let's keep it this way until MC12 is long-dead.
constexpr int SIM_BARCODE_THRESHOLD = 200000;
constexpr int SIM_REGENERATION_INCREMENT = 1000000;

constexpr int PARTONPDGMAX = 43;
constexpr int NPPDGMIN = 1000000;
constexpr int NPPDGMAX = 8999999;
constexpr int PHOTOSMIN = 10000;

/// @brief This barcode is used by objects matched to particles from pile-up
/// interactions in standard MC Production
constexpr int crazyParticleBarcode(std::numeric_limits<int32_t>::max());

template <class T>  inline bool is_simulation_particle(T p){ return (barcode(p)>SIM_BARCODE_THRESHOLD);}
template <>  inline bool is_simulation_particle(int b){ return (b>SIM_BARCODE_THRESHOLD);}

template <class T>  inline bool is_simulation_vertex(T p){ return (barcode(p)<-SIM_BARCODE_THRESHOLD);}
template <>  inline bool is_simulation_vertex(int b){ return (b<-SIM_BARCODE_THRESHOLD);}

}
#endif
