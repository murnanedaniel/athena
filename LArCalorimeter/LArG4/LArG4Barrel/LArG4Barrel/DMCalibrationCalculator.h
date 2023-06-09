/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// LArG4::DM::CalibrationCalculator
//   11-april-2006   G.Unal
//
// this class provides calibration hits for steps into the Dead Matter between
//  barrel and End-Cap cryostats

#ifndef LArG4_DM_CalibrationCalculator_H
#define LArG4_DM_CalibrationCalculator_H

#include "LArG4Code/VCalibrationCalculator.h"
#include "LArG4Code/LArG4Identifier.h"

#include "CaloG4Sim/SimulationEnergies.h"

#include "globals.hh"

#include <vector>

// Forward declaractions:
class G4Step;

// Note the use of nested namespaces (mainly to prevent long names
// like LArG4HECCalibrationCalculator).  This class is contained in
// the namespace LArG4::DM

namespace LArG4 {

  namespace DM {

    class CalibrationCalculator : public VCalibrationCalculator {
    public:
    
      CalibrationCalculator();
      virtual ~CalibrationCalculator();
    
      // The Process method returns a boolean value.  If it's true, the
      // hit can be used by Geant4; if it's false, there's something wrong
      // with the energy deposit and it should be ignored.

      // For calibration work, most of the time we want the calculator
      // to determine both the energy and the identifier.  However,
      // sometimes we want it calculate only the identifier (for
      // escaped energy), or only the energy (no known application
      // yet, but you can never tell).  Use the enum (defined in
      // VCalibrationCalculator.h) to control any special processing.

      virtual G4bool Process (const G4Step* step, 
			      const eCalculatorProcessing p = kEnergyAndID);
    
      // The cell identifier determined by the Process method.
      virtual const LArG4Identifier& identifier() const { return m_identifier; }
    
      // The calibration energies as determined by the Process method for
      // the current G4Step.  Units are the native G4 unit of energy.
      virtual const std::vector<G4double>& energies() const { return m_energies; }

    private:

      // The values calculated by Process().
      LArG4Identifier m_identifier;
      std::vector<G4double> m_energies;

      // Energy calculator
      CaloG4::SimulationEnergies m_energyCalculator;

    };

  } // namespace DM

} // namespace LArG4

#endif // LArG4_DM_CalibrationCalculator_H
