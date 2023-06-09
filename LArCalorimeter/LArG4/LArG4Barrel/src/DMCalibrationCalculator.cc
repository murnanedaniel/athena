/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// LArG4::DM:CalibrationCalculator
// Prepared 11-apr-2006 G.Unal       

// This class calculates the values needed for calibration hits in the
// simulation for steps in the Dead Matter between Barrel and End-Cap cryostats

// #define DEBUG_HITS

#include "LArG4Barrel/DMCalibrationCalculator.h"

#include "LArG4Code/LArG4Identifier.h"
#include "LArG4Code/VCalibrationCalculator.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <map>
#include <algorithm>
#include <cmath>
#include <climits>

namespace LArG4 {

  namespace DM {

     ///////////////////////////////////////////////////////////
    // Methods
    ///////////////////////////////////////////////////////////

    CalibrationCalculator::CalibrationCalculator()
    {
    }


    CalibrationCalculator::~CalibrationCalculator() 
    {
    }


    G4bool CalibrationCalculator::Process( const G4Step* a_step,
		    			      const eCalculatorProcessing a_process )
    {
      // Use the calculators to determine the energies and the
      // identifier associated with this G4Step.  Note that the
      // default is to process both the energy and the ID.

      m_energies.clear();
      if ( a_process == kEnergyAndID  ||  a_process == kOnlyEnergy )
	{
	  m_energyCalculator.Energies( a_step, m_energies );
	}
      else
	for (unsigned int i=0; i != 4; i++) m_energies.push_back( 0. );


      m_identifier.clear();
      if ( a_process == kEnergyAndID  ||  a_process == kOnlyID )
      {
// Calculate the identifier. 

// static ?
	  static const double oneOverDeta = 10.;       //   1/Deta = 1./0.1 = 10.
	  static const double oneOverDphi = 32./M_PI;  //   1/Dphi

// Calculate the mid-point of the step, and the simple geometry variables.

	  G4StepPoint* pre_step_point = a_step->GetPreStepPoint();
	  G4StepPoint* post_step_point = a_step->GetPostStepPoint();

	  G4ThreeVector startPoint = pre_step_point->GetPosition();
	  G4ThreeVector endPoint   = post_step_point->GetPosition();
	  G4ThreeVector p = (startPoint + endPoint) * 0.5;

	  //G4double rho = p.perp();
	  G4double eta = fabs( p.pseudoRapidity() );
	  G4double phi = p.phi(); 
	  if ( phi < 0. ) phi += 2.*M_PI; // Normalize for phiBin calculation

#ifdef DEBUG_HITS
          std::cout << " DMCalibrationCalculator rho,eta,phi " << p.perp() << " " << eta << " " << phi << std::endl;
#endif

	  // subdet = +/-4    "+" or " -" according to sign of Z in World coorinate
	  G4int subdet = ( p.z() > 0.) ? 4 : -4;
	  G4int phiBin = (int) ( phi * oneOverDphi );
          if (phiBin>63) phiBin=63;
          G4int type = 1;

	  // Initialize identifier variables with (invalid) default
	  // values (INT_MIN is defined in <climits>).
	  G4int sampling = INT_MIN;
	  G4int region = INT_MIN;
	  G4int etaBin = INT_MIN;

          if (eta < 1.0) {
              sampling = 2;
              region   = 1;
              etaBin = (int) ( eta * oneOverDeta );         
          }
          else if (eta < 1.5) {
              sampling = 2;
              region   = 2;
              etaBin = (int) ( (eta-1.) * oneOverDeta );
          }
          else if (eta < 1.6) {
              sampling = 1;
              region   = 4;
              etaBin   = 0;
          }
          else if (eta < 1.8) {
               sampling = 1;
               region   = 5;
               etaBin = (int) ( (eta-1.5) * oneOverDeta );
          }
          else if (eta < 3.2) {
               sampling = 1;
               region   = 6;
               etaBin = (int) ( (eta-1.3) * oneOverDeta );
          }
          else {
            std::cout << " DM::CalibrationCalculator  hit above 3.2 in eta !!! " << std::endl;
            return false; 
          }


// Append the cell ID to the (empty) identifier.
          m_identifier 	<< 10         // Calorimeter
   		    	<< subdet     // LAr +/-4 where "+" or " -" according to
	                 	         // the sign of Z in World coorinate
		   	<< type
			<< sampling
			<< region
			<< etaBin
			<< phiBin;

      }
  
#ifdef DEBUG_HITS
      G4double energy = accumulate(m_energies.begin(),m_energies.end(),0.);
      std::cout << "LArG4::DM::CalibrationCalculator::Process"
		<< " ID=" << std::string(m_identifier)
		<< " energy=" << energy
		<< " energies=(" << m_energies[0]
		<< "," << m_energies[1]
		<< "," << m_energies[2]
		<< "," << m_energies[3] << ")"
		<< std::endl;
#endif

      // Check for bad result.
      if ( m_identifier == LArG4Identifier() )
	return false;

      return true;
    }

  } // namespace DM

} // namespace LArG4
