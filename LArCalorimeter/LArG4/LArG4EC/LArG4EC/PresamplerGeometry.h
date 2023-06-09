/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// LArG4::EC::PresamplerGeometry

// This class contains the geometry calculations needed to calculate
// an identifier for a given G4Step.

// 17-Aug-2004: Mikhail Leltchouk wrote the original code, Bill
// Seligman revised it for current use.

// Note: This class is intended for use in calculating identifiers in
// both active and inactive regions of the detector.  For calibration
// studies, it must work properly whether the energy is deposited in
// the liquid argon or the absorber.

#ifndef LArG4_EC_PresamplerGeometry_H
#define LArG4_EC_PresamplerGeometry_H

#include "G4ThreeVector.hh"
#include "globals.hh"

// Forward declarations.
class LArG4Identifier;
class G4Step;

// Note the use of nested namespaces (mainly to prevent long names
// like LArG4HECGeometry).  This class is contained in the namespace
// LArG4::EC.

namespace LArG4 {

  namespace EC {

    class PresamplerGeometry {

    public:
      // Standard implementation of a singleton pattern.
      static const PresamplerGeometry* GetInstance();
      virtual ~PresamplerGeometry();

      // 15-Jan-2002 WGS: A "lookup" function for detector measurements,
      // sizes, and other values.
      enum kValue {
	rMinEndcapPresampler,
	rMaxEndcapPresampler,
	zEndcapPresamplerFrontFace,
	zEndcapPresamplerBackFace,
	EndcapPresamplerHalfThickness,
	EndcapPresamplerZpositionInMother
      };
      G4double GetValue(const kValue) const;
      
      // This is the "meat" of this class: calculate the identifier
      // given a G4Step.
      LArG4Identifier CalculateIdentifier( const G4Step* ) const;

    protected:
      // Constructor is protected according to the singleton pattern.
      PresamplerGeometry();

    private:


      struct Clockwork;
      Clockwork *c;

      PresamplerGeometry (const PresamplerGeometry&);
      PresamplerGeometry& operator= (const PresamplerGeometry&);
    };

  } // namespace EC

} // namespace LArG4

#endif // LArG4_EC_PresamplerGeometry_H
