/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
 *  Lepton flavour violating differential cross sections
 */

#include "Pythia8_i/UserProcessFactory.h"

// calls to fortran routines
#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <vector>

namespace Pythia8{
  class Sigma2qqbar2emu;
}

// This one line registers this process with the factory
// An instance of this process can be created with 
// UserProcessFactory::create("qqbar2emu")
Pythia8_UserProcess::UserProcessFactory::Creator<Pythia8::Sigma2qqbar2emu> qqbar2emuCreator("qqbar2emu");

namespace Pythia8 {
  
  /**
   *
   *  Sigma1ql2LFV class.
   *  Cross section for q qbar -> e mu (t:squark).
   */
  
  class Sigma2qqbar2emu : public Sigma2Process{
    
  public:
    //Constructor.
    Sigma2qqbar2emu()
    : lambdaCoup(0.), sigmaTemp(0.), m_stop(0.),
    lambdaSqddbar(0.), lambdaSqdsbar(0.), lambdaSqsdbar(0.), lambdaSqssbar(0.),
    Nc(0), myRand(0) 
    {}
    
    
    /// Initialize process. 
    void initProc() {
      
      // Store parameters.
      lambdaCoup = 0.2;
      Nc = 3;
      m_stop = particleDataPtr->m0(1000006);   //95.5;
      lambdaSqddbar = 0.065;
      lambdaSqdsbar = 3.8e-5;
      lambdaSqsdbar = 3.8e-5;
      lambdaSqssbar = 0.065;
      myRand = new Rndm(1234567);
      
      return;
    } 
    
    //*********
    
    /// Evaluate sigmaHat(sHat), part independent of incoming flavour. 
    void sigmaKin() {
      
      // The 2 out front is for e- mu+ plus e+ mu -
      sigmaTemp = 2./((float)Nc * 64. * M_PI);
      
      return;
    }
    
    //*********
    
    /// Evaluate sigmaHat(sHat), part dependent of incoming flavour. 
    double sigmaHat() {
      
      double lambdaSq = 0.;
      if((id1 == 1 && id2 == -1) || (id1 == -1 && id2 == 1)) lambdaSq = lambdaSqddbar;
      else if((id1 == 1 && id2 == -3) || (id1 == -3 && id2 == 1)) lambdaSq = lambdaSqdsbar;
      else if((id1 == 3 && id2 == -1) || (id1 == -1 && id2 == 3)) lambdaSq = lambdaSqsdbar;
      else if((id1 == 3 && id2 == -3) || (id1 == -3 && id2 == 3)) lambdaSq = lambdaSqssbar;
      double sigma = pow2(lambdaSq)*sigmaTemp*tH2/(sH2*pow2(tH-pow2(m_stop)));
      return sigma;
    }
    
    //*********
    
    // Select identity, colour and anticolour.
    
    void setIdColAcol() {
      
      // Flavours.
      double x = myRand->flat();
      if(x < 0.5) setId( id1, id2, 11, -13);
      else setId( id1, id2, -11, 13);
      
      // tH defined between f and LQ: must swap tHat <-> uHat if qbar q in.
      swapTU = (id1 < 0); 
      
      //cab I don't understand what's done here.
      // Colour flow topologies.
      if (id1 > 0) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
      else         setColAcol( 0, 1, 1, 0, 0, 0, 0, 0);
      
      return;
    }
    
    // Info on the subprocess.
    string name()    const {return "q qbar -> e mu (t:squark)";}
    //cab Need Codes for q qbar -> e+- mu-+
    int    code()    const {return 10001;}
    string inFlux()  const {return "qq";}
    //cab It looks like here is where the final state is specified
    //cab electrons are +-11, muons are +-13
    int    id3Mass() const {return 11;}
    int    id4Mass() const {return 13;}
    
  private:
    
    // Parameters set at initialization or for current kinematics.
    double lambdaCoup, sigmaTemp, m_stop;
    double lambdaSqddbar, lambdaSqdsbar, lambdaSqsdbar, lambdaSqssbar;
    int Nc;
    Rndm* myRand;
  
  };
}
