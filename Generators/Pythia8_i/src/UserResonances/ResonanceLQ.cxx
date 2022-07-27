/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// The ResonanceLQ class modifies the default Pythia8 setup 
// to allow leptoquarks to decay to top quarks
// This is a copy of existing pythia ResonanceLeptoquark class
//
// author Katharine Leney
// based on ResonanceExcitedCI (Olya Igonkina and James Monk)

#include "Pythia8_i/UserResonanceFactory.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8{
  class ResonanceLQ;
}

Pythia8_UserResonance::UserResonanceFactory::Creator<Pythia8::ResonanceLQ> resonanceLeptoquarkCreator("LQ");

namespace Pythia8{
  
  class ResonanceLQ : public ResonanceWidths {
    
  public:
    
    // Constructor. 
    ResonanceLQ(int idResIn): 
      m_kCoup(0.)
    {    
      initBasic(idResIn);
      std::cout << " ResonanceLQ constructor\n"; 
    } 
  
  private:
    
    void initConstants() {
      
      // Locally stored properties and couplings.
      m_kCoup = settingsPtr->parm("LeptoQuark:kCoup");
      
            // make a copy of shared pointer before usage (starting Py8.307 particlePtr is of a type std::weak_ptr<Pythia8::ParticleDataEntry>
#if PYTHIA_VERSION_INTEGER >= 8307      
      ParticleDataEntryPtr particleSPtr = particlePtr.lock();
#else
      ParticleDataEntry* particleSPtr = particlePtr;
#endif

      // Check that flavour info in decay channel is correctly set.
      int id1Now = particleSPtr->channel(0).product(0);
      int id2Now = particleSPtr->channel(0).product(1);
      
      // ============================================================
      // Modify standard Pythia8 setup to allow decays to top quarks
      // ============================================================
      
      if (id1Now < 1 || id1Now > 6) {
	std::cout << "ERROR in ResonanceLQ::init: unallowed input quark flavour reset to u" << std::endl; 
	id1Now   = 2;
	particleSPtr->channel(0).product(0, id1Now);
      }
      if (std::abs(id2Now) < 11 || std::abs(id2Now) > 16) {
	std::cout << "ERROR in ResonanceLQ::init:unallowed input lepton flavour reset to e-" << std::endl; 
	id2Now   = 11;
	particleSPtr->channel(0).product(1, id2Now);
      }
      
      // Set/overwrite charge and name of particle.
      bool changed  = particleSPtr->hasChanged();
      
      int chargeLQ  = particleDataPtr->chargeType(id1Now) 
	+ particleDataPtr->chargeType(id2Now);
      
      particleSPtr->setChargeType(chargeLQ); 
      
      std::string nameLQ = "LQ_" + particleDataPtr->name(id1Now) + ","
	+ particleDataPtr->name(id2Now);
      
      particleSPtr->setNames(nameLQ, nameLQ + "bar"); 
      if (!changed) particleSPtr->setHasChanged(false);
      
      return;
    }
    
    //--------------------------------------------------------------------------
    
    // Calculate various common prefactors for the current mass.
    
    void calcPreFac(bool) {

#ifdef PYTHIA_VERSION_INTEGER
  #if PYTHIA_VERSION_INTEGER > 8300
      CoupSM* couplingsPtr = infoPtr->coupSMPtr;
  #endif
#endif
      alpEM   = couplingsPtr->alphaEM(mHat * mHat);
      preFac  = 0.25 * alpEM * m_kCoup * mHat; 
      
      return;
    }
    
    //--------------------------------------------------------------------------
    
    // Calculate width for currently considered channel.
    
    void calcWidth(bool) {
      
      // Check that above threshold.
      if (ps == 0.) return;
      
      // Width into lepton plus quark.
      if (id1Abs > 10 && id1Abs < 17 && id2Abs < 7) widNow = preFac * pow3(ps);
      return;
    }
    
    // Locally stored properties and couplings.
    //double m_lambda, m_coupF, m_coupFprime, m_coupFcol, m_sin2tW, m_cos2tW;
    double m_kCoup;
  };
}
