/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "UserHooksUtils.h"
#include "UserSetting.h"
#include "Pythia8_i/UserHooksFactory.h"
#include "boost/lexical_cast.hpp"
#include <stdexcept>
#include <iostream>

namespace Pythia8{
  class main31;
}

Pythia8_UserHooks::UserHooksFactory::Creator<Pythia8::main31> main31Creator("Main31");

namespace Pythia8{
  
  // Use userhooks to veto PYTHIA emissions above the POWHEG scale.
  
  class main31 : public UserHooks {
    
  public:  
    
    // Constructor and destructor.
    main31() : m_nFinal("Main31:NFinal", 2),
    m_pTHardMode("Main31:pTHard", 2),
    m_pTDefMode("Main31:pTdef", 1),
    m_vetoMode(1), m_vetoCount(3),
    m_pTemtMode(0),
    m_emittedMode(0),
    m_MPIvetoMode(0),
    m_pThard(0), m_pTMPI(0),
    m_accepted(0),
    m_nAcceptSeq(0),
    m_nISRveto(0), m_nFSRveto(0) {};
    ~main31() {}
    
    //--------------------------------------------------------------------------
    
    // Routines to calculate the pT (according to pTdefMode) in a splitting:
    //   ISR: i (radiator after)  -> j (emitted after) k (radiator before)
    //   FSR: i (radiator before) -> j (emitted after) k (radiator after)
    // For the Pythia pT definition, a recoiler (after) must be specified.
    
    // Compute the Pythia pT separation. Based on pTLund function in History.cc
    double pTpythia(const Event &e, int RadAfterBranch, int EmtAfterBranch,
                    int RecAfterBranch, bool FSR) {
      
      // Convenient shorthands for later
      Vec4 radVec = e[RadAfterBranch].p();
      Vec4 emtVec = e[EmtAfterBranch].p();
      Vec4 recVec = e[RecAfterBranch].p();
      int  radID  = e[RadAfterBranch].id();
      
      // Calculate virtuality of splitting
      double sign = (FSR) ? 1. : -1.;
      Vec4 Q(radVec + sign * emtVec); 
      double Qsq = sign * Q.m2Calc();
      
      // Mass term of radiator
      double m2Rad = (abs(radID) >= 4 && abs(radID) < 7) ?
      pow2(particleDataPtr->m0(radID)) : 0.;
      
      // z values for FSR and ISR
      double z, pTnow;
      if (FSR) {
        // Construct 2 -> 3 variables
        Vec4 sum = radVec + recVec + emtVec;
        double m2Dip = sum.m2Calc();
        double x1 = 2. * (sum * radVec) / m2Dip;
        double x3 = 2. * (sum * emtVec) / m2Dip;
        z     = x1 / (x1 + x3);
        pTnow = z * (1. - z);
        
      } else {
        // Construct dipoles before/after splitting
        Vec4 qBR(radVec - emtVec + recVec);
        Vec4 qAR(radVec + recVec);
        z     = qBR.m2Calc() / qAR.m2Calc();
        pTnow = (1. - z);
      }
      
      // Virtuality with correct sign
      pTnow *= (Qsq - sign * m2Rad);
      
      // Can get negative pT for massive splittings
      if (pTnow < 0.) {
        cout << "Warning: pTpythia was negative" << endl;
        return -1.;
      }
      
#ifdef DBGOUTPUT
      cout << "pTpythia: rad = " << RadAfterBranch << ", emt = "
      << EmtAfterBranch << ", rec = " << RecAfterBranch
      << ", pTnow = " << sqrt(pTnow) << endl;
#endif
      
      // Return pT
      return sqrt(pTnow);
    }
    
    // Compute the POWHEG pT separation between i and j
    double pTpowheg(const Event &e, int i, int j, bool FSR) {
      
      // pT value for FSR and ISR
      double pTnow = 0.;
      if (FSR) {
        // POWHEG d_ij (in CM frame). Note that the incoming beams have not
        // been updated in the parton systems pointer yet (i.e. prior to any
        // potential recoil).
        int iInA = partonSystemsPtr->getInA(0);
        int iInB = partonSystemsPtr->getInB(0);
        double betaZ = - ( e[iInA].pz() + e[iInB].pz() ) /
        ( e[iInA].e()  + e[iInB].e()  );
                
        Vec4 iVecBst(e[i].p()), jVecBst(e[j].p());
        iVecBst.bst(0., 0., betaZ);
        jVecBst.bst(0., 0., betaZ);
        pTnow = sqrt( (iVecBst + jVecBst).m2Calc() *
                     iVecBst.e() * jVecBst.e() /
                     pow2(iVecBst.e() + jVecBst.e()) );
        
      } else {
        // POWHEG pT_ISR is just kinematic pT
        pTnow = e[j].pT();
      }
      
      // Check result
      if (pTnow < 0.) {
        cout << "Warning: pTpowheg was negative" << endl;
        return -1.;
      }
      
#ifdef DBGOUTPUT
      cout << "pTpowheg: i = " << i << ", j = " << j
      << ", pTnow = " << pTnow << endl;
#endif
      
      return pTnow;
    }
    
    // Calculate pT for a splitting based on pTdefMode.
    // If j is -1, all final-state partons are tried.
    // If i, k, r and xSR are -1, then all incoming and outgoing 
    // partons are tried.
    // xSR set to 0 means ISR, while xSR set to 1 means FSR
    double pTcalc(const Event &e, int i, int j, int k, int r, int xSRin) {
      
      // Loop over ISR and FSR if necessary
      double pTemt = -1., pTnow;
      int xSR1 = (xSRin == -1) ? 0 : xSRin;
      int xSR2 = (xSRin == -1) ? 2 : xSRin + 1;
      for (int xSR = xSR1; xSR < xSR2; xSR++) {
        // FSR flag
        bool FSR = (xSR == 0) ? false : true;
        
        // If all necessary arguments have been given, then directly calculate.
        // POWHEG ISR and FSR, need i and j.
        if ((m_pTDefMode(settingsPtr) == 0 || m_pTDefMode(settingsPtr) == 1) && i > 0 && j > 0) {
          pTemt = pTpowheg(e, i, j, (m_pTDefMode(settingsPtr) == 0) ? false : FSR);
          
          // Pythia ISR, need i, j and r.
        } else if (!FSR && m_pTDefMode(settingsPtr) == 2 && i > 0 && j > 0 && r > 0) {
          pTemt = pTpythia(e, i, j, r, FSR);
          
          // Pythia FSR, need k, j and r.
        } else if (FSR && m_pTDefMode(settingsPtr) == 2 && j > 0 && k > 0 && r > 0) {
          pTemt = pTpythia(e, k, j, r, FSR);
          
          // Otherwise need to try all possible combintations.
        } else {
          // Start by finding incoming legs to the hard system after
          // branching (radiator after branching, i for ISR).
          // Use partonSystemsPtr to find incoming just prior to the
          // branching and track mothers.
          int iInA = partonSystemsPtr->getInA(0);
          int iInB = partonSystemsPtr->getInB(0);
          while (e[iInA].mother1() != 1) { iInA = e[iInA].mother1(); }
          while (e[iInB].mother1() != 2) { iInB = e[iInB].mother1(); }
          
          // If we do not have j, then try all final-state partons
          int jNow = (j > 0) ? j : 0;
          int jMax = (j > 0) ? j + 1 : e.size();
          for (; jNow < jMax; jNow++) {
            
            // Final-state and coloured jNow only
            if (!e[jNow].isFinal() || e[jNow].colType() == 0) continue;
            
            // POWHEG
            if (m_pTDefMode(settingsPtr) == 0 || m_pTDefMode(settingsPtr) == 1) {
              
              // ISR - only done once as just kinematical pT
              if (!FSR) {
                pTnow = pTpowheg(e, iInA, jNow, (m_pTDefMode(settingsPtr) == 0) ? false : FSR);
                if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);
                
                // FSR - try all outgoing partons from system before branching 
                // as i. Note that for the hard system, there is no 
                // "before branching" information.
              } else {
                
                int outSize = partonSystemsPtr->sizeOut(0);
                for (int iMem = 0; iMem < outSize; iMem++) {
                  int iNow = partonSystemsPtr->getOut(0, iMem);
                  
                  // Coloured only, i != jNow and no carbon copies
                  if (iNow == jNow || e[iNow].colType() == 0) continue;
                  if (jNow == e[iNow].daughter1() 
                      && jNow == e[iNow].daughter2()) continue;
                  
                  pTnow = pTpowheg(e, iNow, jNow, (m_pTDefMode(settingsPtr) == 0) 
                                   ? false : FSR);
                  if (pTnow > 0.) pTemt = (pTemt < 0) 
                    ? pTnow : min(pTemt, pTnow);
                } // for (iMem)
                
              } // if (!FSR)
              
              // Pythia
            } else if (m_pTDefMode(settingsPtr) == 2) {
              
              // ISR - other incoming as recoiler
              if (!FSR) {
                pTnow = pTpythia(e, iInA, jNow, iInB, FSR);
                if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);
                pTnow = pTpythia(e, iInB, jNow, iInA, FSR);
                if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);
                
                // FSR - try all final-state coloured partons as radiator
                //       after emission (k).
              } else {
                for (int kNow = 0; kNow < e.size(); kNow++) {
                  if (kNow == jNow || !e[kNow].isFinal() ||
                      e[kNow].colType() == 0) continue;
                  
                  // For this kNow, need to have a recoiler.
                  // Try two incoming.
                  pTnow = pTpythia(e, kNow, jNow, iInA, FSR);
                  if (pTnow > 0.) pTemt = (pTemt < 0) 
                    ? pTnow : min(pTemt, pTnow);
                  pTnow = pTpythia(e, kNow, jNow, iInB, FSR);
                  if (pTnow > 0.) pTemt = (pTemt < 0) 
                    ? pTnow : min(pTemt, pTnow);
                  
                  // Try all other outgoing.
                  for (int rNow = 0; rNow < e.size(); rNow++) {
                    if (rNow == kNow || rNow == jNow ||
                        !e[rNow].isFinal() || e[rNow].colType() == 0) continue;
                    pTnow = pTpythia(e, kNow, jNow, rNow, FSR);
                    if (pTnow > 0.) pTemt = (pTemt < 0) 
                      ? pTnow : min(pTemt, pTnow);
                  } // for (rNow)
                  
                } // for (kNow)
              } // if (!FSR)
            } // if (pTdefMode)
          } // for (j)
        }
      } // for (xSR)
      
#ifdef DBGOUTPUT
      cout << "pTcalc: i = " << i << ", j = " << j << ", k = " << k
      << ", r = " << r << ", xSR = " << xSRin
      << ", pTemt = " << pTemt << endl;
#endif
      
      return pTemt;
    }
    
    //--------------------------------------------------------------------------
    
    // Extraction of m_pThard based on the incoming event.
    // Assume that all the final-state particles are in a continuous block
    // at the end of the event and the final entry is the POWHEG emission.
    // If there is no POWHEG emission, then m_pThard is set to Qfac.
    
    bool canVetoMPIStep()    { return true; }
    int  numberVetoMPIStep() { return 1; }
    bool doVetoMPIStep(int nMPI, const Event &e) {
      // Extra check on nMPI
      if (nMPI > 1) return false;
      
      // Find if there is a POWHEG emission. Go backwards through the
      // event record until there is a non-final particle. Also sum pT and
      // find pT_1 for possible MPI vetoing
      int    count = 0;
      double pT1 = 0., pTsum = 0.;
      for (int i = e.size() - 1; i > 0; i--) {
        if (e[i].isFinal()) {
          count++;
          pT1    = e[i].pT();
          pTsum += e[i].pT();
        } else break;
      }
      // Extra check that we have the correct final state
      if (count != m_nFinal(settingsPtr) && count != m_nFinal(settingsPtr) + 1) {
        cout << "Error: wrong number of final state particles in event" << endl;
        exit(1);
      }
      // Flag if POWHEG radiation present and index
      bool isEmt = (count == m_nFinal(settingsPtr)) ? false : true;
      int  iEmt  = (isEmt) ? e.size() - 1 : -1;
      
      // If there is no radiation or if pThardMode is 0 then set m_pThard to QRen.
      if (!isEmt || m_pTHardMode(settingsPtr) == 0) {
        m_pThard = infoPtr->QRen();
        
        // If pThardMode is 1 then the pT of the POWHEG emission is checked against
        // all other incoming and outgoing partons, with the minimal value taken
      } else if (m_pTHardMode(settingsPtr) == 1) {
        m_pThard = pTcalc(e, -1, iEmt, -1, -1, -1);
        
        // If pThardMode is 2, then the pT of all final-state partons is checked
        // against all other incoming and outgoing partons, with the minimal value
        // taken
      } else if (m_pTHardMode(settingsPtr) == 2) {
        m_pThard = pTcalc(e, -1, -1, -1, -1, -1);
        
      }
      
      // Find MPI veto pT if necessary
      if (m_MPIvetoMode == 1) {
        m_pTMPI = (isEmt) ? pTsum / 2. : pT1;
      }
      
#ifdef DBGOUTPUT
      cout << "doVetoMPIStep: QRen = " << infoPtr->QRen()
      << ", m_pThard = " << m_pThard << endl << endl;
#endif
      
//      std::cout<<"vetoScale = "<<m_pThard<<std::endl;
      
      // Initialise other variables
      m_accepted   = false;
      m_nAcceptSeq = m_nISRveto = m_nFSRveto = 0;
      
      // Do not veto the event
      return false;
    }
    
    //--------------------------------------------------------------------------
    
    // ISR veto
    
    bool canVetoISREmission() { return (m_vetoMode == 0) ? false : true; }
    bool doVetoISREmission(int, const Event &e, int iSys) {
      // Must be radiation from the hard system
      if (iSys != 0) return false;
      
      // If we already have m_accepted 'm_vetoCount' emissions in a row, do nothing
      if (m_vetoMode == 1 && m_nAcceptSeq >= m_vetoCount) return false;
      
      // Pythia radiator after, emitted and recoiler after.
      int iRadAft = -1, iEmt = -1, iRecAft = -1;
      for (int i = e.size() - 1; i > 0; i--) {
        if      (iRadAft == -1 && e[i].status() == -41) iRadAft = i;
        else if (iEmt    == -1 && e[i].status() ==  43) iEmt    = i;
        else if (iRecAft == -1 && e[i].status() == -42) iRecAft = i;
        if (iRadAft != -1 && iEmt != -1 && iRecAft != -1) break;
      }
      if (iRadAft == -1 || iEmt == -1 || iRecAft == -1) {
        e.list();
        cout << "Error: couldn't find Pythia ISR emission" << endl;
        exit(1);
      }
      
      // m_pTemtMode == 0: pT of emitted w.r.t. radiator
      // m_pTemtMode == 1: min(pT of emitted w.r.t. all incoming/outgoing)
      // m_pTemtMode == 2: min(pT of all outgoing w.r.t. all incoming/outgoing)
      int xSR      = (m_pTemtMode == 0) ? 0       : -1;
      int i        = (m_pTemtMode == 0) ? iRadAft : -1;
      int j        = (m_pTemtMode != 2) ? iEmt    : -1;
      int k        = -1;
      int r        = (m_pTemtMode == 0) ? iRecAft : -1;
      double pTemt = pTcalc(e, i, j, k, r, xSR);
      
#ifdef DBGOUTPUT
      cout << "doVetoISREmission: pTemt = " << pTemt << endl << endl;
#endif
      
      // Veto if pTemt > m_pThard
      if (pTemt > m_pThard) {
        m_nAcceptSeq = 0;
        m_nISRveto++;
        return true;
      }
      
      // Else mark that an emission has been m_accepted and continue
      m_nAcceptSeq++;
      m_accepted = true;
      return false;
    }
    
    //--------------------------------------------------------------------------
    
    // FSR veto
    
    bool canVetoFSREmission() { return (m_vetoMode == 0) ? false : true; }
    bool doVetoFSREmission(int, const Event &e, int iSys, bool) {
      // Must be radiation from the hard system
      if (iSys != 0) return false;
      
      // If we already have m_accepted 'm_vetoCount' emissions in a row, do nothing
      if (m_vetoMode == 1 && m_nAcceptSeq >= m_vetoCount) return false;
      
      // Pythia radiator (before and after), emitted and recoiler (after)
      int iRecAft = e.size() - 1;
      int iEmt    = e.size() - 2;
      int iRadAft = e.size() - 3;
      int iRadBef = e[iEmt].mother1();
      if ( (e[iRecAft].status() != 52 && e[iRecAft].status() != -53) ||
          e[iEmt].status() != 51 || e[iRadAft].status() != 51) {
        e.list();
        cout << "Error: couldn't find Pythia FSR emission" << endl;
        exit(1);
      }
      
      // Behaviour based on m_pTemtMode:
      //  0 - pT of emitted w.r.t. radiator before
      //  1 - min(pT of emitted w.r.t. all incoming/outgoing)
      //  2 - min(pT of all outgoing w.r.t. all incoming/outgoing)
      int xSR = (m_pTemtMode == 0) ? 1       : -1;
      int i   = (m_pTemtMode == 0) ? iRadBef : -1;
      int k   = (m_pTemtMode == 0) ? iRadAft : -1;
      int r   = (m_pTemtMode == 0) ? iRecAft : -1;
      
      // When m_pTemtMode is 0 or 1, iEmt has been selected
      double pTemt = 0.;
      if (m_pTemtMode == 0 || m_pTemtMode == 1) {
        // Which parton is emitted, based on m_emittedMode:
        //  0 - Pythia definition of emitted
        //  1 - Pythia definition of radiated after emission
        //  2 - Random selection of emitted or radiated after emission
        //  3 - Try both emitted and radiated after emission
        int j = iRadAft;
        if (m_emittedMode == 0 || (m_emittedMode == 2 && rndmPtr->flat() < 0.5)) j++;
        
        for (int jLoop = 0; jLoop < 2; jLoop++) {
          if      (jLoop == 0) pTemt = pTcalc(e, i, j, k, r, xSR);
          else if (jLoop == 1) pTemt = min(pTemt, pTcalc(e, i, j, k, r, xSR));
          
          // For m_emittedMode == 3, have tried iRadAft, now try iEmt
          if (m_emittedMode != 3) break;
          if (k != -1) swap(j, k); else j = iEmt;
        }
        
        // If m_pTemtMode is 2, then try all final-state partons as emitted
      } else if (m_pTemtMode == 2) {
        pTemt = pTcalc(e, i, -1, k, r, xSR);
        
      }
      
#ifdef DBGOUTPUT
      cout << "doVetoFSREmission: pTemt = " << pTemt << endl << endl;
#endif
      
      // Veto if pTemt > m_pThard
      if (pTemt > m_pThard) {
        m_nAcceptSeq = 0;
        m_nFSRveto++;
        return true;
      }
      
      // Else mark that an emission has been m_accepted and continue
      m_nAcceptSeq++;
      m_accepted = true;
      return false;
    }
    
    //--------------------------------------------------------------------------
    
    // MPI veto
    
    bool canVetoMPIEmission() { return (m_MPIvetoMode == 0) ? false : true; }
    bool doVetoMPIEmission(int, const Event &e) {
      if (m_MPIvetoMode == 1) {
        if (e[e.size() - 1].pT() > m_pTMPI) {
#ifdef DBGOUTPUT
          cout << "doVetoMPIEmission: pTnow = " << e[e.size() - 1].pT()
          << ", m_pTMPI = " << m_pTMPI << endl << endl;
#endif
          return true;
        }
      }
      return false;
    }
    
    //--------------------------------------------------------------------------
    
    // Functions to return information
    
    int    getNISRveto() { return m_nISRveto; }
    int    getNFSRveto() { return m_nFSRveto; }
    
  private:
    
    Pythia8_UserHooks::UserSetting<int> m_nFinal;
    Pythia8_UserHooks::UserSetting<int> m_pTHardMode;
    Pythia8_UserHooks::UserSetting<int> m_pTDefMode;
    
    int  m_vetoMode, m_vetoCount, m_pTemtMode,
    m_emittedMode, m_MPIvetoMode;
    double m_pThard, m_pTMPI;
    bool   m_accepted;
    // The number of m_accepted emissions (in a row)
    int m_nAcceptSeq;
    // Statistics on vetos
    unsigned long int m_nISRveto, m_nFSRveto;
    
  };
}
