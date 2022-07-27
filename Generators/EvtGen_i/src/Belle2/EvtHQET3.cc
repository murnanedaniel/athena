/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2016 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Florian Bernlochner                                      *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <assert.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include <string>

//#include "generators/evtgen/EvtGenModelRegister.h"
#include "EvtGen_i/EvtGenExternal/Belle2/EvtHQET3.h"
#include "EvtGen_i/EvtGenExternal/Belle2/EvtHQET3FF.h"

//B2_EVTGEN_REGISTER_MODEL(EvtHQET3);
//registerModel(new EvtHQET3);

using std::endl;

EvtHQET3::EvtHQET3():
  m_hqetffmodel(0)
  , m_calcamp(0)
{}

EvtHQET3::~EvtHQET3()
{
  delete m_hqetffmodel;
  m_hqetffmodel = 0;
  delete m_calcamp;
  m_calcamp = 0;
}

std::string EvtHQET3::getName()
{

  return "HQET3";

}



EvtDecayBase* EvtHQET3::clone()
{

  return new EvtHQET3;

}


void EvtHQET3::decay(EvtParticle* p)
{

  p->initializePhaseSpace(getNDaug(), getDaugs());
  m_calcamp->CalcAmp(p, _amp2, m_hqetffmodel);

}

void EvtHQET3::initProbMax()
{

  EvtId parnum, mesnum, lnum, nunum;

  parnum = getParentId();
  mesnum = getDaug(0);
  lnum = getDaug(1);
  nunum = getDaug(2);

  double mymaxprob = m_calcamp->CalcMaxProb(parnum, mesnum,
                                          lnum, nunum, m_hqetffmodel);

  // Leptons
  static EvtId EM = EvtPDL::getId("e-");
  static EvtId EP = EvtPDL::getId("e+");
  static EvtId MUM = EvtPDL::getId("mu-");
  static EvtId MUP = EvtPDL::getId("mu+");
  static EvtId TAUM = EvtPDL::getId("tau-");
  static EvtId TAUP = EvtPDL::getId("tau+");

  if (lnum == EP || lnum == EM || lnum == MUP || lnum == MUM) {
    setProbMax(mymaxprob);
    return;
  }
  if (lnum == TAUP || lnum == TAUM) {
    setProbMax(6500);
    return;
  }



}


void EvtHQET3::init()
{

  checkNDaug(3);

  //We expect the parent to be a scalar
  //and the daughters to be X lepton neutrino
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(1, EvtSpinType::DIRAC);
  checkSpinDaughter(2, EvtSpinType::NEUTRINO);

  EvtSpinType::spintype d1type = EvtPDL::getSpinType(getDaug(0));
  if (d1type == EvtSpinType::SCALAR) {
    if (getNArg() == 3) {
      m_hqetffmodel = new EvtHQET3FF(getArg(0), getArg(1), getArg(2));
      m_calcamp = new EvtSemiLeptonicScalarAmp;
    } else {
      EvtGenReport(EVTGEN_ERROR, "EvtGen") << "HQET3 model for scalar meson daughters needs 2 arguments. Sorry." << endl;
      ::abort();
    }
  } else if (d1type == EvtSpinType::VECTOR) {
    if (getNArg() == 5) {
      m_hqetffmodel = new EvtHQET3FF(getArg(0), getArg(1), getArg(2), getArg(3), getArg(4));
      m_calcamp = new EvtSemiLeptonicVectorAmp;
    } else  {
      EvtGenReport(EVTGEN_ERROR, "EvtGen") << "HQET3 model for vector meson daughtersneeds 4 arguments. Sorry." << endl;
      ::abort();
    }
  } else {
    EvtGenReport(EVTGEN_ERROR, "EvtGen") << "HQET3 model handles only scalar and vector meson daughters. Sorry." << endl;
    ::abort();
  }


}

