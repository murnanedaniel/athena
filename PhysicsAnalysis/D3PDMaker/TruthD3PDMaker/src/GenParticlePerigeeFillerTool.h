/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRACKD3PDMAKER_GENPARTICLEPERIGEEFILLERTOOL_H
#define TRACKD3PDMAKER_GENPARTICLEPERIGEEFILLERTOOL_H

#include "D3PDMakerUtils/BlockFillerTool.h"
#include "AtlasHepMC/GenParticle.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkToolInterfaces/ITruthToTrack.h"

namespace D3PD {
#ifdef HEPMC3
class GenParticlePerigeeFillerTool: public D3PD::BlockFillerTool<HepMC::ConstGenParticlePtr>
{
public:
  typedef D3PD::BlockFillerTool<HepMC::ConstGenParticlePtr> Base;
#else
class GenParticlePerigeeFillerTool : public D3PD::BlockFillerTool<HepMC::GenParticle>
{
public:
  typedef D3PD::BlockFillerTool<HepMC::GenParticle> Base;
#endif

  GenParticlePerigeeFillerTool (const std::string& type,
			const std::string& name,
			const IInterface* parent);

  virtual StatusCode initialize() override;
  virtual StatusCode book() final;

#ifdef HEPMC3
  virtual StatusCode fill (const HepMC::ConstGenParticlePtr& p) override;
#else
  virtual StatusCode fill (const HepMC::GenParticle& p) override;
#endif

private:

  void clearData();

  /* Tool to generate perigee parameters for a generated particle **/
  ToolHandle<Trk::ITruthToTrack> m_truthToTrack;

  /* NTuple variables: **/
  // truth track parameters at perigee
  int   *m_ok;
  float *m_d0;
  float *m_z0;
  float *m_phi;
  float *m_theta;
  float *m_qoverp;

}; // class GenParticlePerigeeFillerTool

} // namespace D3PD

#endif
