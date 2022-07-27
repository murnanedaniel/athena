/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// Dear emacs, this is -*-c++-*-

#ifndef ELASTICTRUTHTRAJECTORYBUILDER_H
#define ELASTICTRUTHTRAJECTORYBUILDER_H

#include "TrkToolInterfaces/ITruthTrajectoryBuilder.h"
#include "AthenaBaseComps/AthAlgTool.h"

#include "AtlasHepMC/GenVertex_fwd.h"

namespace Trk {
  
  class ElasticTruthTrajectoryBuilder: virtual public ITruthTrajectoryBuilder,
				       public AthAlgTool 
  {
  public:
    ElasticTruthTrajectoryBuilder(const std::string& type, const std::string& name, const IInterface* parent);

    virtual StatusCode initialize();
    
    void buildTruthTrajectory(TruthTrajectory *result, HepMC::ConstGenParticlePtr input) const;
    
    /** Returns an umambiguous mother of the truth particle on a TruthTrajectory, or 0.
     * Implemented via truthTrajectoryCuts().
     */
    HepMC::ConstGenParticlePtr getMother(HepMC::ConstGenParticlePtr particle) const;

    /** Returns an umambiguous daughter of the truth particle on a TruthTrajectory, or 0.
     * This assumes no one mother can correspond to no more than one daughter.
     * Implemented via truthTrajectoryCuts().
     */
    HepMC::ConstGenParticlePtr getDaughter(HepMC::ConstGenParticlePtr particle) const;

    /** Return type for the next method */
    typedef std::pair<HepMC::ConstGenParticlePtr, HepMC::ConstGenParticlePtr> MotherDaughter;

    /** 
     * Decides if the vertex connects two particles on the same
     * TruthTrajectory.  Returns the (mother,daughter) pair,
     * where either of the two can be NULL.
     *
     * This interface implies that there is one-to-one correspondence
     * between mother and daughter particles; i.e. TruthTrajectories
     * are not allowed to branch.
     */
    static MotherDaughter truthTrajectoryCuts(HepMC::ConstGenVertexPtr vtx) ;
  };
  
} // end namespace Trk

#endif/*ELASTICTRUTHTRAJECTORYBUILDER_H*/
