/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SHOWER_LIB_FCALDIST_ETA_ENERGY_SHOWER_LIB_H
#define SHOWER_LIB_FCALDIST_ETA_ENERGY_SHOWER_LIB_H
// STL includes
#include <map>
#include <vector>

// local includes
#include "LArG4ShowerLib/IShowerLib.h"

// G4 forward declarations
class G4Track;

// Namespace for the ShowerLib related classes
namespace ShowerLib {

  // forward declarations in namespace

  /**
   *
   *   @short Class for shower library shower lib
   *
   *          Store complete shower library (collection of energy bins)
   *
   *  @author Wolfgang Ehrenfeld, University of Hamburg, Germany
   *  @author Sasha Glazov, DESY Hamburg, Germany
   *
   * @version \$Id: FCALDistEtaEnergyShowerLib.h 769594 2016-08-23 13:48:34Z ssnyder $
   *
   */

  class FCALDistEtaEnergyShowerLib : public IShowerLib {

  public:

    //! factory method. create a library from root file. returns NULL if file is invalid.
    static IShowerLib* readFromROOTFile(TFile* source);

    //! factory method. create empty library with the given structure. returns NULL if file is invalid.
    static IShowerLib* createEmptyLib(const std::string& inputFile);

    //! default destructor
    virtual ~FCALDistEtaEnergyShowerLib() {}

    //! get shower for given G4 track
    virtual std::vector<EnergySpot>* getShower(const G4Track* track, ShowerLibStatistics* stats, int randomShift) const;
    //! get average length of showers for the given energy
    virtual double getContainmentZ(const G4Track* track) const;
    //! get average lateral spread of the showers for the given energy
    virtual double getContainmentR(const G4Track* track) const;
    //! store shower in the library
    virtual bool storeShower(const HepMC::GenParticle* genParticle,const Shower* shower);
    //! write library to ROOT file
    virtual bool writeToROOT(TFile* dest);

    virtual ShowerLibStatistics* createStatistics() const;

    inline virtual const std::string getName() const { return "FCALDist Eta Energy ShowerLib"; }

  protected:

    virtual const std::string printParameters() const;

  private:

    FCALDistEtaEnergyShowerLib():m_xrodcent(0.0)
                                ,m_yrodcent(0.0)
                                ,m_step(0.0)
#ifndef __FSLIB_NO_BACKWARD_COMPAT__
                                ,m_compat(false)
#endif
    {}

    //! read library from given TTree
    bool read(TTree* source);
    //! write library to given TTree
    bool write(TTree* dest) const;
    //!
    bool readStructure(std::map<float, std::vector<float> >& structure);

    float distance(double x, double y) const;

    typedef std::map<float,Shower> distbin;
    typedef std::map<float,distbin> etabin;
    typedef std::map<float,etabin> library;

    library m_libData;

    //distance calculator parameters
    double m_xrodcent;
    double m_yrodcent;
    double m_step;

#ifndef __FSLIB_NO_BACKWARD_COMPAT__
    bool m_compat;
#endif

  };

} // namespace ShowerLib

#endif // SHOWER_LIB_FCALDIST_SHOWER_LIB_H
