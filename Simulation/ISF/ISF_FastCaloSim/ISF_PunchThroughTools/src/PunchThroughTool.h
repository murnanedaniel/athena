/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_PUNCHTHROUGHTOOLS_SRC_PUNCHTHROUGHTOOL_H
#define ISF_PUNCHTHROUGHTOOLS_SRC_PUNCHTHROUGHTOOL_H 1

#include "ISF_FastCaloSimInterfaces/IPunchThroughTool.h"
#include <string>

#include "ISF_FastCaloSimInterfaces/IPunchThroughClassifier.h"

// Athena Base
#include "AthenaBaseComps/AthAlgTool.h"

//Barcode
#include "BarcodeInterfaces/IBarcodeSvc.h"

//Geometry
#include "SubDetectorEnvelopes/IEnvelopeDefSvc.h"

#include "ISF_Interfaces/IGeoIDSvc.h"

// Gaudi & StoreGate
#include "GaudiKernel/IPartPropSvc.h"

#include "BarcodeEvent/Barcode.h"
#include "BarcodeEvent/PhysicsProcessCode.h"
#include "GeoPrimitives/GeoPrimitives.h"

#include "ISF_Event/ISFParticleContainer.h"

#include "AtlasHepMC/GenEvent_fwd.h"

//libXML
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

/*-------------------------------------------------------------------------
 *  Forward declarations
 *-------------------------------------------------------------------------*/

class TFile;

namespace HepPDT {
  class ParticleDataTable;
}

namespace ISF {
  class PunchThroughParticle;
  class PDFcreator;
}

namespace ISF {

  class PunchThroughTool : public extends<AthAlgTool, IPunchThroughTool>
  {
  public:
    /** Constructor */
    PunchThroughTool(const std::string&,const std::string&,const IInterface*);

    /** Destructor */
    virtual ~PunchThroughTool () = default;

    /** AlgTool initialize method */
    virtual StatusCode initialize();
    /** AlgTool finalize method */
    virtual StatusCode finalize  ();
    /** interface function: fill a vector with the punch-through particles */
    const ISF::ISFParticleVector* computePunchThroughParticles(const ISF::ISFParticle &isfp, const TFCSSimulationState& simulstate, CLHEP::HepRandomEngine* rndmEngine) const;

  private:
    /*---------------------------------------------------------------------
     *  Private member functions
     *---------------------------------------------------------------------*/
    /** registers a type of punch-through particles which will be simulated */
    StatusCode registerParticle(int pdgID, bool doAntiparticle = false,
                                double minEnergy = 0., int maxNumParticles = -1,
                                double numParticlesFactor = 1., double energyFactor = 1.,
                                double posAngleFactor = 1.,
                                double momAngleFactor = 1.);
    /** register a correlation for the two given types of punch-through particles
        with a given energy threshold above which we will have full correlation */
    StatusCode registerCorrelation(int pdgID1, int pdgID2,double minCorrEnergy = 0., double fullCorrEnergy = 0.);

    /** reads out the lookuptable for the given type of particle */
    std::unique_ptr<ISF::PDFcreator> readLookuptablePDF(int pdgID, const std::string& folderName);

    /** create the right number of punch-through particles for the given pdg
     *  and return the number of particles which was created. also create these
     *  particles with the right distributions (energy, theta, phi).
     *  if a second argument is given, create exactly this number of particles
     *  (also with the right energy,theta,phi distributions */
    int getAllParticles(const ISF::ISFParticle &isfp, ISFParticleVector& isfpCont, CLHEP::HepRandomEngine* rndmEngine, int pdg, double interpEnergy, double interpEta, int numParticles = -1) const;

    /** get the right number of particles for the given pdg while considering
     *  the correlation to an other particle type, which has already created
     *  'corrParticles' number of particles */
    int getCorrelatedParticles(const ISF::ISFParticle &isfp, ISFParticleVector& isfpCont, int doPdg, int corrParticles, CLHEP::HepRandomEngine* rndmEngine, double interpEnergy, double interpEta) const;

    /** create exactly one punch-through particle with the given pdg and the given max energy */
    ISF::ISFParticle *getOneParticle(const ISF::ISFParticle &isfp, int pdg, CLHEP::HepRandomEngine* rndmEngine, double interpEnergy, double interpEta) const;

    /** create a ISF Particle state at the MS entrace containing a particle with the given properties */
    ISF::ISFParticle *createExitPs(const ISF::ISFParticle &isfp, int PDGcode, double energy, double theta, double phi, double momTheta, double momPhi) const;

    /** get the floating point number in a string, after the given pattern */
    double getFloatAfterPatternInStr(const char *str, const char *pattern);
    /** get particle through the calorimeter */
    Amg::Vector3D propagator(double theta, double phi) const;

    //apply the inverse PCA transform
    std::vector<double> inversePCA(std::vector<double> &variables) const;

    //apply the inverse CDF trainsform
    double inverseCdfTransform(double variable, std::map<double, double> inverse_cdf_map) const;

    //dot product between matrix and vector, used to inverse PCA
    std::vector<double> dotProduct(const std::vector<std::vector<double>> &m, const std::vector<double> &v) const;

    //returns normal cdf based on normal gaussian value
    double normal_cdf(double x) const;

    //apply energy interpolation
    double interpolateEnergy(const double &energy, CLHEP::HepRandomEngine* rndmEngine) const;

    //apply eta interpolation
    double interpolateEta(const double &eta, CLHEP::HepRandomEngine* rndmEngine) const;

    //load inverse quantile transformer from XML
    StatusCode initializeInverseCDF(std::string quantileTransformerConfigFile);

    //get CDF mapping for individual XML node
    std::map<double, double> getVariableCDFmappings(xmlNodePtr& nodeParent);

    //load inverse PCA from XML
    StatusCode initializeInversePCA(std::string inversePCAConfigFile);


    /*---------------------------------------------------------------------
     *  Private members
     *---------------------------------------------------------------------*/

     /** energy and eta points in param */
     std::vector<double>                 m_energyPoints;
     std::vector<double>                 m_etaPoints;

    /** calo-MS borders */
    double                               m_R1{0.};
    double                               m_R2{0.};
    double                               m_z1{0.};
    double                               m_z2{0.};

    /** ParticleDataTable needed to get connection pdg_code <-> charge */
    const HepPDT::ParticleDataTable*    m_particleDataTable{nullptr};

    /** ROOT objects */
    TFile*                              m_fileLookupTable{nullptr};   //!< the punch-through lookup table file

    /** needed to create punch-through particles with the right distributions */
    std::map<int, PunchThroughParticle*> m_particles;       //!< store all punch-through information for each particle id

    /*---------------------------------------------------------------------
     *  Properties
     *---------------------------------------------------------------------*/

    /** Properties */
    std::string                          m_filenameLookupTable{"CaloPunchThroughParametrisation.root"};     //!< holds the filename of the lookup table (property)
    std::string                          m_filenameInverseCDF;     //!< holds the filename of inverse quantile transformer config
    std::string                          m_filenameInversePCA;     //!< holds the filename of inverse PCA config

    ToolHandle< IPunchThroughClassifier >     m_punchThroughClassifier;
    std::vector<int>                     m_pdgInitiators;           //!< vector of punch-through initiator pgds
    std::vector<int>                     m_initiatorsMinEnergy;     //!< vector of punch-through initiator min energyies to create punch through
    std::vector<double>                  m_initiatorsEtaRange;      //!< vector of min and max abs eta range to allow punch through initiators
    std::vector<int>                     m_punchThroughParticles;   //!< vector of pdgs of the particles produced in punch-throughs
    std::vector<bool>                    m_doAntiParticles;         //!< vector of bools to determine if anti-particles are created for each punch-through particle type
    std::vector<int>                     m_correlatedParticle;      //!< holds the pdg of the correlated particle for each given pdg
    std::vector<double>                  m_minCorrEnergy;           //!< holds the energy threshold below which no particle correlation is computed
    std::vector<double>                  m_fullCorrEnergy;          //!< holds the energy threshold above which a particle correlation is fully developed
    std::vector<double>                  m_posAngleFactor;          //!< tuning parameter to scale the position deflection angles
    std::vector<double>                  m_momAngleFactor;          //!< tuning parameter to scale the momentum deflection angles
    std::vector<double>                  m_minEnergy;               //!< punch-through particles minimum energies
    std::vector<int>                     m_maxNumParticles;         //!< maximum number of punch-through particles for each particle type
    std::vector<double>                  m_numParticlesFactor;      //!< scale the number of punch-through particles
    std::vector<double>                  m_energyFactor;            //!< scale the energy of the punch-through particles

    /*---------------------------------------------------------------------
     *  Constants
     *---------------------------------------------------------------------*/
     constexpr static double m_SQRT_0p5 = std::sqrt(0.5);
     constexpr static double m_SQRT_2 = std::sqrt(2);

    /*---------------------------------------------------------------------
     *  ServiceHandles
     *---------------------------------------------------------------------*/
    ServiceHandle<IPartPropSvc>          m_particlePropSvc;         //!< particle properties svc
    ServiceHandle<IGeoIDSvc>             m_geoIDSvc;
    ServiceHandle<Barcode::IBarcodeSvc>  m_barcodeSvc;
    ServiceHandle<IEnvelopeDefSvc>       m_envDefSvc;

    /** beam pipe radius */
    double                              m_beamPipe{500.};

    std::vector<std::vector<double>> m_inverse_PCA_matrix;
    std::vector<double> m_PCA_means;

    std::map<double, double>  m_variable0_inverse_cdf;
    std::map<double, double>  m_variable1_inverse_cdf;
    std::map<double, double>  m_variable2_inverse_cdf;
    std::map<double, double>  m_variable3_inverse_cdf;
    std::map<double, double>  m_variable4_inverse_cdf;
  };
}

#endif
