/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#ifndef TRIGT1CTP_CTPSIMULATION_H
#define TRIGT1CTP_CTPSIMULATION_H

// Basic includes:
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/Property.h"
#include "StoreGate/DataHandle.h"

// monitoring from HLT
#include "TrigInterfaces/IMonitoredAlgo.h"
#include "TrigConfigSvc/ITrigConfigSvc.h"

// local includes:
#include "TrigT1CTP/ISpecialTrigger.h"


// Forward includes:
class StoreGateSvc;
class IAtRndmGenSvc;
class IMonitorToolBase;

namespace TrigConf {
  class PIT;
  class ILVL1ConfigSvc;
  class TrigConfCoolFolderSpec;  
  class BunchGroup;
}

namespace LVL1 {
  class MuCTPICTP;
  class EmTauCTP;
  class JetCTP;
  class EnergyCTP;
  class MbtsCTP;
  class BcmCTP;
  class LucidCTP;
  class ZdcCTP;
  class BptxCTP;
  class NimCTP;
}

namespace LVL1MUCTPI{

}

/// Namespace for the CTP related classes
namespace LVL1CTP {

  // Forward include(s) in local namespace:
  class ThresholdMap;
  class ItemMap;
  class ResultBuilder;

  /**
   *
   *   @short CTP simulation algorithm
   *
   * This is the simulation of the final CTP hardware's functions. It is not an actual
   * hardware simulation (unlike that of the MuCTPI), "only" the functioning of the CTP is
   * simulated. The algorithm was developed to be used with the interface to the trigger
   * configuration database. 
   *
   *  @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
   *  @author Wolfgang Ehrenfeld <Wolfgang.Ehrenfeld@desy.de>
   *
   * @version \$Id: CTPSimulation.h,v 1.20 2009-01-29 21:14:24 efeld Exp $
   *
   */

  class CTPSimulation : public AthAlgorithm, public IMonitoredAlgo {

  public:

    CTPSimulation( const std::string& name, ISvcLocator* pSvcLocator );
    ~CTPSimulation();

    virtual StatusCode initialize();
    virtual StatusCode start();
    virtual StatusCode beginRun();
    virtual StatusCode execute();
    virtual StatusCode finalize();
    virtual StatusCode callback(IOVSVC_CALLBACK_ARGS_P(idx, keys));
    //! handle to result builder (for easy data access) 
    const ResultBuilder* resultBuilder() const;
    StoreGateSvc *m_detStore;

  private:

    //! Function pointer to the correct multiplicity extraction function:
    StatusCode ( LVL1CTP::CTPSimulation::*m_extractFunction ) ( void );

    //! Function to extract the multiplicity information sent by the sub-systems
    StatusCode extractMultiplicities();
    StatusCode LoadBunchGroups();   

    // Needed services:
    ServiceHandle<TrigConf::ILVL1ConfigSvc> m_configSvc;    //!< property, see @link CTPSimulation::CTPSimulation @endlink
    ServiceHandle<IAtRndmGenSvc> m_rndmSvc;                 //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_rndmEngineName;                           //!< property, see @link CTPSimulation::CTPSimulation @endlink 

    // Needed tools:
    ToolHandleArray < IMonitorToolBase > m_monitors;        //!< property, see @link CTPSimulation::CTPSimulation @endlink 

    // Inputs to the CTP:
    const DataHandle< LVL1::MuCTPICTP > m_muctpiCTP;        //!< MUCTPI input
    const DataHandle< LVL1::EmTauCTP >  m_emtauCTP;         //!< EmTau input
    const DataHandle< LVL1::JetCTP >    m_jetCTP;           //!< Jet input
    const DataHandle< LVL1::EnergyCTP > m_energyCTP;        //!< Energy input
    const DataHandle< LVL1::MbtsCTP >   m_mbtsACTP;         //!< MBTSA input
    const DataHandle< LVL1::MbtsCTP >   m_mbtsCCTP;         //!< MBTSC input
    const DataHandle< LVL1::BcmCTP >    m_bcmCTP;           //!< BCM input
    const DataHandle< LVL1::LucidCTP >  m_lucidCTP;         //!< LUCID input
    const DataHandle< LVL1::ZdcCTP >    m_zdcCTP;           //!< ZDC input
    const DataHandle< LVL1::NimCTP >    m_nimCTP;           //!< NIM input
    const DataHandle< LVL1::BptxCTP >   m_bptxCTP;           //!< BPTX input

    
    // Maps between the different trigger objects:
    ThresholdMap*       m_decisionMap;                      //!< pointer to threshold map
    ItemMap*            m_itemMap;                          //!< pointer to item map
    InternalTriggerMap  m_internalTrigger;                  //!< pointer to internal trigger map

    // Properties: stearing flags
    bool        m_doCalo;                                   //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doMBTS;                                   //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doMBTSSI;                                 //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doMuon;                                   //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doBCM;                                          //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doLUCID;                                  //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doZDC;                                    //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doBPTX;                                   //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doNIM;                                    //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doRNDM;                                    //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_doPSCL;                                    //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_IsData;                                   //!< property, see @link CTPSimulation::CTPSimulation @endlink
    IntegerProperty m_prescaleMode;                         //!< property, see @link CTPSimulation::CTPSimulation @endlink
    // Properties: StoreGate location of input
    std::string m_jetEnergyConfLoc;                         //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_muonCTPLoc;                               //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_emTauCTPLoc;                              //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_jetCTPLoc;                                //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_energyCTPLoc;                             //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_mbtsACTPLoc;                              //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_mbtsCCTPLoc;                              //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_bcmCTPLoc;                                //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_lucidCTPLoc;                              //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_zdcCTPLoc;                                //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_bptxCTPLoc;                               //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_nimCTPLoc;                                //!< property, see @link CTPSimulation::CTPSimulation @endlink
    // Properties: StoreGate location of output
    std::string m_roiOutputLoc;                             //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_roiOutputLoc_Rerun;                             //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_rdoOutputLoc;                             //!< property, see @link CTPSimulation::CTPSimulation @endlink
    std::string m_rdoOutputLoc_Rerun;                       //!< property, see @link CTPSimulation::CTPSimulation @endlink
    // properties: flags for internal triggers
    bool        m_applyDeadTime;                            //!< property, see @link CTPSimulation::CTPSimulation @endlink
    bool        m_applyBunchGroup;                          //!< property, see @link CTPSimulation::CTPSimulation @endlink
    // Properties: detectorStore location of bunchgroups
    std::string m_BunchGroupLoc;                            //!< property, see @link CTPSimulation::CTPSimulation @endlink
    /**
     * If IntroduceArtificialTriggerOffsets is set to true, and a text
     * file with the configuration of the trigger offsets is defined via
     * OffsetConfigFile, a readout window of m_readoutWindow BCs is
     * introduced to insert artificial trigger offsets in the simulated
     * bytestream data to mimic timing problems that might ocurr in real data.
    */
    uint32_t m_readoutWindow;                               //!< size of readout window to be inserted in BCs
    bool m_introduceArtificialTriggerOffsets;               //!< status flag for introducing artificial trigger offsets
    std::string m_offsetConfigFile;                         //!< path to configuration file for trigger offsets

    ResultBuilder* m_resultBuilder;                         //!< handle to result builder (for easy data access)

    // additional monitoring variables
    std::vector<int> m_prescales;                           //!< prescale set (for monitoring)

  }; // class CTPSimulation

} // namespace LVL1CTP

#endif // TRIGT1CTP_CTPSIMULATION_H
