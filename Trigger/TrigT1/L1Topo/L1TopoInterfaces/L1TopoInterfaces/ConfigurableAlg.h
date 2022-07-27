/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/
//  BaseTOBoAlg.h
//  TopoCore
//  Created by Joerg Stelzer on 11/10/12.

#ifndef TCS_ConfigurableAlg
#define TCS_ConfigurableAlg

#include "TrigConfBase/TrigConfMessaging.h"

#include "L1TopoInterfaces/Parameter.h"
#include "L1TopoInterfaces/ParameterSpace.h"
#include "L1TopoInterfaces/AlgFactory.h"
#include "L1TopoCommon/StatusCode.h"
#include "L1TopoEvent/GenericTOB.h"

#include <vector>
#include <string>
#include <memory>

class TH1;
class TH2;
class IL1TopoHistSvc;

namespace TCS {

   class ConfigurableAlg : public TrigConf::TrigConfMessaging {
   protected:
      enum AlgType { NONE, SORT, DECISION, COUNT };

   public:

      // constructor
      ConfigurableAlg(const std::string & name, AlgType algType);

      // destructor
      virtual ~ConfigurableAlg();

      // setting the class name
      void setClassName(const std::string & className) { m_className = className; }

      void setLegacyMode(bool isLegacyTopo) {m_isLegacyTopo=isLegacyTopo;}

      // accessors
      const std::string & name() const { return m_name; }
      const std::string & className() const { return m_className; }
      std::string fullname() const { return m_className + "/" + m_name; }

      // the index of the algorithm in the hardware
      unsigned int algoId() const { return m_algoId; }

      // reset
      virtual StatusCode initialize() { return StatusCode::SUCCESS;}

      StatusCode reset() { return doReset(); }

      bool isSortingAlg() const { return m_algType == SORT; }

      bool isDecisionAlg() const { return m_algType == DECISION; }

      bool isCountingAlg() const { return m_algType == COUNT; }

      bool isLegacyTopo() const { return m_isLegacyTopo; }

      // Kinematic calculation
      unsigned int calcDeltaPhiBW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcDeltaEtaBW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcInvMassBW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcTMassBW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcDeltaR2BW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned long quadraticSumBW(int i1, int i2);
      unsigned int calcDeltaPhi(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcDeltaEta(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcInvMass(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcTMass(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);
      unsigned int calcDeltaR2(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2);


      // const access to parameter
      const Parameter & parameter(const std::string & parameterName) const;

      const Parameter & parameter(const std::string & parameterName, unsigned int selection) const;

      const ParameterSpace & parameters() const { return m_parameters; }

      // set parameters from configuration
      void setParameters(const ParameterSpace &);
      
      // setters
      void setAlgoId(unsigned int algoId) { m_algoId = algoId; }

      /**
       * sets the external hist service
       */
      void setL1TopoHistSvc(std::shared_ptr<IL1TopoHistSvc>);

      void bookHistMult(std::vector<std::string> &regName, const std::string& name, const std::string& title, const std::string& xtitle, const int binx, const int xmin, const int xmax);
      void bookHist(std::vector<std::string> &regName, const std::string& name, const std::string& title, const int binx, const int xmin, const int xmax);
      void bookHistMult(std::vector<std::string> &regName, const std::string& name, const std::string& title, const std::string& xtitle, const std::string& ytitle, const int binx, const int xmin, const int xmax, const int biny, const int ymin, const int ymax);
      void bookHist(std::vector<std::string> &regName, const std::string& name, const std::string& title, const int binx, const int xmin, const int xmax, const int biny, const int ymin, const int ymax);


   protected:

      // define parameters by developer
      void defineParameter(const std::string & name, TCS::parType_t value);

      void defineParameter(const std::string & name, TCS::parType_t value, unsigned int selection);
      
      virtual StatusCode doReset() = 0;
      
      void registerHist(TH1 *);

      void registerHist(TH2 *);

      void fillHist1D(const std::string & histName, double x);

      void fillHist2D(const std::string & histName, double x, double y);
      
      bool isocut(const std::string& threshold, const unsigned int bit);
      
      bool isocut(const unsigned int threshold, const unsigned int bit);
   private:

      class ConfigurableAlgImpl;
      std::unique_ptr<ConfigurableAlgImpl> m_impl;

      void defineParameter(const Parameter &);
      std::string ToString(const int val);
     
      std::string    m_name {""};
      std::string    m_className {""};

      // just for the hardware
      unsigned int m_algoId {0};

      ParameterSpace m_parameters {""};

      AlgType  m_algType; // stores type of alg (Sorting or Decision)

      bool m_isLegacyTopo;

   };

std::ostream & operator<<(std::ostream &, const TCS::ConfigurableAlg &);

} // end of namespace TCS

#endif /* defined(__TopoCore__TopoAlg__) */
