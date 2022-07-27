/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGCONFDATA_L1CONNECTOR_H
#define TRIGCONFDATA_L1CONNECTOR_H

#include "TrigConfData/ConstIter.h"
#include "TrigConfData/DataStructure.h"

#include <map>
#include <vector>

namespace TrigConf {

   /** @brief a TriggerLine entry describes the location of a threshold multiplicity on a cable (connector)
    *
    * for electrical connections from L1Topo boards it also knows
    * which fpga they come from and which clock signal they have (those signals run on doubled clock)
    */
   class TriggerLine {
   public:
   TriggerLine(const std::string & name, unsigned int startbit, unsigned int nbits, unsigned int flatindex, unsigned int fpga=0, unsigned int clock=0, const std::string & connName="") :
     m_name(name), m_startbit(startbit), m_nbits(nbits), m_flatindex(flatindex), m_fpga(fpga), m_clock(clock), m_connName(connName)
      {}
      const std::string & name() const { return m_name; }
      unsigned int startbit() const { return  m_startbit; }
      unsigned int flatindex() const { return  m_flatindex; }
      unsigned int endbit() const { return  m_startbit + m_nbits - 1; }
      unsigned int nbits() const { return m_nbits; }
      unsigned int fpga() const { return m_fpga; }
      unsigned int clock() const { return m_clock; }
      const std::string & connName() const { return m_connName; }
   private:
      std::string m_name;      // the name of the threshold whose multiplicity is transmitted
      unsigned int m_startbit;  // the location on the cable - first bit
      unsigned int m_nbits;     // the location on the cable - number of bits used to encode the multiplicity
      unsigned int m_flatindex; // position of output bit in topo board for a given fpga/clock - first bit
      unsigned int m_fpga;      // for electrical signals from L1Topo boards only: the fpga the signal is coming from
      unsigned int m_clock;     // for electrical signals from L1Topo boards only: the clock of the signal
      std::string m_connName;  // the name of the connector where the triggerline is allocated
   };

   /** @brief L1 connectors configuration */
   class L1Connector final : public DataStructure {
   public:

      enum class ConnectorType { ELECTRICAL, OPTICAL, CTPIN };

      /** Constructor */
      L1Connector();

      L1Connector(const L1Connector &) = delete;
      L1Connector& operator=(const L1Connector&) = delete;
      L1Connector(L1Connector&&) = delete;

      /** Constructor initialized with configuration data 
       * @param data The data containing the L1 menu 
       */
      L1Connector(const std::string & connName, const ptree & data);

      /** Destructor */
      virtual ~L1Connector() override = default;

      virtual std::string className() const override;

      /** Accessor to the number of trigger lines */
      std::size_t size() const;

      std::string type() const;

      /** Accessor to the connector type */
      ConnectorType connectorType() const;

      /** names of all trigger lines */
      std::vector<std::string> triggerLineNames() const;

      /** Accessor to the triggerlines on the connector
       * 
       * For electrical connectors from the L1Topo boards a triggerline vector holds up to 16 signals, which come from 
       * the same fpga and are transmitted at the same clock flank. So in this case the fpga and clock have to be specified.
       * For all other connectors the default value 0 has to be used for fpga and clock
       *
       * @param fpga - the L1Topo fpga (0 or 1)
       * @param clock - the clock of the signal group (0 or 1)
       */
      const std::vector<TrigConf::TriggerLine> & triggerLines(unsigned int fpga = 0, unsigned int clock = 0) const;

      bool hasLine( const std::string & lineName ) const;

      const TrigConf::TriggerLine & triggerLine( const std::string & lineName ) const;

      bool legacy() const { return m_isLegacy; }
      
      [[deprecated("Use legacy() instead.")]]
      bool isLegacy() const { return m_isLegacy; }
      
      std::size_t maxFpga() const { return m_maxFpga; }

      std::size_t maxClock() const { return m_maxClock; }

   private:

      /** Update the internal members */
      virtual void update() override;

      ConnectorType m_type;
      std::vector<TrigConf::TriggerLine> m_triggerLines[2][2];
      std::map<std::string, TrigConf::TriggerLine*> m_lineByName;
      std::size_t m_maxFpga{1};
      std::size_t m_maxClock{1};

      bool m_isLegacy;
   };

}

#endif
