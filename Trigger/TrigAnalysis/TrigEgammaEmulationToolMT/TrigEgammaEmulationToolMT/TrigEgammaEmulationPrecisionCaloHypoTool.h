/*
 *   Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
 *   */


#ifndef TrigEgammaEmulationPrecisionCaloHypoTool_h
#define TrigEgammaEmulationPrecisionCaloHypoTool_h

#include "AsgTools/AsgTool.h"
#include "TrigEgammaEmulationToolMT/TrigEgammaEmulationToolMT.h"
#include "TrigEgammaEmulationToolMT/TrigEgammaEmulationBaseHypoTool.h"
#include "TrigEgammaEmulationToolMT/ITrigEgammaEmulationBaseHypoTool.h"


namespace Trig{

  class TrigEgammaEmulationPrecisionCaloHypoTool: public TrigEgammaEmulationBaseHypoTool,
                                                virtual public ITrigEgammaEmulationBaseHypoTool
  { 
    ASG_TOOL_CLASS( TrigEgammaEmulationPrecisionCaloHypoTool , ITrigEgammaEmulationBaseHypoTool)

    public:

      TrigEgammaEmulationPrecisionCaloHypoTool(const std::string& myname);
      ~TrigEgammaEmulationPrecisionCaloHypoTool()=default;

      virtual bool emulate( const TrigData &input, bool &pass) const override;
      
    private:

      int findCutIndex( float eta ) const;

      Gaudi::Property< float > m_detacluster { this, "dETACLUSTERthr", 0. , "" };
      Gaudi::Property< float > m_dphicluster { this, "dPHICLUSTERthr", 0. , "" };  
      Gaudi::Property< float > m_et2thr      { this, "ET2Thr"        , 0. , "" };  
      Gaudi::Property< std::vector<float> > m_etabin { this, "EtaBins", {} , "Bins of eta" }; 
      Gaudi::Property< std::vector<float> > m_eTthr { this, "ETthr", {}, "ET Threshold" };

  };

}//namespace
#endif
