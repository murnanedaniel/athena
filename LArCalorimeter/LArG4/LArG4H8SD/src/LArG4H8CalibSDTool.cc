/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "LArG4H8CalibSDTool.h"

#include "LArG4Code/VCalibrationCalculator.h"

#include "LArG4Code/LArG4CalibSD.h"
#include "LArG4Code/EscapedEnergyProcessing.h"
#include "LArG4Code/CalibrationDefaultCalculator.h"
#include "CaloG4Sim/VEscapedEnergyProcessing.h"
#include "CaloG4Sim/EscapedEnergyRegistry.h"

#include "LArG4Barrel/PresamplerCalibrationCalculator.h"
#include "LArG4Barrel/LArBarrelCalibrationCalculator.h"

#include "LArTBCryostatCalibrationCalculator.h"
#include "H8CalibrationDefaultCalculator.h"

#include "CaloG4Sim/CalibrationDefaultProcessing.h"
#include "CxxUtils/make_unique.h"


LArG4H8CalibSDTool::LArG4H8CalibSDTool(const std::string& type, const std::string& name, const IInterface *parent)
  : LArG4SDTool(type,name,parent)
  , m_HitColl("LArH8CalibrationHits")
  , m_barCrySD (nullptr)
  , m_bpInSD (nullptr)
  , m_bpDeadSD (nullptr)
  //, m_bpModSD (nullptr)
  , m_barInSD (nullptr)
  , m_barDeadSD (nullptr)
  , m_uninstSD (nullptr)
{
    declareProperty("BarrelCryDeadVolumes",m_barCryVolumes);
    declareProperty("BarrelPreInactiveVolumes",m_bpInVolumes);
    declareProperty("BarrelPreDeadVolumes",m_bpDeadVolumes);
    declareProperty("BarrelPreCalibVolumes",m_bpCalibVolumes);
    declareProperty("BarrelInVolumes",m_barInVolumes);
    declareProperty("BarrelDeadVolumes",m_barDeadVolumes);
}

StatusCode LArG4H8CalibSDTool::initializeSD()
{
    m_barCrySD  = new LArG4CalibSD( "LAr::TBBarrelCryostat::Dead" , new LArTBCryostatCalibrationCalculator() );
    m_bpInSD    = new LArG4CalibSD( "LAr::TBBarrel::Presampler::Inactive" , new LArG4::BarrelPresampler::CalibrationCalculator() );
    m_bpDeadSD  = new LArG4CalibSD( "LAr::TBBarrel::Presampler::Dead" , new LArG4::BarrelPresampler::CalibrationCalculator() );
//    m_bpModSD   = new LArG4CalibSD( "LAr::TBBarrel::Presampler::Module::Calibration" , new LArG4::BarrelPresampler::CalibrationCalculator() );
    m_barInSD   = new LArG4CalibSD( "LAr::TBBarrel::Inactive" , new LArG4::Barrel::CalibrationCalculator() );
    m_barDeadSD = new LArG4CalibSD( "LAr::TBBarrel::Dead" , new LArG4::Barrel::CalibrationCalculator() );

    m_uninstSD  = new LArG4CalibSD( "Default::Dead::Uninstrumented::Calibration::Region" , new H8CalibrationDefaultCalculator() );
    
    std::map<G4VSensitiveDetector*,std::vector<std::string>*> configuration;
    configuration[m_barCrySD]  = &m_barCryVolumes;
    configuration[m_bpInSD]    = &m_bpInVolumes;
    configuration[m_bpDeadSD]  = &m_bpDeadVolumes;
//    configuration[m_bpModSD]   = &m_bpModVolumes;
    configuration[m_barInSD]   = &m_barInVolumes;
    configuration[m_barDeadSD] = &m_barDeadVolumes;
    configuration[m_uninstSD]    = new std::vector<std::string>; // No volumes for this guy
    setupAllSDs(configuration);

    setupHelpers(m_barCrySD);
    setupHelpers(m_bpInSD);
    setupHelpers(m_bpDeadSD);
//    setupHelpers(m_bpModSD);
    setupHelpers(m_barInSD);
    setupHelpers(m_barDeadSD);
    setupHelpers(m_uninstSD);

    ATH_MSG_DEBUG( " **** Call  LArG4H8CalibSD::initializeProcessing() " );
 
    // ADS not needed with new user actions.
//    ATH_MSG_INFO("CalibrationDefaultProcessing");
//    CaloG4::CalibrationDefaultProcessing * cdp = new CaloG4::CalibrationDefaultProcessing(); // Registers itself
//    ATH_MSG_INFO("setting default SD");
//    cdp->SetDefaultSD(m_uninstSD);
 
    // Initialize the escaped energy processing for LAr volumes.
    CaloG4::VEscapedEnergyProcessing* eep = new EscapedEnergyProcessing(m_uninstSD);
    CaloG4::EscapedEnergyRegistry* registry = CaloG4::EscapedEnergyRegistry::GetInstance();
    registry->AddAndAdoptProcessing( "LAr::", eep );
    
    return StatusCode::SUCCESS;
}


StatusCode LArG4H8CalibSDTool::Gather()
{
    // In this case, *unlike* other SDs, the *tool* owns the collection
#ifdef ATHENAHIVE
    // Temporary fix for Hive until isValid is fixed
    m_HitColl = CxxUtils::make_unique<CaloCalibrationHitContainer>(m_HitColl.name());
#else
    if (!m_HitColl.isValid()) m_HitColl = CxxUtils::make_unique<CaloCalibrationHitContainer>(m_HitColl.name());
#endif
    m_barCrySD  ->EndOfAthenaEvent( &*m_HitColl );
    m_bpInSD    ->EndOfAthenaEvent( &*m_HitColl );
    m_bpDeadSD  ->EndOfAthenaEvent( &*m_HitColl );
//    m_bpModSD   ->EndOfAthenaEvent( &*m_HitColl );
    m_barInSD   ->EndOfAthenaEvent( &*m_HitColl );
    m_barDeadSD ->EndOfAthenaEvent( &*m_HitColl );
    m_uninstSD  ->EndOfAthenaEvent( &*m_HitColl );
    return StatusCode::SUCCESS;
}
