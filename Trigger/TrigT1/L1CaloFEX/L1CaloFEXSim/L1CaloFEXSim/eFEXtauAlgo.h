/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

//***************************************************************************
//                           eFEXtauAlgo.h  -  
//                              -------------------
//     begin                : 06 05 2020
//     email                : nicholas.andrew.luongo@cern.ch
//  ***************************************************************************/


#ifndef eFEXtauAlgo_H
#define eFEXtauAlgo_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "L1CaloFEXToolInterfaces/IeFEXtauAlgo.h"
#include "AthenaKernel/CLASS_DEF.h"
#include "L1CaloFEXSim/eFEXtauTOB.h"
#include "L1CaloFEXSim/eTowerContainer.h"

namespace LVL1 {
  
  //Doxygen class description below:
  /** The eFEXtauAlgo class calculates the tau TOB variables
  */
  
  class eFEXtauAlgo : public AthAlgTool, virtual public IeFEXtauAlgo {
    
  public:
    /** Constructors */
    eFEXtauAlgo(const std::string& type, const std::string& name, const IInterface* parent);

    /** standard Athena-Algorithm method */
    virtual StatusCode initialize() override;
    
    /** Destructor */
    virtual ~eFEXtauAlgo();

    virtual StatusCode safetyTest() override;
    virtual void setup(int inputTable[3][3], int efex_id, int fpga_id, int central_eta) override;

    /** standard Athena-Algorithm method */
    //virtual StatusCode initialize();

    /** standard Athena-Algorithm method */
    //virtual StatusCode finalize();

    virtual bool isCentralTowerSeed() override;
    virtual std::unique_ptr<eFEXtauTOB> getTauTOB() override;
    virtual unsigned int rCoreCore() override;
    virtual unsigned int rCoreEnv() override;
    virtual void getRCore(std::vector<unsigned int> & rCoreVec) override;
    virtual float getRealRCore() override;
    virtual unsigned int rHadCore() override;
    virtual unsigned int rHadEnv() override;
    virtual void getRHad(std::vector<unsigned int> & rHadVec) override;
    virtual float getRealRHad() override;
    virtual unsigned int getEt() override;
    virtual unsigned int getBitwiseEt() override;

  protected:

  private:
    int m_eFexalgoTowerID[3][3];

    void buildLayers(int efex_id, int fpga_id, int central_eta);
    void setSupercellSeed();
    void setUnDAndOffPhi();
    virtual bool getUnD() override;
    virtual unsigned int getSeed() override;
	
    unsigned int m_em0cells[3][3];
    unsigned int m_em1cells[12][3];
    unsigned int m_em2cells[12][3];
    unsigned int m_em3cells[3][3];
    unsigned int m_hadcells[3][3];
    unsigned int m_twrcells[3][3];
    unsigned int m_seed;
    bool m_cellsSet = false;
    bool m_und;
    unsigned int m_offPhi;

    SG::ReadHandleKey<LVL1::eTowerContainer> m_eTowerContainerKey {this, "MyETowers", "eTowerContainer", "Input container for eTowers"};

  };
  
} // end of namespace

//CLASS_DEF( LVL1::eFEXtauAlgo , 140708609 , 1 )

#endif
