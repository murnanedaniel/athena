/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRKDUMPALGS_DUMPOBJECTS_H
#define TRKDUMPALGS_DUMPOBJECTS_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"

#include "InDetIdentifier/PixelID.h"
#include "PixelReadoutGeometry/PixelDetectorManager.h"
#include "InDetIdentifier/SCT_ID.h"
#include "SCT_ReadoutGeometry/SCT_DetectorManager.h"

#include "GaudiKernel/IPartPropSvc.h"

#include <fstream>
#include <iostream>


class PixelID;
class SCT_ID;
class TTree;

namespace InDetDD {
  class PixelDetectorManager;
  class SCT_DetectorManager;
}

namespace HepPDT{
  class ParticleDataTable;
}

namespace HepMC {
  class GenParticle;
}

class DumpObjects : public AthAlgorithm {
  public:
    DumpObjects(const std::string& name, ISvcLocator* pSvcLocator);
    ~DumpObjects(){}
   
    virtual StatusCode initialize() override final;
    virtual StatusCode execute() override final;
    virtual StatusCode finalize() override final;

  private:
    bool isPassed(const HepMC::GenParticle* particle,
                  float& px, float& py, float& pz, float& pt, float& eta, 
                  float& vx, float& vy, float& vz, float& radius, 
                  float& status, float& charge,
		  std::vector<int> &vParentID, std::vector<int> &vParentBarcode,
                  int &vProdNin, int &vProdNout, int &vProdStatus, int &vProdBarcode);
    const PixelID * m_pixelID;
    const SCT_ID  * m_SCT_ID;
    const InDetDD::PixelDetectorManager * m_pixelManager;
    const InDetDD::SCT_DetectorManager  * m_SCT_Manager;    
    mutable int m_event;
    mutable int m_selected;
    std::ofstream * m_outfile;
    ServiceHandle<IPartPropSvc>                 m_particlePropSvc;          
    const HepPDT::ParticleDataTable*            m_particleDataTable;         
    std::string m_name;
    
    int m_offset;
    
    // Accepting truth only if passes these requirements
    float m_max_eta = 4.0;
    float m_min_pt = 1000.;
    int   m_max_barcode = 200e3;
    float m_maxProdVertex = 260.; 

    // options to dump 
    bool m_dumpTruthParticles;
    bool m_dumpClusters;
    bool m_dumpSpacePoints;

    // Information related to the TTree
    //
    std::string m_ntupleFileName;       //!< jobOption: Ntuple file name
    std::string m_ntupleDirName;        //!< jobOption: Ntuple directory name
    std::string m_ntupleTreeName;       //!< jobOption: Ntuple tree name
    int m_maxCL;                        //!< jobOption: maximum number of clusters
    bool m_csvFile;                     //!< jobOption: save data in csv format
    bool m_rootFile;                    //!< jobOption: save data in root format
    TTree* m_nt;

    int m_nCL_ID;
    int *m_CLID;

    int m_nCL;
    int *m_CLindex;
    std::vector<std::string> *m_CLhardware;
    double* m_CLx;
    double *m_CLy;
    double *m_CLz;
    int *m_CLbarrel_endcap;
    int *m_CLlayer_disk;
    int *m_CLeta_module;
    int *m_CLphi_module;
    int *m_CLside;
    uint64_t *m_CLmoduleID;
    std::vector<std::vector<int>> *m_CLparticleLink_eventIndex;
    std::vector<std::vector<int>> *m_CLparticleLink_barcode;
    std::vector<std::vector<bool>> *m_CLbarcodesLinked;
    std::vector<std::vector<int>> *m_CLphis, *m_CLetas, *m_CLtots;
    double *m_CLloc_direction1, *m_CLloc_direction2,*m_CLloc_direction3;
    double *m_CLJan_loc_direction1, *m_CLJan_loc_direction2, *m_CLJan_loc_direction3;
    int *m_CLpixel_count;
    float *m_CLcharge_count;
    float *m_CLloc_eta, *m_CLloc_phi;
    float *m_CLglob_eta, *m_CLglob_phi;
    double *m_CLeta_angle, *m_CLphi_angle;
    float *m_CLnorm_x, *m_CLnorm_y, *m_CLnorm_z;
    std::vector<std::vector<double>> *m_CLlocal_cov;

    int m_nPartEVT;
    int* m_CLevent_number;
    int* m_CLbarcode;
    float *m_CLpx, *m_CLpy, *m_CLpz;
    float* m_CLpt;
    float *m_CLeta;
    float *m_CLvx, *m_CLvy, *m_CLvz;
    float* m_CLradius;
    float* m_CLstatus;
    float* m_CLcharge;
    int* m_CLpdg_id;
    int* m_CLpassed;
    int *m_CLvProdNin, *m_CLvProdNout, *m_CLvProdStatus, *m_CLvProdBarcode;
    std::vector<std::vector<int>> *m_CLvParentID, *m_CLvParentBarcode;

    int m_nSP;
    int* m_SPindex;
    double *m_SPx, *m_SPy, *m_SPz;
    int *m_SPCL1_index, *m_SPCL2_index;

    int m_nTRK;
    int* m_TRKindex;
    int *m_TRKtrack_fitter, *m_TRKparticle_hypothesis;
    std::vector<std::vector<int>> *m_TRKproperties, *m_TRKpattern;
    int *m_TRKndof, *m_TRKmot, *m_TRKoot;
    float* m_TRKchiSq;
    std::vector<std::vector<int>> *m_TRKmeasurementsOnTrack_pixcl_sctcl_index, *m_TRKoutliersOnTrack_pixcl_sctcl_index;
    int* m_TRKcharge;
    std::vector<std::vector<double>> *m_TRKperigee_position, *m_TRKperigee_momentum;
    int *m_TTCindex, *m_TTCevent_index, *m_TTCparticle_link;
    float *m_TTCprobability;

    int m_nDTT;
    int *m_DTTindex, *m_DTTsize;
    std::vector<std::vector<int>> *m_DTTtrajectory_eventindex, *m_DTTtrajectory_barcode, *m_DTTstTruth_subDetType, *m_DTTstTrack_subDetType, *m_DTTstCommon_subDetType;
};

#endif // TRKDUMPALGS_DUMPOBJECTS_H
