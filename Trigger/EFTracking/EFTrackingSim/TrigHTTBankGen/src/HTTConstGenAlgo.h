// Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#ifndef HTTConstGenAlgo_h
#define HTTConstGenAlgo_h

/**
 * @file HTTConstGenAlgo.h
 * @author Unknown; major rewrite Riley Xu - riley.xu@cern.ch
 * @date May 28th, 2020
 * @brief Algorithm to generate fit constants.
 *
 * This algorithm inputs matrix files generated by HTTMatrixGenAlgo,
 * and outputs ROOT/text files containing the fit constants used by HTT to fit
 * tracks. Each sector contains the information in geo_constants struct declared
 * below.
 *
 * See ATLAS-TDR-FTK-021 section 5.2.2. for a detailed description of constant generation.
 * The constants are calculated via equation 9,
 *
 *      C_ij = Inv(V_jk) . ( <x_k p_i> - <x_k> <p_i> )
 *
 * where x_k are the hit coordinates, p_i the helix parameters, and V_jk the covariance matrix
 *
 *      V_jk = <x_j x_k> - <x_j> <x_k>
 */

#include "GaudiKernel/ITHistSvc.h"
#include "AthenaBaseComps/AthAlgorithm.h"
#include "TrigHTTObjects/HTTVectors.h"
#include "TrigHTTMaps/ITrigHTTMappingSvc.h"
#include "TrigHTTMaps/HTTPlaneMap.h"
#include "TTree.h"
#include "HTTMatrixAccumulator.h"

#include <string>
#include <vector>

#include <TMatrixD.h>

#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"

class ITrigHTTMappingSvc;


// data structure that contain definition of geometrical constants for linear fits
struct geo_constants
{
    //these are the matrices themselves
    std::vector<double> Vd0; // impact paramers coefs
    std::vector<double> Vcurvature; // curvature coefs
    std::vector<double> Vphi; // phi coefs
    std::vector<double> Vz0; // z0 coefs
    std::vector<double> Veta; // eta coefs
    vector2D<double>    kernel; // kernel (covariance matrix). Size (ndim - npar, ndim)
    std::vector<double> kaverages; // averages, useful to keep track of

    // these are the constant/offset values
    HTTTrackPars pars;

    int real; // this value is greater than 0 if these constants are correctly evaulated

    geo_constants(size_t nCoords) :
        Vd0(nCoords),
        Vcurvature(nCoords),
        Vphi(nCoords),
        Vz0(nCoords),
        Veta(nCoords),
        kernel(nCoords - HTTTrackPars::NPARS, nCoords),
        kaverages(nCoords - HTTTrackPars::NPARS),
        pars(0),
        real(0)
    {
    }
};


class HTTConstGenAlgo : public AthAlgorithm
{
    public:

        HTTConstGenAlgo(const std::string& name, ISvcLocator* pSvcLocator);
        virtual ~HTTConstGenAlgo() = default;
        StatusCode initialize() override;
        StatusCode finalize() override;

	// Execute does not do anything for this alg. This class does not process events, everything is done in initalize
	StatusCode execute() override;

        StatusCode bookHistograms();

    private:

        ///////////////////////////////////////////////////////////////////////
        // Handles

        ServiceHandle<ITrigHTTMappingSvc> m_HTTMapping{this, "TrigHTTMappingSvc","TrigHTTMappingSvc"};
        ServiceHandle<ITHistSvc> m_tHistSvc{this, "THistSvc","THistSvc"};

        const HTTPlaneMap* m_pmap = nullptr;


	
        ///////////////////////////////////////////////////////////////////////
        // Configuration
	Gaudi::Property<std::string> m_cfpath{this, "merged_file_path", "", "merged file"};
	Gaudi::Property<std::string> m_skipFile{this, "skip_sectors", "File with list of sectors to skip"};
	Gaudi::Property<bool> m_Monitor{this,"Monitor",false,"flag to enable the monitor"};
	Gaudi::Property<int> m_region{this, "region",0,"region to run"};
	Gaudi::Property<bool> m_CheckGood2ndStage{this,"CheckGood2ndStage","Check goodness of 2nd stage fit constants?"};    
	Gaudi::Property<bool> m_useHitScaleFactor{this,"UseHitScaleFactor",false,"Scale factor for hits"};
	Gaudi::Property<bool> m_isSecondStage{this,"IsSecondStage",false,"If false, we're doing a 1st stage fit, otherwise 2nd stage"};
	Gaudi::Property<bool> m_dumpMissingHitsConstants{this, "missHitsConsts", false, "if this is true we dump constants assuming a missing hit in each layer, too"};

        ///////////////////////////////////////////////////////////////////////
        // ROOT Objects

        TFile *m_mafile = nullptr;
        TTree *m_ctree = nullptr;
        TTree *m_matrix_tree = nullptr;
        TTree *m_good_tree = nullptr;

        TH1F *m_h_vc = nullptr;
        TH1F *m_h_vd = nullptr;
        TH1F *m_h_vf = nullptr;
        TH1F *m_h_vz = nullptr;
        TH1F *m_h_veta = nullptr;

        ///////////////////////////////////////////////////////////////////////
        // Slice Info

        HTTTrackParsI m_sliceNBins;
        HTTTrackPars m_sliceMin;
        HTTTrackPars m_sliceMax;

        ///////////////////////////////////////////////////////////////////////
        // Sizes

        int m_nLayers = 0;
        int m_nKernel = 0;
        int m_nKAverages = 0;
        int m_nCoords = 0;
        int m_nCoords_2 = 0; // m_nCoords^2


        ///////////////////////////////////////////////////////////////////////
        // Main Storage Objects

        // These have size = # of good sectors.
        std::vector<geo_constants> m_geo_consts;

        // These are the constants for missing hits, first index is missing hit
        // The second index is the same as above (good sector number)
        std::vector<std::vector<geo_constants>> m_geo_consts_with_missinghit;

        // Size = # of sectors. Which sectors to skip for generating constants.
        std::vector<bool> m_skipList;

        ///////////////////////////////////////////////////////////////////////
        // Helper Functions

        StatusCode copySliceTree(TFile *file);
        StatusCode prepareOutputTree();
        void readSkipList(size_t nEntries);
        void generate_constants();
        void fillConstTree(std::vector<module_t> & modules, HTTMatrixAccumulator & acc, geo_constants & geo);
        bool isNAN(double value, const char* name);
        bool failedConstants(geo_constants const & geo, std::vector<bool> const & usable);
        void DumpConstants(std::vector<geo_constants> &geo_consts, std::string & filename);
        void writeSectors();
        bool GetConstants(HTTMatrixAccumulator const &acc_norm, geo_constants &geo, int entryNumber); // use values in acc_norm
        bool GetConstants(HTTMatrixAccumulator const &acc_norm, geo_constants &geo, int entryNumber, std::vector<bool> const &coordsToUse, unsigned int nusable); // full method with different number of usable coordinates
        void createMissingHitsConstants(HTTMatrixAccumulator const & acc_norm, size_t entry);
	HTTMatrixAccumulator normalize(HTTMatrixAccumulator const & acc_raw);
	geo_constants makeConsts(HTTMatrixAccumulator const & acc, std::vector<bool> const & usable,
				 std::vector<double> const & inv_covariance,
				 std::vector<double> const & eigvals, vector2D<double> const & eigvecs);
	std::vector<double> matrix_multiply(std::vector<double> const & A, std::vector<double> const & b);
	void eigen(size_t n_redu, size_t n_full, TMatrixD &mtx, std::vector<bool> const & usable, std::vector<double> & eigvals_v, vector2D<double> & eigvecs_v);
	std::vector<double> invert(size_t n_full, TMatrixD mtx, std::vector<bool> const & usable);
	TMatrixD getReducedMatrix(size_t n, std::vector<double> const & mtx_v, std::vector<bool> const & usable, size_t nDimToUse);
	bool isSingular(TMatrixD mtx);
	double dot(const double* vec1, const double* vec2, size_t size);
	geo_constants calculate_gcorth(geo_constants geo, int nCoords, std::vector<bool> const & usable);



};

#endif // HTTConstGenAlgo_h

