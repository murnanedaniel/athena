/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunction_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunction_h

// STL includes
#include <vector>

#if defined(__FastCaloSimStandAlone__)
  #define ATLAS_NOT_THREAD_SAFE
#else
  #include "CxxUtils/checker_macros.h"
#endif

#include "ISF_FastCaloSimEvent/TFCSFunction.h"

class TH1;

class TFCS1DFunction:public TFCSFunction
{
  public:
    TFCS1DFunction() {};
    ~TFCS1DFunction() {};
    
    virtual int ndim() const {return 1;};

    ///Function gets array of random numbers rnd[] in the range [0,1) as arguments 
    ///and returns function value in array value.
    ///For a n-dimensional function, value and rnd should both have n elements.
    virtual void rnd_to_fct(float value[],const float rnd[]) const;
   
    ///Function gets random number rnd in the range [0,1) as argument and returns function value
    virtual double rnd_to_fct(double rnd) const = 0;

    ///The == operator compares the content of instances. 
    ///The implementation in the base class only returns true for a comparison with itself
    virtual bool operator==(const TFCS1DFunction& ref) const {return this==&ref;};

    static double get_maxdev(TH1* , TH1* );

    static double CheckAndIntegrate1DHistogram(const TH1* hist, std::vector<double>& integral_vec,int& first,int& last);

    static TH1* generate_histogram_random_slope ATLAS_NOT_THREAD_SAFE (int nbinsx=50,double xmin=0,double xmax=1,double zerothreshold=0.1);
    static TH1* generate_histogram_random_gauss ATLAS_NOT_THREAD_SAFE (int nbinsx=50,int ntoy=100000,double xmin=1,double xmax=5,double xpeak=1.5,double sigma=0.6);
    static void unit_test ATLAS_NOT_THREAD_SAFE (TH1* hist,TFCS1DFunction* rtof,int nrnd=1000000,TH1* histfine=nullptr);

  private:

  ClassDef(TFCS1DFunction,2)  //TFCS1DFunction

};

#endif
