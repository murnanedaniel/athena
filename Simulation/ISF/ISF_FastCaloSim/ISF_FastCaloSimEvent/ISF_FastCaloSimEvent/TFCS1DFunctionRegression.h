/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegression_h
#define ISF_FASTCALOSIMEVENT_TFCS1DFunctionRegression_h

#include "ISF_FastCaloSimEvent/TFCS1DFunction.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>

using namespace std;

class TFCS1DFunctionRegression:public TFCS1DFunction
{
  public:

    TFCS1DFunctionRegression() {};
    ~TFCS1DFunctionRegression() {};

    virtual double rnd_to_fct(double rnd);
    double  regression_value(double uniform);
    void    set_weights(vector<vector<double> > fWeightMatrix0to1, vector<vector<double> > fWeightMatrix1to2);
    double  sigmoid(double);

  private:

    vector<vector<double> > m_fWeightMatrix0to1;
    vector<vector<double> > m_fWeightMatrix1to2;


  ClassDef(TFCS1DFunctionRegression,1)  //TFCS1DFunctionRegression

};

#endif
