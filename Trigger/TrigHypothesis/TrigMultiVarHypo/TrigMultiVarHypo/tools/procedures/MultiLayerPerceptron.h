/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/



#ifndef TRIGMULTIVARHYPO_MULTILAYERPERCEPTRON_H
#define TRIGMULTIVARHYPO_MULTILAYERPERCEPTRON_H

#include <string>
#include <vector>
#include "TrigMultiVarHypo/tools/procedures/IModel.h"


namespace Ringer{

  class MultiLayerPerceptron: public IModel {
  
    private:
  
  
      std::vector<unsigned int>   m_nodes; 
      double                   ***m_weights;
      double                    **m_bias;
      std::vector<std::string>    m_tfnames;
  
    
    public:
  
      MultiLayerPerceptron( std::vector<unsigned int>  &, 
                            std::vector<double>        &, 
                            std::vector<double>        &,
                            std::vector<std::string>   &,
                            double etmin,  
                            double etmax, 
                            double etamin,  
                            double etamax,
                            double mumin,   
                            double mumax);
  
      ~MultiLayerPerceptron();

      // no copy
      MultiLayerPerceptron(const MultiLayerPerceptron&) = delete;
      MultiLayerPerceptron& operator=(const MultiLayerPerceptron&) = delete;

      Ringer::RnnOutInfo propagate(std::vector<float> &input) const;
  };

}
#endif
