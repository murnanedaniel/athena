/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file   GsfBetheHeitlerEffects.h
 * @date   Friday 18th February 2005
 * @author Tom Athkinson, Anthony Morley, Christos Anastopoulos
 * @brief  Bethe-Heitler material effects for electrons, used in the GSF
 */

#ifndef Trk_GsfBetheHeitlerEffects_H
#define Trk_GsfBetheHeitlerEffects_H

#include "TrkGaussianSumFilter/IBetheHeitlerEffects.h"

#include "AthenaBaseComps/AthAlgTool.h"
#include "TrkEventPrimitives/ParticleHypothesis.h"
#include "TrkEventPrimitives/PropDirection.h"
#include "TrkMultiComponentStateOnSurface/MultiComponentState.h"

namespace Trk {

class GsfBetheHeitlerEffects
  : public AthAlgTool
  , virtual public IBetheHeitlerEffects
{

public:
  /// Maximum number of Gaussian components in the mixture
  static constexpr int maxNumberofComponents = 8;
  /// Polynomail coefficients , order is this -1 
  static constexpr int polynomialCoefficients = 6;
  struct ComponentValues
  {
    // Default ctors/dtor/assignment operators
    ComponentValues() = default;
    ~ComponentValues() = default;
    ComponentValues(const ComponentValues&) = default;
    ComponentValues& operator=(const ComponentValues&) = default;
    ComponentValues(ComponentValues&&) = default;
    ComponentValues& operator=(ComponentValues&&) = default;
    // Constructor with arguments
    ComponentValues(double aWeight, double aMean, double aVariance)
      : weight(aWeight)
      , mean(aMean)
      , variance(aVariance)
    {}
    double weight;
    double mean;
    double variance;
  };

  using MixtureParameters = std::array<ComponentValues, maxNumberofComponents>;

  /** Helper class for construction and evaluation of polynomial */
  class Polynomial
  {
  public:
    Polynomial() = default;
    ~Polynomial() = default;
    Polynomial(const Polynomial&) = default;
    Polynomial& operator=(const Polynomial&) = default;
    Polynomial(Polynomial&&) = default;
    Polynomial& operator=(Polynomial&&) = default;
    
    /** Constructor from a vector of coefficients 
     * (in decreasing order of powers of x) */
    Polynomial(const std::vector<double>& coefficients)
      : m_coefficients(coefficients){};

    // Evaluation of the polynomial for given material thickness (t)
    double operator()(const double& t) const
    {
      double sum(0.);
      std::vector<double>::const_iterator coefficient = m_coefficients.begin();

      for (; coefficient != m_coefficients.end(); ++coefficient) {
        sum = t * sum + (*coefficient);
      }

      return sum;
    }

  private:
    std::vector<double> m_coefficients;
  };

  GsfBetheHeitlerEffects(const std::string&,
                         const std::string&,
                         const IInterface*);

  virtual ~GsfBetheHeitlerEffects() = default;

  /** AlgTool initialise method */
  virtual StatusCode initialize() override final;

  virtual void compute(Trk::GSFEnergyLossCache& cache,
                       const ComponentParameters& componentParameters,
                       const MaterialProperties& materialProperties,
                       double pathLenght,
                       PropDirection direction = anyDirection,
                       ParticleHypothesis particleHypothesis =
                         nonInteracting) const override final;

private:
  // Read polynomial fit parameters from a specified file
  bool readParameters();

  // Read coeffients for a single polynomial fit
  Polynomial readPolynomial(std::ifstream&, const int);

  std::array<Polynomial,maxNumberofComponents> m_polynomialWeights;
  std::array<Polynomial,maxNumberofComponents> m_polynomialMeans;
  std::array<Polynomial,maxNumberofComponents> m_polynomialVariances;
  std::array<Polynomial,maxNumberofComponents> m_polynomialWeightsHighX0;
  std::array<Polynomial,maxNumberofComponents> m_polynomialMeansHighX0;
  std::array<Polynomial,maxNumberofComponents> m_polynomialVariancesHighX0;

  int m_numberOfComponents;
  int m_transformationCode;
  int m_correctionFlag;
  int m_numberOfComponentsHighX0;
  int m_transformationCodeHighX0;

  double m_singleGaussianRange;
  double m_lowerRange;
  double m_xOverRange;
  double m_upperRange;
  double m_componentMeanCut;

  bool m_useHighX0;
  std::string m_parameterisationFileName;
  std::string m_parameterisationFileNameHighX0;
};

}

#endif
