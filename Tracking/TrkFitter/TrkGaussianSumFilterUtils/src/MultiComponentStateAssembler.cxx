/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file   MultiComponentStateAssembler.cxx
 * @date   Monday 20th December 2004
 * @author Atkinson,Anthony Morley, Christos Anastopoulos
 *
 * Implementation code for MultiComponentStateAssembler
 */

#include "TrkGaussianSumFilterUtils/MultiComponentStateAssembler.h"
#include "TrkGaussianSumFilterUtils/ComponentParameters.h"
//
#include <limits>

namespace {

using namespace Trk::MultiComponentStateAssembler;

/** @brief Helper for ordering by larger to smaller weight*/
class SortByLargerSimpleComponentWeight
{
public:
  SortByLargerSimpleComponentWeight() = default;
  bool operator()(const Trk::ComponentParameters& firstComponent,
                  const Trk::ComponentParameters& secondComponent) const
  {
    return firstComponent.second > secondComponent.second;
  }
};

/** @bried Method to check the validity of of the cached state */
inline bool
isStateValid(const Cache& cache)
{
  return !cache.multiComponentState.empty();
}

/** @bried Method to assemble state with correct weightings */
Trk::MultiComponentState
doStateAssembly(Cache&& cache, const double newWeight)
{
  if (!isStateValid(cache)) {
    return {};
  }
  const size_t cacheSize = cache.multiComponentState.size();
  if (cache.validWeightSum <= 0.) {
    if (!cache.multiComponentState.empty()) {
      const double fixedWeights = 1. / static_cast<double>(cacheSize);
      for (auto& component : cache.multiComponentState) {
        component.second = fixedWeights;
      }
    }
    Trk::MultiComponentState assembledState =
      std::move(cache.multiComponentState);
    // Reset the cache before leaving
    return assembledState;
  }

  Trk::MultiComponentState assembledState =
    std::move(cache.multiComponentState);
  const double scalingFactor =
    cache.validWeightSum > 0. ? newWeight / cache.validWeightSum : 1.;
  for (auto& component : assembledState) {
    component.second *= scalingFactor;
  }
  // Reset the cache before leaving
  return assembledState;
}

/** @bried Method to Check component entries before full assembly */
bool
prepareStateForAssembly(Cache& cache)
{
  // Protect against empty state
  if (!isStateValid(cache)) {
    return false;
  }

  // Check for minimum fraction of valid states
  double den = cache.validWeightSum + cache.invalidWeightSum;
  double validWeightFraction = den > 0 ? cache.validWeightSum / den : 0;
  if (cache.invalidWeightSum > 0. &&
      validWeightFraction < cache.minimumValidFraction) {
    return false;
  }
  // Check to see assembly has not already been done
  if (cache.assemblyDone) {
    return true;
  }
  // Sort Multi-Component State by weights
  std::sort(cache.multiComponentState.begin(),
            cache.multiComponentState.end(),
            SortByLargerSimpleComponentWeight());

  double totalWeight(cache.validWeightSum + cache.invalidWeightSum);
  if (totalWeight != 0.) {

    // ordered in descending order
    // return the 1st element where (element<value)

    const double minimumWeight =
      std::max(cache.minimumFractionalWeight * totalWeight,
               std::numeric_limits<double>::min());

    const Trk::ComponentParameters dummySmallestWeight(nullptr, minimumWeight);

    auto lower_than = std::upper_bound(cache.multiComponentState.begin(),
                                       cache.multiComponentState.end(),
                                       dummySmallestWeight,
                                       SortByLargerSimpleComponentWeight());

    // reverse iterate , so as to delete removing the last
    auto lower_than_reverse = std::make_reverse_iterator(lower_than);
    for (auto itr = cache.multiComponentState.rbegin();
         itr != lower_than_reverse;
         ++itr) {
      cache.multiComponentState.erase(itr.base() - 1);
    }
  }
  // Now recheck to make sure the state is now still valid
  if (!isStateValid(cache)) {
    return false;
  }

  // Set assembly flag
  cache.assemblyDone = true;
  return true;
}

} // end anonymous namespace

void
Trk::MultiComponentStateAssembler::reset(Cache& cache)
{
  cache.assemblyDone = false;
  if (!cache.multiComponentState.empty()) {
    cache.multiComponentState.clear();
  }
  cache.validWeightSum = 0.;
  cache.invalidWeightSum = 0.;
}

bool
Trk::MultiComponentStateAssembler::addComponent(
  Cache& cache,
  ComponentParameters&& componentParameters)
{
  if (cache.assemblyDone) {
    return false;
  }
  cache.validWeightSum += componentParameters.second;
  cache.multiComponentState.emplace_back(std::move(componentParameters.first),
                                         componentParameters.second);
  return true;
}

bool
Trk::MultiComponentStateAssembler::addMultiState(
  Cache& cache,
  Trk::MultiComponentState&& multiComponentState)
{
  if (cache.assemblyDone) {
    return false;
  }
  double sumW(0.);
  for (auto& component : multiComponentState) {
    sumW += component.second;
    cache.multiComponentState.emplace_back(std::move(component.first),
                                           component.second);
  }
  multiComponentState.clear();
  cache.validWeightSum += sumW;
  return true;
}

Trk::MultiComponentState
Trk::MultiComponentStateAssembler::assembledState(
  MultiComponentStateAssembler::Cache&& cache)
{
  if (!prepareStateForAssembly(cache)) {
    return {};
  }
  double totalWeight = cache.validWeightSum;
  if (cache.invalidWeightSum > 0. || cache.validWeightSum <= 0.) {
    totalWeight = cache.validWeightSum + cache.invalidWeightSum;
  }

  return doStateAssembly(std::move(cache), totalWeight);
}

