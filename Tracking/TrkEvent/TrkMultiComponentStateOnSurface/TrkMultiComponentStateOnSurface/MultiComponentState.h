/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/*********************************************************************************
                        MultiComponentState.h  -  description
                        -------------------------------------
begin                : Friday 31st December 2004
author               : atkinson, amorley,anastopoulos
email                : Anthony.Morley@cern.ch
decription           : Basic definitions for a track state described by more
                       than one set of Track Parameters. The resulting state is
                       a mixture of components. Each component is described by
                       a ComponentParameters object which is of the type
                       std::pair< const TrackParameters*, double>
                       The double describes the weighting of the component -
                       or its relative importance in the mixture.
*********************************************************************************/

#ifndef TrkMultiComponentState
#define TrkMultiComponentState

#include "TrkMultiComponentStateOnSurface/ComponentParameters.h"
#include <vector>

class MsgStream;

namespace Trk {

/** 
 * Prefer SimpleMultiComponentState
 */
class MultiComponentState : public std::vector<ComponentParameters>
{
public:
  /** Default constructor */
  MultiComponentState();
  /** Constructor from a single Component*/
  MultiComponentState(const ComponentParameters&);
  /** destructor */
  ~MultiComponentState();

  /*
   * Since this is practically a vector of pair<pointer,value>
   * disable copy.
   */
  MultiComponentState(const MultiComponentState& other) = delete;
  MultiComponentState& operator=(const MultiComponentState& other) = delete;
  
  MultiComponentState(MultiComponentState&& other) = default;
  MultiComponentState& operator=(MultiComponentState&& other) = default;

  /** Clone method */
  MultiComponentState* clone() const;

  /** Clone with covariance matricies scaled by a factor */
  MultiComponentState* cloneWithScaledError(double) const;

  /** Clone with covariance matrix componants scaled by individual factors
      This will only work if there are 5 track parameters in each componant
  */
  MultiComponentState* cloneWithScaledError(double, double, double, double, double) const;

  /** Check to see if all components in the state have measured track parameters */
  bool isMeasured() const;

  /** Performing renormalisation of total state weighting to one */
  void renormaliseState(double norm = 1);

  /** Dump methods */
  MsgStream& dump(MsgStream&) const;
  std::ostream& dump(std::ostream&) const;
};

/** Overload of << operator for MsgStream and std::ostream */
MsgStream&
operator<<(MsgStream&, const MultiComponentState&);
std::ostream&
operator<<(std::ostream&, const MultiComponentState&);

} // end Trk namespace

inline Trk::MultiComponentState::~MultiComponentState()
{
  Trk::MultiComponentState::const_iterator component = this->begin();
  for (; component != this->end(); ++component) {
    delete component->first;
  }
  this->clear();
}

#endif
