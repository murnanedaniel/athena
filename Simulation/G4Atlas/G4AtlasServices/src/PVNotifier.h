/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef PVNotifier_H
#define PVNotifier_H

// Notifier class to prepend detector name to all G4 Physical Volumes
// Only to be used by the G4GeometryNotifierSvc

#include "G4VNotifier.hh"

class G4GeometryNotifierSvc;

  class PVNotifier: public G4VNotifier {
  friend class G4GeometryNotifierSvc;
  private:
    PVNotifier(G4GeometryNotifierSvc*);

    G4GeometryNotifierSvc* m_notifierSvc;
  public:
    void NotifyRegistration();
    void NotifyDeRegistration();
  };

#endif
