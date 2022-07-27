/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef G4UserActions_LengthIntegrator_H
#define G4UserActions_LengthIntegrator_H

#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ServiceHandle.h"

#include "G4Pow.hh"
#include "TString.h"

#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"

#include <string>
#include <map>

// Forward declarations
class TProfile;
class TProfile2D;


namespace G4UA
{

  /// @class LengthIntegrator
  /// @brief A user action used to evaluate thickness of all detectors
  ///        traversed by outgoing particles.
  ///
  /// This user action is currently used only in special runs with geantinos.
  /// Thickness is recorded in terms of both rad length and int length.
  ///
  /// NOTE: the current design is safe for multi-threading, but _not_
  /// performant due to sharing of the histograms and excessive locking. If
  /// this action needs to be used in multi-threaded jobs, we can rewrite it so
  /// that each instance has its own copy of the histograms which get merged in
  /// finalization of the LengthIntegratorTool.
  ///
  class LengthIntegrator final : public G4UserEventAction,
                                 public G4UserSteppingAction
  {

    public:

      /// Constructor takes the name of the histogram service as argument.
    LengthIntegrator(const std::string& histSvcName, bool doHistos);

      /// Called at beginning of G4 event to cache some details about the
      /// current primary vertex and particle. Also resets some measurements.
      virtual void BeginOfEventAction(const G4Event*) override;

      /// Called at end of G4 event to finalize measurements and fill hists
      virtual void EndOfEventAction(const G4Event*) override;

      /// Called at every particle step to accumulate thickness.
      virtual void UserSteppingAction(const G4Step*) override;

    private:

      // Holder for G4 math tools
      G4Pow* m_g4pow;

      //Tree for material information
      TTree* m_tree;

      //Tree Branches
      int   m_genNPart = 0;
      float m_genEta = 0.0F;
      float m_genPhi = 0.0F;
      float m_genZ = 0.0F;
      float m_genR = 0.0F;
      
      //X0 Branches
      float m_total_X0 = 0.0F;
      float m_total_L0 = 0.0F;

      std::vector<double> m_collected_X0;
      std::vector<double> m_collected_L0;

      std::vector<float> m_collected_hitr;
      std::vector<float> m_collected_hitz;

      std::vector<float> m_collected_outhitr;
      std::vector<float> m_collected_outhitz;

      std::vector<float> m_collected_density;
      std::vector<std::string> m_collected_material;
      std::vector<std::string> m_collected_volume;
      
      std::vector<std::string> m_collected_groupedmaterial;
      std::vector<std::string> m_collected_volumetype;

      bool m_splitModerator;
      bool m_splitPP1;  

      void fillNtuple();
      std::string getMaterialClassification(const std::string& name);
      std::string getVolumeType(const std::string& s);

      // Add elements and values into the map
      void addToDetThickMap(const std::string&, double, double);

      /// Setup one set of measurement hists for a detector name.
      void regAndFillHist(const std::string&, const std::pair<double, double>&);

      /// this method checks if a histo is on THsvc already and caches a local pointer to it
      /// if the histo is not present, it creates and registers it
      TProfile2D* getOrCreateProfile(const std::string& regName, const TString& histoname, const TString& xtitle, int nbinsx, float xmin, float xmax,const TString& ytitle, int nbinsy,float ymin, float ymax,const TString& ztitle);

      /// Handle to the histogram service
      ServiceHandle<ITHistSvc> m_hSvc;
      
      //Do we create histograms
      bool m_doHistos;

      /// Cached eta of the current primary
      double m_etaPrimary;
      /// Cached phi of the current primary
      double m_phiPrimary;

      /// Map of detector thickness measurements for current event
      std::map<std::string, std::pair<double, double> > m_detThickMap;

      /// Rad-length profile hist in R-Z
      TProfile2D* m_rzProfRL;
      /// Rad-length profile hist in eta
      std::map<std::string, TProfile*> m_etaMapRL;
      /// Rad-length profile hist in phi
      std::map<std::string, TProfile*> m_phiMapRL;



      /// Int-length profile hist in R-Z
      TProfile2D* m_rzProfIL;
      /// Int-length profile hist in eta
      std::map<std::string, TProfile*> m_etaMapIL;
      /// Int-length profile hist in phi
      std::map<std::string, TProfile*> m_phiMapIL;

      // 2D plots of rad-length and int-length
      std::map<std::string,TProfile2D*,std::less<std::string> > m_rzMapRL;
      std::map<std::string,TProfile2D*,std::less<std::string> > m_xyMapRL;

      std::map<std::string,TProfile2D*,std::less<std::string> > m_rzMapIL;
      std::map<std::string,TProfile2D*,std::less<std::string> > m_xyMapIL;

  }; // class LengthIntegrator

} // namespace G4UA

#endif
