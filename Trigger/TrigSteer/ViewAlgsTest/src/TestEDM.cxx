/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "./TestEDM.h"
namespace {
  void set(TestCluster* c, const char* n, float x ) {
    c->setDetail(n, x);
  }
  float get(const TestCluster* c, const char* n ) {
    float x{};
    c->getDetail(n, x);
    return x;
  }

}
namespace TestEDM {


  void setClusterEt(TestCluster* c, float et) {
    set(c, "et", et);
  }

  float getClusterEt(const TestCluster* c) {
    return get(c, "et");
  }
  
  void setClusterEta(TestCluster* c, float eta) {
    set(c, "eta", eta);
  }
  float getClusterEta(const TestCluster* c) {
    return get(c, "eta");
  }
  
  void setClusterPhi(TestCluster* c, float phi) {
    set(c, "phi", phi); 
  }
  float getClusterPhi(const TestCluster* c) {
    return get(c, "phi");
  }


  void setTrackPt(TestTrack* c, float pt) {
    set(c, "pt", pt);
  }
  float getTrackPt(const TestTrack* c) {
    return get(c, "pt");  
  }
  
  void setTrackEta(TestTrack* c, float eta) {
    set(c, "eta", eta);
  }
  float getTrackEta(const TestTrack* c) {
    return get(c, "eta");  
  }
  
  void setTrackPhi(TestTrack* c, float phi) {
    set(c, "phi", phi);
  }
  float getTrackPhi(const TestTrack* c) {
    return get(c, "phi");  
  }



}
