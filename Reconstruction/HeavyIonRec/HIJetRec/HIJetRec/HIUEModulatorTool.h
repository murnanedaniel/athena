/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef __HIJETREC_HIUEMODULATORTOOL_H__
#define __HIJETREC_HIUEMODULATORTOOL_H__

#include "HIJetRec/IHIUEModulatorTool.h"
#include "AsgTools/AsgTool.h"

class HIUEModulatorTool : public asg::AsgTool, virtual public IHIUEModulatorTool
{
  ASG_TOOL_CLASS(HIUEModulatorTool,IHIUEModulatorTool)
public:

  HIUEModulatorTool(const std::string& myname);

  float getModulation(float phi) const;
  StatusCode setEventShapeForModulation(const xAOD::HIEventShape* shape);
  void setHarmonics(const std::vector<unsigned int>& v);
  StatusCode retrieveShape();
  StatusCode initialize();

  static float modulate(const std::vector<unsigned int>& nh_vector, const xAOD::HIEventShape* shape, float phi);
  
private:
  StatusCode checkQVectorSize(const xAOD::HIEventShape* shape, unsigned int n) const;
  StatusCode checkCompatibility() const;
  std::string m_shape_key;
  const xAOD::HIEventShape* m_shape;
  std::vector<unsigned int> m_nh_vector;
  bool m_do_v2;
  bool m_do_v3;
  bool m_do_v4;


};
#endif
