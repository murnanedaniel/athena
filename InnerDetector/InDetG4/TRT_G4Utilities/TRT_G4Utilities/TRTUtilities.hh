/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRTUtilities_hh
#define TRTUtilities_hh

#include "globals.hh"
#include <fstream>
#include "AthenaKernel/MsgStreamMember.h"

class TRTOutputFile;


class TRTUtilities
{
public:
    ~TRTUtilities();
  
  static TRTUtilities* GetPointer()
  {
    if (!pUtilities)
      pUtilities = new TRTUtilities;
    return pUtilities;
  }
  
  std::string NumberToString(int);
  
  void PrintTubeParameters(double, double);
  void PrintTubeParameters(double, double, double);
  void PrintTubeParameters(double, double, double, double);

  MsgStream& msg (MSG::Level lvl) const { return m_msg << lvl; }
  bool msgLevel (MSG::Level lvl)    { return m_msg.get().level() <= lvl; }

private:
  TRTUtilities();
  
  std::ofstream& GetReference();

  int printMessages;
  
  TRTOutputFile* pOutputFile;
  
  static TRTUtilities* pUtilities;

  mutable Athena::MsgStreamMember m_msg;

};

#endif
