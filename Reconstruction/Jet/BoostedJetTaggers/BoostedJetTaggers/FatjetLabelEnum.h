#ifndef __FATJETLABELENUM_H
#define __FATJETLABELENUM_H

namespace FatjetTruthLabel
{
  enum TypeEnum
  {
    UNKNOWN=0, // Not tagged yet
    tqqb,      // full-contained top->qqb
    Wqq,       // full-contained W->qq
    Zqq,       // full-contained Z->qq
    Wqq_From_t,// full-contained W->qq (also mathced to top)
    other_From_t, // matched to top
    other_From_V, // matched to W/Z
    notruth,   // failed to truth-jet matching (pileup)
    qcd,       // not matched to top or W/Z (background jet)
  };  
  inline int enumToInt(const TypeEnum type)
  {
    switch (type)
      {
      case tqqb:         return 1;
      case Wqq:          return 2;
      case Zqq:          return 3;
      case Wqq_From_t:   return 4;
      case other_From_t: return 5;
      case other_From_V: return 6;
      case notruth:      return 7;
      case qcd:          return 8;
      default:           return 0;
      }
  }  
  inline TypeEnum intToEnum(const int type)
  {
    if ( type==1 ){
      return tqqb;
    }else if ( type==2 ){
      return Wqq;
    }else if ( type==3 ){
	return Zqq;
    }else if ( type==4 ){
      return Wqq_From_t;
    }else if ( type==5 ){
      return other_From_t;
    }else if ( type==6 ){
      return other_From_V;
    }else if ( type==7 ){
      return notruth;
    }else if ( type==8 ){
      return qcd;
    }
    
    return UNKNOWN;
  }
  inline TypeEnum stringToEnum(const TString& name)
  {
    if (name.EqualTo("tqqb",TString::kIgnoreCase))
        return tqqb;
    if (name.EqualTo("Wqq",TString::kIgnoreCase))
        return Wqq;
    if (name.EqualTo("Zqq",TString::kIgnoreCase))
        return Zqq;
    if (name.EqualTo("Wqq_From_t",TString::kIgnoreCase))
        return Wqq_From_t;
    if (name.EqualTo("other_From_t",TString::kIgnoreCase))
        return other_From_t;
    if (name.EqualTo("other_From_V",TString::kIgnoreCase))
        return other_From_V;
    if (name.EqualTo("notruth",TString::kIgnoreCase))
        return notruth;
    if (name.EqualTo("qcd",TString::kIgnoreCase))
        return qcd;
    return UNKNOWN;
  }
}


namespace FatjetCutResult
{
  enum TypeEnum
  {
    UNKNOWN=0, // not tagged yet
    passMpassD2_2Var,
    passMfailD2_2Var,
    failMpassD2_2Var,
    failMfailD2_2Var
  };
  inline int enumToInt(const TypeEnum type)
  {
    switch (type)
      {
      case passMpassD2_2Var: return 1;
      case passMfailD2_2Var: return 2;
      case failMpassD2_2Var: return 3;
      case failMfailD2_2Var: return 4;
      default:               return 0;
      }
  }
  inline TypeEnum intToEnum(const int type)
  {
    if ( type==1 ){
      return passMpassD2_2Var;
    }else if ( type==2 ){
      return passMfailD2_2Var;
    }else if ( type==3 ){
      return failMpassD2_2Var;
    }else if ( type==4 ){
      return failMfailD2_2Var;
    }
    
    return UNKNOWN;
  }
}

#endif
