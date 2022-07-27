/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
//   Implementation file for class SiDetElementLink_xk
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////
// Version 1.0 21/04/2004 I.Gavrilenko
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include "SiDetElementsRoadTool_xk/SiDetElementLink_xk.h"


///////////////////////////////////////////////////////////////////
// Set parameters
///////////////////////////////////////////////////////////////////

void InDet::SiDetElementLink_xk::set(const double* P, bool isITk)
{
  m_z            = float(P[ 1])                        ; //  Z
  m_dz           = float(std::abs(P[12]-P[11])*.5)         ; // dZ
  m_phi          = float(P[ 2])                        ; // azimuthal angle
  m_geo   [0]    = float(P[ 3])                        ; // min. distance
  m_geo   [1]    = float(P[ 4])                        ; // phi
  m_geo   [2]    = float(P[ 5])                        ; // sin(phi)
  m_geo   [3]    = float(P[ 6])                        ; // cos(phi)
  m_geo   [4]    = float(P[ 7])                        ; // sin(polar)
  m_geo   [5]    = float(P[ 8])                        ; // cos(polar)
  m_center[0]    = float(P[18])                        ; //
  m_center[1]    = float(P[19])                        ; //
  m_bound [0][2] = float(sqrt(P[20]*P[20]+P[21]*P[21])); // -F
  m_bound [0][0] = float(P[20]/double(m_bound[0][2]))  ; // 
  m_bound [0][1] = float(P[21]/double(m_bound[0][2]))  ; //
  m_bound [1][2] = float(sqrt(P[22]*P[22]+P[23]*P[23])); // +ZR
  m_bound [1][0] = float(P[22]/double(m_bound[1][2]))  ; //
  m_bound [1][1] = float(P[23]/double(m_bound[1][2]))  ; //
  m_bound [2][2] = float(sqrt(P[24]*P[24]+P[25]*P[25])); // +F
  m_bound [2][0] = float(P[24]/double(m_bound[2][2]))  ; //
  m_bound [2][1] = float(P[25]/double(m_bound[2][2]))  ; //
  m_bound [3][2] = float(sqrt(P[26]*P[26]+P[27]*P[27])); // -ZR
  m_bound [3][0] = float(P[26]/double(m_bound[3][2]))  ; //
  m_bound [3][1] = float(P[27]/double(m_bound[3][2]))  ; //
  if (isITk) {
    m_bound [0][2]+=P[40];
    m_bound [1][2]+=P[40];
    m_bound [2][2]+=P[40];
    m_bound [3][2]+=P[40];
  }
}

///////////////////////////////////////////////////////////////////
// Detector element intersection using cashed information
// Input  parameters: r[0] - X    a[0] - Ax 
//                    r[1] - Y    a[1] - Ay
//                    r[2] - Z    a[2] - Az
// Output parameters: O[0] - close distance in azimuthal direction
//                    O[1] - close distance in r or z    direction
//                    O[2] - step to detector element 
///////////////////////////////////////////////////////////////////

void InDet::SiDetElementLink_xk::intersect
(const float* r ,const float* a,float* O) const
{
  const float* g = &m_geo[0];
  float     s[3] = {g[3]*g[4],g[2]*g[4],g[5]}                         ;
  float     S    = a[0]*s[0]+a[1]*s[1]+a[2]*s[2]                      ;
  if(S!=0.) S    = (g[0]-(r[0]*s[0]+r[1]*s[1]+r[2]*s[2]))/S           ;
  float    rn[3] = {r[0]+S*a[0],r[1]+S*a[1],r[2]+S*a[2]}              ;
  float  d0      =       rn[1]*g[3]-rn[0]*g[2]            -m_center[0];
  float  d1      = g[5]*(rn[0]*g[3]+rn[1]*g[2])-g[4]*rn[2]-m_center[1];
  O[0]           = m_bound[0][0]*d0+m_bound[0][1]*d1-m_bound[0][2]    ;
  float  daz     = m_bound[2][0]*d0+m_bound[2][1]*d1-m_bound[2][2]    ; 
  O[1]           = m_bound[1][0]*d0+m_bound[1][1]*d1-m_bound[1][2]    ;
  float  drz     = m_bound[3][0]*d0+m_bound[3][1]*d1-m_bound[3][2]    ;
  if(O[0] < daz) O[0] = daz; 
  if(O[1] < drz) O[1] = drz;
  O[2] = S;
}


///////////////////////////////////////////////////////////////////
// Detector element intersection using cashed information
// Input  parameters: r[0] - X    a[0] - Ax
//                    r[1] - Y    a[1] - Ay
//                    r[2] - Z    a[2] - Az
// Output parameters: Step - step to detector element
///////////////////////////////////////////////////////////////////
bool InDet::SiDetElementLink_xk::intersectITk
(const float* r,const float* a,float& Step) const
{
  const float* g = &m_geo[0];

  float     S    = a[0]*g[6]+a[1]*g[7]+a[2]*g[5]                      ;
  if(S!=0.) S    = (g[0]-(r[0]*g[6]+r[1]*g[7]+r[2]*g[5]))/S           ;
  float r0 = r[0]+S*a[0]                          ;
  float r1 = r[1]+S*a[1]                          ;
  float r2 = r[2]+S*a[2]                          ;
  float d0 =       r1*g[3]-r0*g[2]          -g[8] ;
  float d1 = g[5]*(r0*g[3]+r1*g[2])-(g[4]*r2+g[9]);

  if((m_bound[1][0]*d0+m_bound[1][1]*d1) > m_bound[1][2] ||
     (m_bound[3][0]*d0+m_bound[3][1]*d1) > m_bound[3][2] ||
     (m_bound[0][0]*d0+m_bound[0][1]*d1) > m_bound[0][2] ||
     (m_bound[2][0]*d0+m_bound[2][1]*d1) > m_bound[2][2]) return false;
  Step = S+r[5];
  return true;
}


