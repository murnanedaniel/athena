/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include <Riostream.h> 

#include "BeamSpotPdf.h" 
#include <RooAbsReal.h> 
#include <RooAbsCategory.h> 
#include <math.h> 
#include <TMath.h> 



 BeamSpotPdf::BeamSpotPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _y,
                        RooAbsReal& _z,
                        RooAbsReal& _vxx,
                        RooAbsReal& _vyy,
                        RooAbsReal& _vxy,
                        RooAbsReal& _mx,
                        RooAbsReal& _sx,
                        RooAbsReal& _ax,
                        RooAbsReal& _my,
                        RooAbsReal& _sy,
                        RooAbsReal& _ay,
                        RooAbsReal& _mz,
                        RooAbsReal& _sz,
                        RooAbsReal& _k,
			  RooAbsReal& _rho
			 ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   y("y","y",this,_y),
   z("z","z",this,_z),
   vxx("vxx","vxx",this,_vxx),
   vyy("vyy","vyy",this,_vyy),
   vxy("vxy","vxy",this,_vxy),
   mx("mx","mx",this,_mx),
   sx("sx","sx",this,_sx),
   ax("ax","ax",this,_ax),
   my("my","my",this,_my),
   sy("sy","sy",this,_sy),
   ay("ay","ay",this,_ay),
   mz("mz","mz",this,_mz),
   sz("sz","sz",this,_sz),
   k("k","k",this,_k),
  rho("rho","rho",this,_rho)
 { 
 } 


 BeamSpotPdf::BeamSpotPdf(const BeamSpotPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   y("y",this,other.y),
   z("z",this,other.z),
   vxx("vxx",this,other.vxx),
   vyy("vyy",this,other.vyy),
   vxy("vxy",this,other.vxy),
   mx("mx",this,other.mx),
   sx("sx",this,other.sx),
   ax("ax",this,other.ax),
   my("my",this,other.my),
   sy("sy",this,other.sy),
   ay("ay",this,other.ay),
   mz("mz",this,other.mz),
   sz("sz",this,other.sz),
   k("k",this,other.k),
   rho("rho",this,other.rho)
 { 
 } 



 Double_t BeamSpotPdf::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   double dx = x - (mx + ax*(z-mz));
   double dy = y - (my + ay*(z-mz));
   double dz = z - mz;
   double k2 = k*k;
   double wxx = sx*sx + k2*vxx;
   double wyy = sy*sy + k2*vyy;
   double wxy = rho*sx*sy + k2*vxy;
   double detw = wxx*wyy - wxy*wxy;
   double iwxx = wyy/detw;
   double iwyy = wxx/detw;
   double iwxy = -wxy/detw;
   double arg = dx*iwxx*dx + 2*dx*iwxy*dy + dy*iwyy*dy;
   arg += dz*dz/(sz*sz);
   return exp(-0.5*arg);
 }

Int_t BeamSpotPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
 { 
   //Previously, we checked to see if the integration range was larger than 5 sigma from the mean, and used 
   //the integral from -inf to +inf to approximate it. Uncomment the next block of code to implement that procedure:
   
   
 //  int mult = 5;
// 
  // if (matchArgs(allVars,analVars,x,y,z)){
    // if ((x.min(rangeName) < mx - mult*sx) && (x.max(rangeName) > mx + mult*sx)){
      // if ((y.min(rangeName) < my - mult*sy) && (y.max(rangeName) > my + mult*sy)){
	// if ((z.min(rangeName) < mz - mult*sz) && (z.max(rangeName) > mz + mult*sz)){
//	   return 1 ;
//	 }
//       }
//     }
//     
//   }
//   if (matchArgs(allVars,analVars,x,y)){
//     if ((x.min(rangeName) < mx - mult*sx ) && (x.max(rangeName) > mx + mult*sx)){
//       if ((y.min(rangeName) < my - mult*sy) && (y.max(rangeName) > my + mult*sy)){
//	 return 2 ; 
//       }
//     }
//     }	 
   

   //Since we don't care about the range of integration, we aren't using this rangeName variable.
   //The following line does nothing but silence the compiler warning.
   //(void) rangeName;

   //Now we can use the Erf approximation over the entire integration range:
     if (matchArgs(allVars,analVars,x,y,z)) return 1;
     if (matchArgs(allVars,analVars,x,y) ) return 2;
     if (matchArgs(allVars,analVars,x,z) ) return 3;
     if (matchArgs(allVars,analVars,y,z) ) return 4;
   
   return 0 ; 
 } 

 Double_t BeamSpotPdf::analyticalIntegral(Int_t code, const char* rangeName) const  
 { 
   // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
   // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
   // BOUNDARIES FOR EACH OBSERVABLE x

   double pi = TMath::Pi();
   double k2 = k*k;
   double wxx = sx*sx + k2*vxx;
   double wyy = sy*sy + k2*vyy;
   double wxy = rho*sx*sy + k2*vxy;
   double wzz = sz*sz;
   double detw = wxx*wyy - wxy*wxy;
  
   if (code==1) {
     //integrating over x,y,z
    
  //   return sqrt((0.5*pi)*(0.5*pi)*(0.5*pi))*sqrt(wxx)*sqrt(wyy)*sz*
  //     (TMath::Erf((mx-x.min(rangeName))/sqrt(2*wxx))-TMath::Erf((mx-x.max(rangeName))/sqrt(2*wxx)))
  //     *(TMath::Erf((my-y.min(rangeName))/sqrt(2*wyy))-TMath::Erf((my-y.max(rangeName))/sqrt(2*wyy)))
  //     *(TMath::Erf((mz-z.min(rangeName))/sqrt(2*wzz))-TMath::Erf((mz-z.max(rangeName))/sqrt(2*wzz)));
  //
     //The next line is used if we want to approximate the integral with that over -inf to +inf     
     return 2*pi*sqrt(detw)*sqrt(2*pi)*sz;
   } 
   
   
   if (code==2) {
     //integrating over x,y only
     double dz = z - mz;
     //erf approximation
    //
    //   return 0.5*pi*sqrt(wxx)*sqrt(wyy)*
    //   (TMath::Erf((mx-x.min(rangeName))/sqrt(2*wxx))-TMath::Erf((mx-x.max(rangeName))/sqrt(2*wxx)))
    //   *(TMath::Erf((my-y.min(rangeName))/sqrt(2*wyy))-TMath::Erf((my-y.max(rangeName))/sqrt(2*wyy)))
    //   *exp(-0.5*dz*dz/(sz*sz));
    //
       
     //inf approximation
     return 2*pi*sqrt(detw)*exp(-0.5*dz*dz/(sz*sz));
   }
   if (code==3){
     //integrating over x,z
     double dy = y-my;
     return 0.5*pi*sqrt(wxx)*sqrt(wzz)*
       (TMath::Erf((mx-x.min(rangeName))/sqrt(2*wxx))-TMath::Erf((mx-x.max(rangeName))/sqrt(2*wxx)))
       *(TMath::Erf((mz-z.min(rangeName))/sqrt(2*wzz))-TMath::Erf((mz-z.max(rangeName))/sqrt(2*wzz)))
       *exp(-0.5*dy*dy/(wyy));;
   }
   if (code==4){
     //integrating over y,z
     double dx = x-mx;
     return 0.5*pi*sqrt(wyy)*sqrt(wzz)*
       (TMath::Erf((my-y.min(rangeName))/sqrt(2*wyy))-TMath::Erf((my-y.max(rangeName))/sqrt(2*wyy)))
       *(TMath::Erf((mz-z.min(rangeName))/sqrt(2*wzz))-TMath::Erf((mz-z.max(rangeName))/sqrt(2*wzz)))
       *exp(-0.5*dx*dx/(wxx));;
   }
   
   return 0 ; 
 } 
