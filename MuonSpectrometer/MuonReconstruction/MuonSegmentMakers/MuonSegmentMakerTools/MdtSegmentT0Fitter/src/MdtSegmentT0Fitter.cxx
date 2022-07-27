/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MdtSegmentT0Fitter/MdtSegmentT0Fitter.h"

#include "MdtCalibSvc/MdtCalibrationSvcSettings.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"

#include "MdtCalibData/IRtRelation.h"
#include "MdtCalibData/IRtResolution.h"
#include "MdtCalibData/MdtRtRelation.h"

#include "MuonRIO_OnTrack/MdtDriftCircleOnTrack.h"
#include "MuonPrepRawData/MdtPrepData.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <functional>

#include <iostream>
#include <fstream>
#include <atomic>
#include <mutex>

namespace {
  // number of fit parameters
  constexpr unsigned int NUMPAR=3;
 
  // prints a message if a radius is bigger than this
  constexpr double MAX_RAD=16.;

  // time corresponding to r=15 mm for internal rt
  //constexpr double TUBE_TIME = 757.22;

  //constexpr double MAX_DRIFT= 855;

   // garbage time value to return when radius isn't wihin rt range
   constexpr double R2TSPURIOUS = 50000;
   
   constexpr int WEAK_TOPO_T0ERROR = 10;
   
   constexpr int STRONG_TOPO_T0ERROR = 50;

   struct HitCoords {
        public:
            HitCoords(const double z_coord, const double t_coord, 
                      const double y_coord, const double w_coord, 
                      const double r_coord, const MuonCalib::IRtRelation * rt_rel):
            z(z_coord),
            t(t_coord),
            y(y_coord),
            w(w_coord),
            r(r_coord),
            rt(rt_rel){}
        double z;
        double t;
        double y;
        double w;
        double r;
        const MuonCalib::IRtRelation *rt;
  };
  
  class FunctionToMinimize : public ROOT::Math::IMultiGenFunction {
    public:      
      FunctionToMinimize(const int used) : m_data(),m_used(used),m_t0Error(-1) {}
     
      double DoEval(const double* xx) const override {
        const double ang = xx[0];
        const double b = xx[1];
        const double t0 = xx[2];
        
        const double cosin = std::cos(ang);
        const double sinus = std::sin(ang);
        
        double fval = 0.;
        // Add t0 constraint
        if (m_t0Error == WEAK_TOPO_T0ERROR ) {
         fval += xx[2]*xx[2]/(1.0 *m_t0Error*m_t0Error);
        }        
        for(int i=0;i<m_used;++i) {
          const double t = m_data[i].t - t0;
          const double z = m_data[i].z;
          const double y = m_data[i].y;
          const double w = m_data[i].w;
          const double dist = std::abs(b*cosin + z*sinus - y*cosin); // same thing as fabs(a*z - y + b)/sqrt(1. + a*a);
          const double uppercut = m_data[i].rt->tUpper();
          const double lowercut = m_data[i].rt->tLower();
                   
          // Penalty for t<lowercut and t >uppercut
          if (t> uppercut ) { // too large
            fval += (t-uppercut)* (t-uppercut)*0.1;
          } else if (t < lowercut) {// too small
            fval += (t-lowercut)*(t-lowercut)*0.1;
          }
          const double r = t< lowercut ?  m_data[i].rt->radius(lowercut) : t > uppercut ? m_data[i].rt->radius(uppercut) :  m_data[i].rt->radius(t);
          fval += (dist - r)*(dist - r)*w;
        }
        
        return fval;
      }
      ROOT::Math::IBaseFunctionMultiDim* Clone() const override {return new FunctionToMinimize(m_used);}
      unsigned int NDim() const override {return 3;}
      void setT0Error(const int t0Error){m_t0Error=t0Error;}
      void addCoords(const double z, const double t, const double y, const double w, const double r, const MuonCalib::IRtRelation *rt){
        m_data.emplace_back(z,t,y,w,r,rt);
      }
    private:
      std::vector<HitCoords> m_data;
      int m_used;
      int m_t0Error;
  };
  
   /***********************************************************************************/
  /// RT function from Craig Blocker
  /// ok for now, possibly replace with actual RT function used to calibrate run

  //constexpr double T2R_A[] = {1.184169e-1, 3.32382e-2, 4.179808e-4, -5.012896e-6, 2.61497e-8, -7.800677e-11, 1.407393e-13, -1.516193e-16, 8.967997e-20, -2.238627e-23};
  //constexpr double RCORR_A[] = {234.3413, -5.803375, 5.061677e-2, -1.994959e-4, 4.017433e-7, -3.975037e-10, 1.522393e-13};

  

  double r2t_ext(std::vector<const MuonCalib::IRtRelation*> *rtpointers, double r, int i) {
    const MuonCalib::IRtRelation* rtrel = rtpointers->at(i);
    double ta = rtrel->tLower();
    double tb = rtrel->tUpper();
    if(r<rtrel->radius(ta) ) {
      return -1*R2TSPURIOUS;
    } else if(r>rtrel->radius(tb)) {
      return R2TSPURIOUS;
    }

    int itr = 0;
    while (ta <= tb) {
      double tm  = (ta + tb) / 2;  // compute mid point.
      double rtm = rtrel->radius(tm);
      if(std::abs(rtm - r) < 0.001 ) {
        return tm;
      }
      else if (r > rtm) {
        ta = tm;  // repeat search in top half.
      }
      else if (r < rtm ) {
        tb = tm; // repeat search in bottom half.
      }

      itr++;
      if(itr>50) return -1;
    }
    return -1;    // failed to find key
  }
  int sign(double a) {
    return a > 0 ? 1 : a < 0 ? -1 : 0;
  }
}

namespace TrkDriftCircleMath {

  MdtSegmentT0Fitter::MdtSegmentT0Fitter(const std::string& ty,const std::string& na,const IInterface* pa)
  : AthAlgTool(ty,na,pa),
    DCSLFitter(),
    m_ntotalCalls(0),
    m_npassedNHits(0),
    m_npassedSelectionConsistency(0),
    m_npassedNSelectedHits(0),
    m_npassedMinHits(0),
    m_npassedMinuitFit(0) {
    declareInterface <IDCSLFitProvider> (this);
  }

  StatusCode MdtSegmentT0Fitter::initialize() {
    ATH_CHECK(m_calibrationDbTool.retrieve());
    return StatusCode::SUCCESS;
  }

  StatusCode MdtSegmentT0Fitter::finalize() {

    double scaleFactor = m_ntotalCalls != 0 ? 1./(double)m_ntotalCalls : 1.;

    ATH_MSG_INFO( "Summarizing fitter statistics " << "\n"
                  << " Total fits       " << std::setw(10) << m_ntotalCalls << "   " << scaleFactor*m_ntotalCalls << "\n"
                  << " hits > 2         " << std::setw(10) << m_npassedNHits << "   " << scaleFactor*m_npassedNHits << "\n"
                  << " hit consis.      " << std::setw(10) << m_npassedSelectionConsistency << "   " << scaleFactor*m_npassedSelectionConsistency << "\n"
                  << " sel. hits > 2    " << std::setw(10) << m_npassedNSelectedHits << "   " << scaleFactor*m_npassedNSelectedHits << "\n"
                  << " Hits > min hits  " << std::setw(10) << m_npassedMinHits << "   " << scaleFactor*m_npassedMinHits << "\n"
                  << " Passed Fit       " << std::setw(10) << m_npassedMinuitFit << "   " << scaleFactor*m_npassedMinuitFit  );
    return StatusCode::SUCCESS;
  }

  bool MdtSegmentT0Fitter::fit( Segment& result, const Line& line, const DCOnTrackVec& dcs, const HitSelection& selection, double t0Seed ) const {
    ++m_ntotalCalls;

    ATH_MSG_DEBUG("New seg: ");

    const DCOnTrackVec& dcs_keep = dcs;

    unsigned int N = dcs_keep.size();

    result.setT0Shift(-99999,-99999);

    if(N<2) {
      return false;
    }
    ++m_npassedNHits;
    if( selection.size() != N ) {
      ATH_MSG_ERROR("MdtSegmentT0Fitter.cxx:fit with t0 <bad HitSelection>");
      return false;
    }
    ++m_npassedSelectionConsistency;
    int used=0;
    for(unsigned int i=0;i<N;++i){
      if( selection[i] == 0 ) ++used;
    }
    if(used < 2){
      ATH_MSG_DEBUG("TOO FEW HITS SELECTED");
      return false;
    }
    ++m_npassedNSelectedHits;
    
    if(used < m_minHits) {
      ATH_MSG_DEBUG("FEWER THAN Minimum HITS N " << m_minHits << " total hits " <<N<<" used " << used);

      //
      //     Copy driftcircles and reset the drift radii as they might have been overwritten
      //     after a succesfull t0 fit
      //

      DCOnTrackVec dcs_new;
      dcs_new.reserve(dcs.size());
      
      double chi2p = 0.;    
      int n_elements = dcs.size();
      for(int i=0; i< n_elements; ++i ){
          const DriftCircle* ds  = & dcs[i];
          if(std::abs(ds->r()-ds->rot()->driftRadius())>m_dRTol) ATH_MSG_DEBUG("Different radii on dc " << ds->r() << " rot " << ds->rot()->driftRadius());
          
          DriftCircle dc_keep(ds->position(), ds->rot()->driftRadius(), ds->dr(), ds->drPrecise(), ds->driftState(), ds->id(), ds->index(),ds->rot() );
          DCOnTrack dc_new(dc_keep, 0., 0.);
          
          dc_new.state(dcs[i].state());
          dcs_new.push_back( dc_new );
          if( selection[i] == 0 ){
            double t = ds->rot()->driftTime();
            const MuonCalib::MdtRtRelation *rtInfo = m_calibrationDbTool->getRtCalibration(ds->rot()->identify());
            
            double tUp = rtInfo->rt()->tUpper();
            double tLow = rtInfo->rt()->tLower();
            
            if(t<tLow) chi2p += (t-tLow)*(t-tLow)*0.1;
            else if(t>tUp) chi2p += (t-tUp)*(t-tUp)*0.1;
            }
      }
      
      if(chi2p>0) ATH_MSG_DEBUG("NO Minuit Fit TOO few hits Chi2 penalty " << chi2p);
      
      bool oldrefit = DCSLFitter::fit( result, line, dcs_new, selection );
      
      chi2p += result.chi2();
      // add chi2 penalty for too large or too small driftTimes  t < 0 or t> t upper
      result.set(chi2p, result.ndof(),  result.dtheta(),  result.dy0());
      int iok = 0;
      if(oldrefit) iok = 1;
      ATH_MSG_DEBUG(" chi2 total " << result.chi2() << " angle " << result.line().phi() << " y0 " << result.line().y0()  << " nhits "<< selection.size() << " refit ok " << iok);
      return oldrefit;
    
    }
    
    ATH_MSG_DEBUG("FITTING FOR T0 N "<<N<<" used " << used);
    

    ++m_npassedMinHits;

   
    ATH_MSG_DEBUG(" in  MdtSegmentT0Fitter::fit with N dcs "<< N << " hit selection size " <<  selection.size());
    ATH_MSG_DEBUG("in fit "<<result.hasT0Shift()<< " " <<result.t0Shift());
    

    double Zc(0);
    double Yc(0);
    double S(0);
    double Sz(0);
    double Sy(0);
    std::vector<double> y(N);
    std::vector<double> z(N);
    std::vector<double> w(N);
    std::vector<double> r(N);
    std::vector<double> dr(N);
    std::vector<double> t(N);
    std::vector<const MuonCalib::IRtRelation*> rtpointers(N);

    FunctionToMinimize minFunct(used);

    {
      DCOnTrackVec::const_iterator it = dcs_keep.begin();
      DCOnTrackVec::const_iterator it_end = dcs_keep.end();
      for(int ii=0 ;it!=it_end; ++it, ++ii ){
        const Muon::MdtDriftCircleOnTrack *roto = it->rot();
        if (!roto) {
          ATH_MSG_ERROR("MdtSegmentT0Fitter: NO DC ROT pointer found");
          return false;
        }
        ATH_MSG_DEBUG("hit # "<<ii );
        y[ii] = it->y();
        z[ii] = it->x();
        r[ii] = std::abs(roto->driftRadius());
        dr[ii] = it->dr();
        const Muon::MdtPrepData *peerd;
        peerd = dynamic_cast<const Muon::MdtPrepData*>(roto->prepRawData());
        if(!peerd) {
          ATH_MSG_ERROR("MdtSegmentT0Fitter: Can't convert to MdtPrepData* !! Not fitting for t0");
          return false;
        }
        Identifier id = roto->identify();
        const MuonCalib::MdtRtRelation *rtInfo = m_calibrationDbTool->getRtCalibration(id);
        rtpointers[ii] = rtInfo->rt();
        t[ii] = roto->driftTime();

        double newerror = m_scaleErrors ? it->drPrecise() : it->dr();

        if( newerror > 0.) w[ii] = 1./(newerror);
        else w[ii] = 0.;
        w[ii]*=w[ii];
        if(r[ii]<0){
          r[ii] = 0.;
          ATH_MSG_DEBUG("MdtSegmentT0Fitter (using times) ERROR: <Negative r> " << r[ii]);
        }

        ATH_MSG_DEBUG("DC:  (" << y[ii] << "," << z[ii] << ") R = " << r[ii] << " W " << w[ii] <<" t " <<t[ii]<< " id: "<<it->id()<<" sel " << selection[ii]);

        if( selection[ii] == 0 ){
          S+=w[ii];
          Sz+= w[ii]*z[ii];
          Sy+= w[ii]*y[ii];
          if(r[ii] > MAX_RAD ) {
            ATH_MSG_DEBUG("Radius is too big");
          }
        }
      }
    }
    
    const double inv_S = 1. / S;
    Zc = Sz*inv_S;
    Yc = Sy*inv_S;

    ATH_MSG_DEBUG("Yc " << Yc << " Zc " << Zc);

    /// go to coordinates centered at the average of the hits
    for(unsigned int i=0;i<N;++i) {
      y[i]  -= Yc;
      z[i]  -= Zc;
    }

    int selcount(0);
    DCOnTrackVec::const_iterator it = dcs_keep.begin();
    DCOnTrackVec::const_iterator it_end = dcs_keep.end();

    // replicate for the case where the external rt is used...
    // each hit has an rt function with some range...we want to fit such that
    // tlower_i < ti - t0 < tupper_i
    double min_tlower=1e10;
    double max_tupper=-1e10;

    double t0seed=0; // the average t0 of the hit
    double st0 = 0; // the std deviation of the hit t0s
    double min_t0 = 1e10; // the smallest t0 seen
    double tee0, tl, th;


    for(int ii=0 ;it!=it_end; ++it, ++ii ){
      if( selection[ii] == 0 ) {
        double r2tval = r2t_ext(&rtpointers,  r[ii], ii) ;
        tl = rtpointers[ii]->tLower();
        th = rtpointers[ii]->tUpper();
        if(t[ii] - tl < min_tlower) min_tlower = t[ii] - tl;
        if(t[ii] - th > max_tupper) max_tupper = t[ii] - th;
        tee0 = t[ii] - r2tval;

        
        ATH_MSG_DEBUG(" z "<<z[ii]
             <<" y "<<y[ii]
             <<" r "<<r[ii]
             <<" t "<<t[ii]
             <<" t0 "<<tee0
             <<" tLower "<<tl
             <<" tUpper "<<th);
          
        
        t0seed += tee0;
        st0 += tee0*tee0;
        if(tee0 < min_t0 && std::abs(r2tval) < R2TSPURIOUS) min_t0 = tee0;

        minFunct.addCoords(z[ii], t[ii], y[ii], w[ii], r[ii], rtpointers[ii]);

        selcount++;
      }
    }
    t0seed /= selcount;
    st0 = st0/selcount - t0seed*t0seed;
    st0 = st0 > 0. ? std::sqrt(st0) : 0.;

    ATH_MSG_DEBUG(" t0seed "<<t0seed<<" sigma "<<st0<< " min_t0 "<<min_t0);

    // ************************* seed the parameters
    const double theta = line.phi();
    double cosin = std::cos(theta);
    double sinus = std::sin(theta);

    if ( sinus < 0.0 ) {
      sinus = -sinus;
      cosin = -cosin;
    } else if ( sinus == 0.0 && cosin < 0.0 ) {
      cosin = -cosin;
    }
    
    ATH_MSG_DEBUG("before fit theta "<<theta<<" sinus "<<sinus<< " cosin "<< cosin);

    double d = line.y0() + Zc*sinus-Yc*cosin;

    
    ATH_MSG_DEBUG(" line x y "<<line.position().x()<<" "<<line.position().y());
    ATH_MSG_DEBUG(" Zc Yc "<< Zc <<" "<<Yc);
    ATH_MSG_DEBUG(" line x0 y0 "<<line.x0()<<" "<<line.y0());
    ATH_MSG_DEBUG(" hit shift " << -Zc*sinus+Yc*cosin);
    
// Calculate signed radii

    int nml1p = 0;
    int nml2p = 0;
    int nml1n = 0;
    int nml2n = 0;
    double sdist;
    it = dcs_keep.begin();
    for(int ii=0 ;it!=it_end; ++it, ++ii ){
      if( selection[ii] != 0 ) continue;
      sdist = d*cosin + z[ii]*sinus - y[ii]*cosin; // same thing as |a*z - y + b|/sqrt(1. + a*a);
      if(it->id().ml()==0&&sdist > 0) nml1p++;
      if(it->id().ml()==0&&sdist < 0) nml1n++;
      if(it->id().ml()==1&&sdist > 0) nml2p++;
      if(it->id().ml()==1&&sdist < 0) nml2n++;
    }

// Define t0 constraint in Minuit
    int t0Error = STRONG_TOPO_T0ERROR;
    if (nml1p+nml2p < 2 || nml1n+nml2n < 2) t0Error = WEAK_TOPO_T0ERROR;

    minFunct.setT0Error(t0Error);

// Reject topologies where in one of the Multilayers no +- combination is present
    if((nml1p<1||nml1n<1)&&(nml2p<1||nml2n<1)&&m_rejectWeakTopologies) {
       ATH_MSG_DEBUG("Combination rejected for positive radii ML1 " <<  nml1p << " ML2 " <<  nml2p << " negative radii ML1 " << nml1n << " ML " << nml2n << " used hits " << used << " t0 Error " << t0Error);
      it = dcs.begin();
      it_end = dcs.end();
      double chi2p = 0.;
      DCOnTrackVec dcs_new;
      dcs_new.reserve(dcs.size());
      for(int i=0; it!=it_end; ++it, ++i ){
	const DriftCircle* ds  = & dcs[i];
        if(std::abs(ds->r()-ds->rot()->driftRadius())>m_dRTol) ATH_MSG_DEBUG("Different radii on dc " << ds->r() << " rot " << ds->rot()->driftRadius());
        DriftCircle dc_keep(ds->position(), ds->rot()->driftRadius(), ds->dr(), ds->drPrecise(), ds->driftState(), ds->id(), ds->index(),ds->rot() );
        DCOnTrack dc_new(dc_keep, 0., 0.);
        dc_new.state(dcs[i].state());
        dcs_new.push_back( dc_new );
        if( selection[i] == 0 ){
          double t = ds->rot()->driftTime();
          const MuonCalib::MdtRtRelation *rtInfo = m_calibrationDbTool->getRtCalibration(ds->rot()->identify());
          double tUp = rtInfo->rt()->tUpper();
          double tLow = rtInfo->rt()->tLower();
          if(t<tLow) chi2p += (t-tLow)*(t-tLow)*0.1;
          if(t>tUp) chi2p += (t-tUp)*(t-tUp)*0.1;
        }
      }
      if(chi2p>0) ATH_MSG_DEBUG(" Rejected weak topology Chi2 penalty " << chi2p);
      bool oldrefit = DCSLFitter::fit( result, line, dcs_new, selection );
      chi2p += result.chi2();
// add chi2 penalty for too large or too small driftTimes  t < 0 or t> t upper
      result.set( chi2p, result.ndof(),  result.dtheta(),  result.dy0() );
      return oldrefit;
    }  // end rejection of weak topologies

    ATH_MSG_DEBUG("positive radii ML1 " <<  nml1p << " ML2 " <<  nml2p << " negative radii ML1 " << nml1n << " ML " << nml2n << " used hits " << used << " t0 Error " << t0Error);

    constexpr Double_t step[3] = {0.01 , 0.01 , 0.1 };
    // starting point
    Double_t variable[3] = {theta,d,0};
    // if t0Seed value from outside use this
    if(t0Seed > -999.) variable[2] = t0Seed;

    ROOT::Minuit2::Minuit2Minimizer minimum("algoName");
    minimum.SetMaxFunctionCalls(1000000);
    minimum.SetTolerance(0.001);
    minimum.SetPrintLevel(-1);
    if(msgLvl(MSG::VERBOSE)) minimum.SetPrintLevel(1);

    
    minimum.SetFixedVariable(0,"a", variable[0]);
    minimum.SetFixedVariable(1,"b", variable[1]);
    minimum.SetVariable(2,"t0",variable[2], step[2]);
    
    minimum.SetFunction(minFunct);

    // do the minimization
    minimum.Minimize();

    const double *results = minimum.X();
    const double *errors = minimum.Errors();
    ATH_MSG_DEBUG("Minimum: f(" << results[0] << "+-" << errors[0] << "," << results[1]<< "+-" << errors[1]<< "," << results[2] << "+-" << errors[2]<< "): " << minimum.MinValue());

    ++m_npassedMinuitFit;

    // Get the fit values
    double aret=results[0];
    double aErr=errors[0];
    double dtheta = aErr;
    double tana = std::tan(aret); // tangent of angle
    double ang = aret;  // between zero and pi
    cosin = std::cos(ang);
    sinus = std::sin(ang);
    if ( sinus < 0.0 ) {
      sinus = -sinus;
      cosin = -cosin;
    } else if ( sinus == 0.0 && cosin < 0.0 ) {
      cosin = -cosin;
    }
    ang = std::atan2(sinus, cosin);
    double b=results[1];
    double bErr=errors[1];
    double t0=results[2];
    double t0Err=errors[2];
    double dy0 = cosin * bErr - b * sinus * aErr;

    double del_t;
    del_t = std::abs(rtpointers[0]->radius((t0+t0Err)) - rtpointers[0]->radius(t0)) ;

    
    ATH_MSG_DEBUG("____________FINAL VALUES________________" );
    ATH_MSG_DEBUG("Values: a "<<tana<<" d "<<b * cosin <<" t0 "<<t0);
    ATH_MSG_DEBUG("Errors: a "<<aErr<<" b "<<dy0 <<" t0 "<<t0Err);
    
    d = b * cosin;
    if(msg().level() <=MSG::DEBUG) {
      msg() << MSG::DEBUG <<"COVAR  ";
      for(int it1=0; it1<3; it1++) {
        for(int it2=0; it2<3; it2++) {
          msg() << MSG::DEBUG <<minimum.CovMatrix(it1,it2)<<" ";
        }
        msg() << MSG::DEBUG << endmsg;
      }
    }

    result.dcs().clear();
    result.clusters().clear();
    result.emptyTubes().clear();

     ATH_MSG_DEBUG("after fit theta "<<ang<<" sinus "<<sinus<< " cosin "<< cosin);

    double chi2 = 0;
    unsigned int nhits(0);
    double yl;

    // calculate predicted hit positions from track parameters
    it = dcs_keep.begin();
    it_end = dcs_keep.end();
    ATH_MSG_DEBUG("------NEW HITS------");

    for(int i=0; it!=it_end; ++it, ++i ){
      double rad, drad;

      double uppercut = rtpointers[i]->tUpper();
      double lowercut = rtpointers[i]->tLower();
      rad = rtpointers[i]->radius(t[i]-t0);
      if(t[i]-t0<lowercut) rad = rtpointers[i]->radius(lowercut);
      if(t[i]-t0>uppercut) rad = rtpointers[i]->radius(uppercut);
      if (w[i]==0) {
        ATH_MSG_WARNING("w[i]==0, continuing");
        continue;
      }
      drad = 1.0/std::sqrt(w[i]) ;

      yl = (y[i] -  tana*z[i] - b);
      
      ATH_MSG_DEBUG("i "<<i<<" ");
      

      double dth = -(sinus*y[i] + cosin*z[i])*dtheta;
      double residuals = std::abs(yl)/std::sqrt(1+tana*tana) - rad;
      
      ATH_MSG_DEBUG(" dth "<<dth<<" dy0 "<<dy0<<" del_t "<<del_t);
      

      double errorResiduals = std::hypot(dth, dy0, del_t);

      // derivatives of the residual 'R'
      double deriv[3];
      // del R/del theta
      double dd = z[i] * sinus + b *cosin - y[i] * cosin;
      deriv[0] = sign(dd) * (z[i] * cosin - b * sinus + y[i] * sinus);
      // del R / del b
      deriv[1] = sign(dd) * cosin ;
      // del R / del t0

      deriv[2] = -1* rtpointers[i]->driftvelocity(t[i]-t0);

      double covsq=0;
      for(int rr=0; rr<3; rr++) {
        for(int cc=0; cc<3; cc++) {
          covsq += deriv[rr]*minimum.CovMatrix(rr,cc)* deriv[cc];
        }
      }
      ATH_MSG_DEBUG(" covsquared " << covsq);
      if( covsq < 0. && msg().level() <=MSG::DEBUG){
        for(int rr=0; rr<3; rr++) {
            for(int cc=0; cc<3; cc++) {
                double dot = deriv[rr]*minimum.CovMatrix(rr,cc)* deriv[cc];
                ATH_MSG_DEBUG(" adding term " << dot << " dev1 " << deriv[rr] << " cov " << minimum.CovMatrix(rr,cc) << " dev2 " << deriv[cc]);
            }
        }
      }
      
      covsq = covsq > 0. ? std::sqrt(covsq) : 0.;
      const DriftCircle* ds  = & dcs_keep[i];
      if (m_propagateErrors) drad = dr[i];
      
      DriftCircle dc_newrad(dcs_keep[i].position(), rad, drad, ds->driftState(), dcs_keep[i].id(), dcs_keep[i].index(),ds->rot() );
      DCOnTrack dc_new(dc_newrad, residuals, covsq);
      dc_new.state(dcs_keep[i].state());

      ATH_MSG_DEBUG("T0 Segment hit res "<<residuals<<" eres "<<errorResiduals<<" covsq "<<covsq<<" ri " << r[i]<<" ro "<<rad<<" drad "<<drad << " sel "<<selection[i]<< " inv error " << w[i]);

      if( selection[i] == 0 ) {
        ++nhits;
        if (!m_propagateErrors) {
          chi2 += residuals*residuals*w[i];
        } else {
          chi2 += residuals*residuals/(drad*drad);
        }
        ATH_MSG_DEBUG("T0 Segment hit res "<<residuals<<" eres "<<errorResiduals<<" covsq "<<covsq<<" ri " << r[i]<<" radius after t0 "<<rad<<" radius error "<< drad <<  " original error " << dr[i]);
// Put chi2 penalty for drift times outside window
        if (t[i]-t0> uppercut ) { // too large
	  chi2  += (t[i]-t0-uppercut)* (t[i]-t0-uppercut)*0.1;
        }else if (t[i]-t0 < lowercut ) {// too small
	  chi2 += ((t[i]-t0-lowercut)*(t[i]-t0-lowercut))*0.1;
        }
      }
      result.dcs().push_back( dc_new );
    }

    double oldshift = result.t0Shift();
    ATH_MSG_DEBUG("end fit old "<<oldshift<< " new " <<t0);
    // Final Minuit Fit result
    if(nhits==NUMPAR) {
      nhits++;
      chi2 += 1.;
    }
    result.set( chi2, nhits-NUMPAR, dtheta, -1.*dy0 );
    result.line().set( LocVec2D( Zc - sinus*d, Yc + cosin*d ), ang );
    if(t0==0.) t0=0.00001;
    result.setT0Shift(t0,t0Err);
   
    ATH_MSG_DEBUG("Minuit Fit complete: Chi2 " << chi2 << " angle " << result.line().phi() << " nhits "<< nhits  << " t0result " << t0);
    ATH_MSG_DEBUG("Minuit Fit complete: Chi2 " << chi2 << " angle " << result.line().phi() << " nhits "<<nhits<<" numpar "<<NUMPAR << " per dof " << chi2/(nhits-NUMPAR));
    ATH_MSG_DEBUG("Fit complete: Chi2 " << chi2 <<" nhits "<<nhits<<" numpar "<<NUMPAR << " per dof " << chi2/(nhits-NUMPAR)<<(chi2/(nhits-NUMPAR) > 5 ? " NOT ":" ")<< "GOOD");
    ATH_MSG_DEBUG("chi2 "<<chi2<<" per dof "<<chi2/(nhits-NUMPAR));
    
    return true;
  }

}
