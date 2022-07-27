/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TRT_SegmentToTrackTool/TRT_SegmentToTrackTool.h"
#include "TrkPseudoMeasurementOnTrack/PseudoMeasurementOnTrack.h"
#include "InDetRIO_OnTrack/TRT_DriftCircleOnTrack.h"
#include "InDetIdentifier/TRT_ID.h"
#include "TrkSurfaces/Surface.h"

///Needed for the track refitter
//Extrapolator tool
#include "TrkExInterfaces/IExtrapolator.h"
//Scoring tool
#include "TrkToolInterfaces/ITrackScoringTool.h"
#include <cmath>

using Amg::Vector3D;
using CLHEP::mm;

namespace InDet {

  // Constructor with parameters:
  TRT_SegmentToTrackTool::TRT_SegmentToTrackTool(const std::string &type,
							     const std::string &name,
							     const IInterface *parent) :
    AthAlgTool(type,name,parent),
    m_fieldUnitConversion(1000.),
    m_extrapolator("Trk::Extrapolator/InDetExtrapolator"),
    m_scoringTool("Trk::TrackScoringTool/TrackScoringTool"),
    m_trtId(nullptr)
  {
    declareInterface<InDet::ITRT_SegmentToTrackTool>( this );
    m_doRefit            = false                                ;       //Do a final careful refit of tracks
    m_suppressHoleSearch = false                                ;       //Suppress hole search
    m_sharedFrac         = 0.3                                  ;       //Maximum fraction of shared hits !!!!!!!!!!!!!!!!!!!!!! offline 0.3!!!!!!!!!!!!!!!!!!!!!!!
    declareProperty("FinalRefit"                 ,m_doRefit           ); //Do a final careful refit of tracks
    declareProperty("SuppressHoleSearch"         ,m_suppressHoleSearch); //Suppress hole search during the track summary creation
    declareProperty("MaxSharedHitsFraction"      ,m_sharedFrac        ); //Maximum fraction of shared drift circles
    declareProperty("Extrapolator"               ,m_extrapolator      ); //Extrapolator tool
    declareProperty("ScoringTool"                ,m_scoringTool       ); //Track scoring tool
  }

  TRT_SegmentToTrackTool::~TRT_SegmentToTrackTool()
  = default;


  StatusCode TRT_SegmentToTrackTool::initialize() {

    ATH_CHECK( AthAlgTool::initialize() );

    ATH_MSG_DEBUG( "Initializing TRT_SegmentToTrackTool" );

    ATH_CHECK(m_extrapolator.retrieve() );

    ATH_CHECK(m_fitterTool.retrieve(DisableTool{ !m_doRefit }));
    ATH_CHECK(m_assoTool.retrieve( DisableTool{m_assoTool.name().empty()} ));
    ATH_CHECK(m_trackSummaryTool.retrieve( DisableTool{m_trackSummaryTool.name().empty()} ));

    // Get the scoring tool
    ATH_CHECK( m_scoringTool.retrieve() );

    ATH_CHECK( detStore()->retrieve(m_trtId, "TRT_ID") );
    ////////////////////////////////////////////////////////////////////////////////
    ATH_CHECK( m_fieldCacheCondObjInputKey.initialize());
    ////////////////////////////////////////////////////////////////////////////////

    // Get output print level
    //
    ATH_MSG_DEBUG( *this );

    return StatusCode::SUCCESS;

  }

  StatusCode InDet::TRT_SegmentToTrackTool::finalize(){
    StatusCode sc = AthAlgTool::finalize();
    return sc;
  }


  ///////////////////////////////////////////////////////////////////
  // Dumps relevant information into the MsgStream
  ///////////////////////////////////////////////////////////////////
  MsgStream&  InDet::TRT_SegmentToTrackTool::dump( MsgStream& out ) const {
    out<<std::endl;
    return dumpconditions(out);
  }

  ///////////////////////////////////////////////////////////////////
  // Dumps conditions information into the MsgStream
  ///////////////////////////////////////////////////////////////////

  MsgStream& InDet::TRT_SegmentToTrackTool::dumpconditions( MsgStream& out ) const
  {
    int n = 65-m_fitterTool.type().size();
    std::string s1; for(int i=0; i<n; ++i) s1.append(" "); s1.append("|");
    n     = 65-m_assoTool.type().size();
    std::string s2; for(int i=0; i<n; ++i) s2.append(" "); s2.append("|");
    n     = 65-m_scoringTool.type().size();
    std::string s5; for(int i=0; i<n; ++i) s5.append(" "); s5.append("|");


    out<<std::endl
       <<"|----------------------------------------------------------------------"
       <<"-------------------|"<<std::endl;
    out<<"| Tool for final track refitting    | "<<m_fitterTool.type()   <<s1<<std::endl;
    out<<"| Tool for track scoring            | "<<m_scoringTool.type()  <<s5<<std::endl;
    out<<"| Association service               | "<<m_assoTool.type()     <<s2<<std::endl;
    out<<"|----------------------------------------------------------------------"
       <<"-------------------|";
    return out;
  }


  ///////////////////////////////////////////////////////////////////
  // Dumps event information into the MsgStream
  ///////////////////////////////////////////////////////////////////
  MsgStream& InDet::TRT_SegmentToTrackTool::dumpevent( MsgStream& out ) {
    return out;
  }

  ///////////////////////////////////////////////////////////////////
  // Dumps relevant information into the ostream
  ///////////////////////////////////////////////////////////////////

  std::ostream& InDet::TRT_SegmentToTrackTool::dump( std::ostream& out ) const
  {
    return out;
  }


  Trk::Track* TRT_SegmentToTrackTool::segToTrack(const EventContext& ctx, const Trk::TrackSegment& tS) const {

    ATH_MSG_DEBUG("Transforming the TRT segment into a track...");

    //
    // Get the track segment fit quality. If not there drop segment
    //
    if (!tS.fitQuality()) {
      ATH_MSG_DEBUG("Segment has no fit quality ! Discard...");
      return nullptr;
    }
    const Trk::FitQuality* fq = tS.fitQuality()->clone();

    //
    // Get the track segment information about the initial track parameters
    //
    const AmgVector(5)& par = tS.localParameters();
    AmgSymMatrix(5) ep = AmgSymMatrix(5)(tS.localCovariance());
    const Trk::TrackParameters* segPar =
      tS.associatedSurface()
        .createUniqueTrackParameters(par[Trk::loc1],
                                     par[Trk::loc2],
                                     par[Trk::phi],
                                     par[Trk::theta],
                                     par[Trk::qOverP],
                                     std::move(ep))
        .release();

    if (segPar) {
      ATH_MSG_VERBOSE("Initial TRT Segment Parameters : " << (*segPar));
    } else {
      ATH_MSG_DEBUG("Could not get initial TRT segment parameters! ");
      // clean up
      delete fq;
      fq = nullptr;
      return nullptr;
    }

    // --- create new track state on surface vector
    auto ntsos = DataVector<const Trk::TrackStateOnSurface>();

    //
    // if no refit, make it a perigee
    //
    if (!m_doRefit) {

      const Trk::TrackStateOnSurface* par_tsos = nullptr;

      // --- create surface at perigee
      Amg::Vector3D perigeePosition(0., 0., 0.);
      Trk::PerigeeSurface perigeeSurface(perigeePosition);
      // --- turn parameters into perigee...
      std::unique_ptr<const Trk::TrackParameters> tmp =
        m_extrapolator->extrapolate(ctx, *segPar, perigeeSurface);
      std::unique_ptr<const Trk::Perigee> perParm = nullptr;
      //pass ownership if of the right type
      if (tmp && tmp->associatedSurface().type() == Trk::SurfaceType::Perigee) {
        perParm.reset(static_cast<const Trk::Perigee*>(tmp.release()));
      }
      if (perParm) {
        ATH_MSG_VERBOSE("Perigee version of Parameters : " << (*segPar));
      } else {
        ATH_MSG_DEBUG("Failed to build perigee parameters.Discard...");
        ntsos.clear();
        delete segPar;
        segPar = nullptr;
        delete fq;
        fq = nullptr;
        return nullptr;
      }

      // now create a perigee TSOS
      std::bitset<Trk::TrackStateOnSurface::NumberOfTrackStateOnSurfaceTypes>
        typePattern;
      typePattern.set(Trk::TrackStateOnSurface::Perigee);
      par_tsos = new Trk::TrackStateOnSurface(
        nullptr, std::move(perParm), nullptr, nullptr, typePattern);
      // push new TSOS into the list
      ntsos.push_back(par_tsos);
    }

    //
    // now loop over the TSOS and copy them in
    //

    // psuedo measurement information
    int npseudo = 0;
    double pseudotheta = 0;
    const Trk::MeasurementBase* pseudo = nullptr;
    // other measurement information
    const Trk::Surface *firstsurf = nullptr, *firstecsurf = nullptr,
                       *lastsurf = nullptr;
    const Trk::MeasurementBase* firstmeas = nullptr;
    // counters for barrel and endcap
    int nbarrel = 0, nendcap = 0;
    // some variables for endcaps
    std::vector<std::pair<double, double>> points;
    points.reserve(40);
    double oldphi = 0;

    // loop over the measurements in track segment (tS)
    for (int it = 0; it < int(tS.numberOfMeasurementBases()); it++) {

      // the track state on service we like to constuct ...
      const Trk::TrackStateOnSurface* seg_tsos = nullptr;

      // is this ROT a psuedo-measurement ?
      if (dynamic_cast<const Trk::PseudoMeasurementOnTrack*>(
            tS.measurement(it))) {
        // did we have a speudo-measurement before ?
        if (pseudo) {
          // get theta from pseudo measurements
          pseudotheta =
            std::atan2(tS.measurement(it)->associatedSurface().center().perp(),
                  tS.measurement(it)->associatedSurface().center().z());
        }
        // keep this pseudo measurement
        pseudo = tS.measurement(it);
        // update counter
        npseudo++;

        if (m_doRefit) {
          // refit means we can simply copy the state, otherwise we skip it
          seg_tsos =
            new Trk::TrackStateOnSurface(tS.measurement(it)->uniqueClone(), nullptr);
        }

      } else {
        //
        // normal measurement, not a pseudo measurement
        //
        // copy measurement
        seg_tsos =
          new Trk::TrackStateOnSurface(tS.measurement(it)->uniqueClone(), nullptr);

        //
        // --- following is for the hack below
        //

        // remember first measurement
        if (!firstmeas)
          firstmeas = tS.measurement(it);
        if (!firstsurf)
          firstsurf = &tS.measurement(it)->associatedSurface();
        // it is always the last one
        lastsurf = &tS.measurement(it)->associatedSurface();

        // this is a rubbish way to find out it is endcap
        if (std::abs(tS.measurement(it)
                   ->associatedSurface()
                   .transform()
                   .rotation()
                   .col(2)
                   .z()) < .5) {
          // increase counter and keep some information
          nendcap++;
          if (!firstecsurf)
            firstecsurf = &tS.measurement(it)->associatedSurface();

          double tmpphi =
            tS.measurement(it)->associatedSurface().center().phi();
          if (!points.empty() &&
              std::abs(tmpphi - oldphi) > M_PI) { // correct for boundary at +/- pi
            if (tmpphi < 0)
              tmpphi += 2 * M_PI;
            else
              tmpphi -= 2 * M_PI;
          }
          // remember oldphi
          oldphi = tmpphi;

          // copy the points
          points.emplace_back(
            tS.measurement(it)->associatedSurface().center().z(), tmpphi);

        } else
          nbarrel++;

        //
        // --- end of hack stuff
        //
      }

      // push new TSOS into the list
      if (seg_tsos)
        ntsos.push_back(seg_tsos);
    }

    // Construct the new track
    Trk::TrackInfo info;
    info.setPatternRecognitionInfo(Trk::TrackInfo::TRTStandalone);

    // create new track candidate
    if (!m_doRefit) {
      return new Trk::Track(info, std::move(ntsos), fq);
    } else {
      //
      // ----------------------------- this is a horrible hack to make the
      // segments fittable
      //

      // in case of only 1 pseudo measurement, use the theta from it.
      if (npseudo == 1)
        pseudotheta = pseudo->localParameters()[Trk::theta];

      // we need new perigee parameters
      double myqoverp = 0, myphi = 0, myd0 = 0, myz0 = 0, mytheta = pseudotheta;

      if (nendcap < 4) {
        //
        // --- are we in the barrel mostly
        //

        // momentum
        myqoverp = par[4] * std::sin(pseudotheta) / std::sin(par[3]);

        // --- create surface at perigee
        Amg::Vector3D perigeePosition(0., 0., 0.);
        Trk::PerigeeSurface perigeeSurface(perigeePosition);
        // -- get perigee
        std::unique_ptr<const Trk::TrackParameters> tmp =
          m_extrapolator->extrapolate(ctx, *segPar, perigeeSurface);
        std::unique_ptr<const Trk::Perigee> tempper = nullptr;
        if (tmp && tmp->associatedSurface().type() == Trk::SurfaceType::Perigee) {
           tempper.reset(static_cast<const Trk::Perigee*>(tmp.release()));
        }
        if (!tempper) {
          ATH_MSG_DEBUG("Could not produce perigee");
          delete segPar;
          segPar = nullptr;
          return nullptr;
        }

        // keep some values
        myd0 = tempper->parameters()[Trk::d0];
        myphi = tempper->parameters()[Trk::phi0];

      } else {
        //
        // --- endcap or transition track
        //

        // get estimate of parameters
        double sx = 0, sy = 0, sxx = 0, sxy = 0, d = 0;
        float zmin = 0, zmax = 0;
        // loop over all points
        for (unsigned int i = 0; i < points.size(); i++) {
          sx += points[i].first;
          sy += points[i].second;
          sxy += points[i].first * points[i].second;
          sxx += points[i].first * points[i].first;
          if (fabs(points[i].first) > fabs(zmax)) {
            zmax = points[i].first;
          }
          if (fabs(points[i].first) < fabs(zmin)) {
            zmin = points[i].first;
          }
        }

        if (std::abs(pseudotheta) < 1.e-6) {
          ATH_MSG_DEBUG("pseudomeasurements missing on the segment?");
          const float Rinn = 644., Rout = 1004.;
          if (zmax * zmin > 0.) {
            pseudotheta = std::atan2(Rout - Rinn, zmax - zmin);
          } else if (std::abs(zmax * zmin) < 1.e-6) {
            if (std::abs(zmax) > 1.e-6) {
              pseudotheta = std::atan2(Rout, zmax);
            } else {
              ATH_MSG_DEBUG("no points in endcap?");
            }
          } else {
            pseudotheta = std::atan2(2. * Rout, zmax - zmin);
          }
        }

        // get q/p
        d = (points.size() * sxx - sx * sx);
        double dphidz = ((points.size() * sxy - sy * sx) / d);
        myqoverp = (std::abs(pseudotheta) > 1e-6)
                     ? -dphidz / (0.6 * std::tan(pseudotheta))
                     : 1000.;

        // some geometry stuff to estimate further paramters...
        double halfz = 200.;
        const Trk::CylinderBounds* cylb =
          dynamic_cast<const Trk::CylinderBounds*>(&firstsurf->bounds());
        if (!cylb)
          ATH_MSG_ERROR("Cast of bounds failed, should never happen");
        else
          halfz = cylb->halflengthZ();
        const Trk::CylinderBounds* cylb2 =
          dynamic_cast<const Trk::CylinderBounds*>(&lastsurf->bounds());
        double halfz2 = 200.;
        if (!cylb2)
          ATH_MSG_ERROR("Cast of bounds failed, should never happen");
        else
          halfz2 = cylb2->halflengthZ();
        Amg::Vector3D strawdir1 = -firstsurf->transform().rotation().col(2);
        Amg::Vector3D strawdir2 = -lastsurf->transform().rotation().col(2);
        Amg::Vector3D pos1;
        Amg::Vector3D pos2;

        // ME: this is hardcoding, not nice and should be fixed
        if (std::abs(lastsurf->center().z()) < 2650 * mm) {
          pos2 = lastsurf->center() + halfz2 * strawdir2;
          if (nbarrel == 0) {
            double dr = std::abs(std::tan(pseudotheta) * (lastsurf->center().z() -
                                                 firstsurf->center().z()));
            pos1 = firstsurf->center() + (halfz - dr) * strawdir1;
          } else {
            double dz = std::abs((pos2.perp() - firstsurf->center().perp()) /
                             std::tan(pseudotheta));
            if (pos2.z() > 0)
              dz = -dz;
            double z1 = pos2.z() + dz;
            pos1 = Amg::Vector3D(
              firstsurf->center().x(), firstsurf->center().y(), z1);
          }
        } else {
          double dr = std::abs(std::tan(pseudotheta) *
                           (lastsurf->center().z() - firstsurf->center().z()));
          pos2 = lastsurf->center() + (dr - halfz2) * strawdir2;
          pos1 = firstsurf->center() - halfz * strawdir1;
        }

        // ME: I don't understand this yet, why is this done only if barrel ==
        // 0, while above this nendcap < 4 ?
        if (nbarrel == 0 &&
            std::abs(std::tan(pseudotheta) * (firstsurf->center().z() -
                                         lastsurf->center().z())) < 250 * mm &&
            std::abs(firstsurf->center().z()) > 1000 * mm) {

          // ME: wow this is hacking the vector ...
          const Trk::MeasurementBase* firstmeas =
            (**ntsos.begin()).measurementOnTrack();
          Amg::MatrixX newcov(2, 2);
          newcov.setZero();
          newcov(0, 0) = (firstmeas->localCovariance())(0, 0);
          newcov(1, 1) = (myqoverp != 0) ? .0001 * myqoverp * myqoverp : 1.e-8;
          Trk::LocalParameters newpar(std::make_pair(0, Trk::locZ),
                                      std::make_pair(myqoverp, Trk::qOverP));
          auto newpseudo =
            std::make_unique<Trk::PseudoMeasurementOnTrack>(
              newpar, newcov, firstmeas->associatedSurface());
          // hack replace first measurement with pseudomeasurement
          ntsos.erase(ntsos.begin());
          ntsos.insert(ntsos.begin(),
                        new Trk::TrackStateOnSurface(std::move(newpseudo), nullptr));
        }

        Amg::Vector3D field1;

        MagField::AtlasFieldCache fieldCache;

        // Get field cache object
        SG::ReadCondHandle<AtlasFieldCacheCondObj> readHandle{
          m_fieldCacheCondObjInputKey, ctx
        };
        const AtlasFieldCacheCondObj* fieldCondObj{ *readHandle };
        if (fieldCondObj == nullptr) {
          ATH_MSG_ERROR(
            "segToTrack: Failed to retrieve AtlasFieldCacheCondObj with key "
            << m_fieldCacheCondObjInputKey.key());
          return nullptr;
        }
        fieldCondObj->getInitializedCache(fieldCache);

        //   MT version uses cache
        fieldCache.getField(Amg::Vector3D(.5 * (pos1 + pos2)).data(),
                            field1.data());

        field1 *= m_fieldUnitConversion; // field in Tesla

        double phideflection =
          -.3 * (pos2 - pos1).perp() * field1.z() * myqoverp / std::sin(pseudotheta);
        double precisephi = (nbarrel == 0)
                              ? (pos2 - pos1).phi() - .5 * phideflection
                              : (pos2 - pos1).phi() + .5 * phideflection;
        double radius = (myqoverp != 0. && field1.z() != 0.)
                          ? -std::sin(pseudotheta) / (.3 * field1.z() * myqoverp)
                          : 1.e6;
        double precisetheta =
          (myqoverp != 0.)
            ? std::atan2(std::abs(radius * phideflection), pos2.z() - pos1.z())
            : pseudotheta;
        if (precisetheta < 0)
          precisetheta += M_PI;
        if (precisephi < -M_PI)
          precisephi += 2 * M_PI;
        if (precisephi > M_PI)
          precisephi -= 2 * M_PI;

        // now extrapolate backwards from the first surface to get starting
        // parameters
        const Trk::StraightLineSurface* surfforpar = nullptr;
        if (nbarrel == 0)
          surfforpar = dynamic_cast<const Trk::StraightLineSurface*>(firstsurf);
        else
          surfforpar = dynamic_cast<const Trk::StraightLineSurface*>(lastsurf);
        if (!surfforpar)
          ATH_MSG_ERROR("Cast of surface failed, should never happen");

        Trk::AtaStraightLine ataline(((nbarrel == 0) ? pos1 : pos2),
                                     precisephi,
                                     precisetheta,
                                     myqoverp,
                                     *surfforpar);
        Trk::PerigeeSurface persurf;
        const Trk::TrackParameters* extrappar =
          m_extrapolator->extrapolateDirectly(ctx, ataline, persurf).release();

        // now get parameters
        if (extrappar) {
          if (nendcap >= 4) {
            myphi = extrappar->parameters()[Trk::phi0];
            myd0 = extrappar->parameters()[Trk::d0];
          }

          // construct theta again
          double z0 = extrappar->parameters()[Trk::z0];
          if (nbarrel == 0)
            mytheta = std::atan(std::tan(extrappar->parameters()[Trk::theta]) *
                           std::abs((z0 - pos1.z()) / pos1.z()));
          else
            mytheta = std::atan(std::tan(extrappar->parameters()[Trk::theta]) *
                           std::abs((z0 - pos2.z()) / pos2.z()));

          if (mytheta < 0)
            mytheta += M_PI;

          delete extrappar;
          extrappar = nullptr;
        }
      }
      while (myphi > M_PI)
        myphi -= 2 * M_PI;
      while (myphi < -M_PI)
        myphi += 2 * M_PI;

      double P[5] = { myd0, myz0, myphi, mytheta, myqoverp };

      // create perigee TSOS and add as first (!) TSOS
      
      auto per =
        std::make_unique<Trk::Perigee>(P[0], P[1], P[2], P[3], P[4], Trk::PerigeeSurface());
      std::bitset<Trk::TrackStateOnSurface::NumberOfTrackStateOnSurfaceTypes>
        typePattern;
      typePattern.set(Trk::TrackStateOnSurface::Perigee);
      Trk::TrackStateOnSurface* seg_tsos = new Trk::TrackStateOnSurface(
        nullptr, std::move(per), nullptr, nullptr, typePattern);
      ntsos.insert(ntsos.begin(), seg_tsos);

      ATH_MSG_VERBOSE("Constructed perigee at input to fit : " << (*per));

      //
      // ------------------------------------------------------- now refit the
      // track
      //

      Trk::Track newTrack (info, std::move(ntsos), fq);
      Trk::Track* fitTrack =
        m_fitterTool->fit(ctx,newTrack, true, Trk::nonInteracting).release();

      // cleanup
      if (segPar) {
        delete segPar;
        segPar = nullptr;
      }

      if (!fitTrack) {
        ATH_MSG_DEBUG("Refit of TRT track segment failed!");
        return nullptr;
      }

      //
      // -------------------------------------- hack the covarinace back to
      // something reasonable
      //
      const Trk::TrackParameters* firstmeaspar = nullptr;
      DataVector<const Trk::TrackParameters>::const_iterator parit =
        fitTrack->trackParameters()->begin();
      do {
        // skip pesudo measurements on perigee
        if ((*parit)->covariance() &&
            ((*parit)->associatedSurface() == tS.associatedSurface()))
          firstmeaspar = *parit;
        ++parit;
      } while (firstmeaspar == nullptr &&
               parit != fitTrack->trackParameters()->end());

      // Create new perigee starting from the modified first measurement that
      // has a more reasonable covariance matrix
      // const Trk::Perigee* perTrack=dynamic_cast<const
      // Trk::Perigee*>(fitTrack->perigeeParameters());
      const Trk::Perigee* perTrack = fitTrack->perigeeParameters();

      if (!perTrack || !perTrack->covariance()) {
        ATH_MSG_ERROR("Cast of perigee fails, should never happen !");
        return nullptr;
      } else {
        ATH_MSG_VERBOSE ("Perigee after refit with fudges to make it converge : " << (*perTrack) );

	if(firstmeaspar && firstmeaspar->position().perp()<2000*mm && std::abs(firstmeaspar->position().z())<3000*mm){

	  // Modify first measurement so that it has reasonable errors on z and theta
	  AmgSymMatrix(5) fcovmat = AmgSymMatrix(5)(*(firstmeaspar->covariance()));
	  // factors by which we like to scale the cov, this takes the original segment errors into account
	  double scaleZ     = std::sqrt(tS.localCovariance()(1,1))/std::sqrt( (fcovmat)(1,1));
	  double scaleTheta = std::sqrt(tS.localCovariance()(3,3))/std::sqrt( (fcovmat)(3,3));
	  // now do it
	  fcovmat(1,0)=scaleZ*((fcovmat)(1,0));
	  fcovmat(0,1) = (fcovmat)(1,0);
	  fcovmat(1,1)=tS.localCovariance()(1,1);
	  fcovmat(2,1)=scaleZ*((fcovmat)(2,1));
	  fcovmat(1,2) = (fcovmat)(2,1);
	  fcovmat(3,1)=scaleZ*scaleTheta*((fcovmat)(3,1));
	  fcovmat(1,3) = (fcovmat)(3,1);
	  fcovmat(4,1)=scaleZ*((fcovmat)(4,1));
	  fcovmat(1,4) = (fcovmat)(4,1);
	  fcovmat(3,0)=scaleTheta*((fcovmat)(3,0));
	  fcovmat(0,3) = (fcovmat)(3,0);
	  fcovmat(3,2)=scaleTheta*((fcovmat)(3,2));
	  fcovmat(2,3) = (fcovmat)(3,2);
	  fcovmat(3,3)=tS.localCovariance()(3,3);
	  fcovmat(4,3)=scaleTheta*((fcovmat)(4,3));
	  fcovmat(3,4) = (fcovmat)(4,3);

	  // const Amg::VectorX& par = firstmeaspar->parameters();
	  const AmgVector(5)& par = firstmeaspar->parameters();
          const Trk::TrackParameters* updatedPars =
            firstmeaspar->associatedSurface().createUniqueTrackParameters(
              par[Trk::loc1],
              par[Trk::loc2],
              par[Trk::phi],
              par[Trk::theta],
              par[Trk::qOverP],
              std::move(fcovmat)).release();

          // now take parameters at first measurement and exptrapolate to perigee
          const Trk::TrackParameters* newperpar =
            m_extrapolator->extrapolate(ctx,
                                        *updatedPars,
                                        perTrack->associatedSurface(),
                                        Trk::anyDirection,
                                        false,
                                        Trk::nonInteracting).release();
          delete updatedPars; updatedPars = nullptr;

	  if (!newperpar || !newperpar->covariance()) {
	    ATH_MSG_WARNING ("Can not hack perigee parameters, extrapolation failed");
	    delete newperpar; newperpar = nullptr;
	  } else {
	    // this is a HACK !!!
            // perTrack is owned by fitTrack which is not const here
            // thus the const-ness is only removed from somthting which is not strictly const here.
	    AmgSymMatrix(5)& errmat ATLAS_THREAD_SAFE = const_cast<AmgSymMatrix(5)&>(*perTrack->covariance());
	    // overwrite cov in perTrack
	    errmat = *newperpar->covariance();
	    delete newperpar; newperpar = nullptr;
	    // check that new cov makes sense !
	    const AmgSymMatrix(5)& CM = *perTrack->covariance();
	    if( CM(1,1)==0.||CM(3,3)==0. ) {
	      ATH_MSG_DEBUG ("Hacked perigee covariance is CRAP, reject track");
	      delete fitTrack; return nullptr;
	    } else {
	      ATH_MSG_VERBOSE ("Perigee after fit with scaled covariance matrix : " << *perTrack);
	    }
	  }
	}
      }
      // return fitted track
      return fitTrack;
    }
  }


  bool TRT_SegmentToTrackTool::segIsUsed(const Trk::TrackSegment& tS,
                                         const Trk::PRDtoTrackMap *prd_to_track_map) const {

    ATH_MSG_DEBUG ("Checking whether the TRT segment has already been used...");

    // some counters to be handled
    int nShared = 0;  // Shared drift circles in segment
    int nHits   = 0;  // Number of TRT measurements

    if (!m_assoTool.name().empty() && !prd_to_track_map) ATH_MSG_ERROR("PRDtoTrackMap to be used but not provided by the client");
    // loop over the track states
    for(int it=0; it<int(tS.numberOfMeasurementBases()); ++it){

      // remove pseudo measurements
      if ( dynamic_cast<const Trk::PseudoMeasurementOnTrack*>(tS.measurement(it)) )
        continue;

      // get the measurment
      const InDet::TRT_DriftCircleOnTrack* trtcircle = dynamic_cast<const InDet::TRT_DriftCircleOnTrack*>(tS.measurement(it));
      if (!trtcircle) continue;

      // get PRD measurement
      const InDet::TRT_DriftCircle* RawDataClus=dynamic_cast<const InDet::TRT_DriftCircle*>(trtcircle->prepRawData());
      if(!RawDataClus) continue;

      // count up number of hits
      nHits++;

      if(!m_assoTool.name().empty() && prd_to_track_map && prd_to_track_map->isUsed(*RawDataClus)) nShared++;
    }

    if(nShared >= int(m_sharedFrac * nHits)) {
      ATH_MSG_DEBUG ("Too many shared hits.Will drop the TRT segment");
      return true;
    } else {
      return false;
    }


  }

  bool TRT_SegmentToTrackTool::toLower(const Trk::TrackSegment& tS) const {

    ATH_MSG_DEBUG ("Try to recover low TRT DC segments in crack...");

    // counters
    int nEC = 0; int nBRL = 0; int firstWheel = -999; int lastLayer = -999;

    // loop over the track states
    for(int it=0; it<int(tS.numberOfMeasurementBases()); ++it){

      //test if it is a pseudo measurement
      if ( dynamic_cast<const Trk::PseudoMeasurementOnTrack*>(tS.measurement(it)) )
	continue;

      // get measurement
      const InDet::TRT_DriftCircleOnTrack* trtcircle = dynamic_cast<const InDet::TRT_DriftCircleOnTrack*>(tS.measurement(it));
      if(!trtcircle) continue;

      // get identifier
      Identifier id = trtcircle->detectorElement()->identify();
      // barrel or endcap
      int isB = m_trtId->barrel_ec(id);
      if      (isB==2 || isB==-2) {
	nEC++;
	if(nEC == 1)
	  firstWheel = m_trtId->layer_or_wheel(id);
      }
      else if (isB==1 || isB==-1) {
	nBRL++;
	lastLayer  = m_trtId->layer_or_wheel(id);
      }

    }

    // now the logic
    return (nEC>0  && nBRL>0) ||

       (nEC==0 && nBRL>0  && lastLayer<2) ||

       (nEC>0  && nBRL==0 && (firstWheel>10 || firstWheel<2));

  }

  void TRT_SegmentToTrackTool::addNewTrack(Trk::Track* trk, ITRT_SegmentToTrackTool::EventData &event_data) const {
    // @TODO avoid non const member m_trackScoreTrackMap
    ATH_MSG_DEBUG ("Add track to the scoring multimap...");
    if (m_trackSummaryTool.isEnabled()) {
       m_trackSummaryTool->computeAndReplaceTrackSummary(*trk,
                                                         nullptr,
                                                         m_suppressHoleSearch);
    }

    //Score the track under investigation
    Trk::TrackScore score = m_scoringTool->score(*trk,m_suppressHoleSearch);
    ATH_MSG_DEBUG ("TRT-only: score is " << score);

    if (score==0) {
      // statistics...
      ATH_MSG_DEBUG ("Track score is zero, reject it");
      event_data.m_counter[ITRT_SegmentToTrackTool::EventData::knTrkScoreZero]++;
      // clean up
      delete trk;
    } else {
      // add track to map, map is sorted small to big !
      event_data.m_trackScores.emplace_back(-score, trk );
    }

    
  }


  TrackCollection* TRT_SegmentToTrackTool::resolveTracks(const Trk::PRDtoTrackMap *prd_to_track_map_in,
                                                          ITRT_SegmentToTrackTool::EventData &event_data) const {

    ATH_MSG_DEBUG ("Resolving the TRT tracks in score map...");

    //    if (m_assoTool.name().empty() && !prd_to_track_map_in) ATH_MSG_ERROR("PRDtoTrackMap to be used but not provided by the client");

    std::unique_ptr<Trk::PRDtoTrackMap> prd_to_track_map;
    if (!m_assoTool.name().empty()) {
      prd_to_track_map=m_assoTool->createPRDtoTrackMap();
      if (prd_to_track_map_in) {
        *prd_to_track_map = *prd_to_track_map_in;
      }
    }

    //final copy - ownership is passed out of algorithm
    std::unique_ptr<TrackCollection> final_tracks = std::make_unique<TrackCollection>();
    final_tracks->reserve( event_data.m_trackScores.size());
    std::stable_sort(event_data.m_trackScores.begin(),
                     event_data.m_trackScores.end(),
                     []( const std::pair< Trk::TrackScore, Trk::Track* > &a,
                         const std::pair< Trk::TrackScore, Trk::Track* > &b)
                     {  return a.first < b.first; });

    for (std::pair< Trk::TrackScore, Trk::Track* > &track_score : event_data.m_trackScores) {

      ATH_MSG_DEBUG ("--- Trying next track "<<track_score.second<<"\t with score "<<-track_score.first);

      // some counters to be handled
      int nShared = 0;  // Shared drift circles in segment
      int nHits   = 0;  // Number of TRT measurements

      // get vector of TSOS
      const DataVector<const Trk::TrackStateOnSurface>* tsos = (track_score.second)->trackStateOnSurfaces();

      // loop over vector of TSOS
      for ( const Trk::TrackStateOnSurface *a_tsos : *tsos) {

	// get measurment from TSOS
	const Trk::MeasurementBase* meas = a_tsos->measurementOnTrack();
	if (!meas) continue;

	// make sure it is a TRT_DC and not a pseudo measurement
	const InDet::TRT_DriftCircleOnTrack* rot = dynamic_cast <const InDet::TRT_DriftCircleOnTrack*> (meas);
	if( !rot ) continue;

	// get to the PRD object
	const InDet::TRT_DriftCircle* RawDataClus=dynamic_cast<const InDet::TRT_DriftCircle*>(rot->prepRawData());
	if(!RawDataClus) continue;

	// count up number of hits
	nHits++;

	// count up number of shared hits
	if(!m_assoTool.name().empty() && prd_to_track_map->isUsed(*RawDataClus)) nShared++;
      }

      ATH_MSG_DEBUG ("TRT-only has " << nHits << " hits and " << nShared << " of them are shared");

      // cut on the number of shared hits with the max fraction
      if(nShared >=  int(m_sharedFrac * nHits)) {
         // statistics
         event_data.m_counter[ITRT_SegmentToTrackTool::EventData::knTrkSegUsed]++;
         ATH_MSG_DEBUG ("Too many shared hits, remove it !");
         delete track_score.second;
         continue;
      }

      // ok, this seems like a useful track
      final_tracks->push_back(track_score.second);

      event_data.m_counter[ITRT_SegmentToTrackTool::EventData::knTRTTrk]++;
      ATH_MSG_DEBUG ("TRT-only is accepted");

      //Register the track with the association tool
      if(!m_assoTool.name().empty()) {
        if(m_assoTool->addPRDs(*prd_to_track_map,*(track_score.second)).isFailure()) {
	  ATH_MSG_WARNING ("addPRDs() failed!");
	}
      }

    }

    return final_tracks.release();

  }

}
