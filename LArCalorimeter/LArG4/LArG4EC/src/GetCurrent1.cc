/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include <cassert>

#include "GaudiKernel/MsgStream.h"

#include "G4ThreeVector.hh"
#include "globals.hh"

#include "LArG4Code/LArG4BirksLaw.h"

#include "LArG4EC/EnergyCalculator.h"

#include "CLHEP/Units/SystemOfUnits.h"

static inline int signof(double a){ return (a == 0.) ? 0: (a < 0 ? -1: 1); }

using namespace LArG4::EC;

G4double EnergyCalculator::GetCurrent1(const G4ThreeVector &P1, const G4ThreeVector &P2, G4double edep) // need to make const
{
/*
	FILE *F = fopen("test1.dat", "w");
	if(F == 0) abort();
	G4ThreeVector a(0., 1200., 0.);
	for(double z = 0; z <= WheelThickness; z += 0.1){
		fprintf(F, "%f", z);
//		for(double x = -20.; x <= 20.; x += 1.){
//			a[0] = x;
			a[2] = z;
//			std::cout << "(" << z << ", " << x << ") " << DistanceToTheNeutralFibre(a) << std::endl;
			for(int s = -1; s <= 1; s ++){
//				std::cout << "\t" << s << " " << elc()->AmplitudeOfSurface(a, s) << std::endl;
				fprintf(F, " %f", elc()->AmplitudeOfSurface(a, s));
			}
//		}
		fprintf(F, "\n");
	}
	fclose(F);

	abort();
*/

//std::cout << "GetCurrent1 -------------------------------------------------------------------------" << std::endl;
//std::cout << "input: (" << P1.x() << ", " << P1.y() << ", " << P1.z() << "), ("
//          << P2.x() << ", " << P2.y() << ", " << P2.z() << "), " << edep << std::endl;

	G4int gaperr = 0;

	const std::pair<G4int, G4int> gap2 = lwc()->GetPhiGapAndSide(P2);
	const std::pair<G4int, G4int> gap1 = lwc()->GetPhiGapAndSide(P1);
//	gap1.first = PhiGapNumberForWheel(gap1.first);
//	gap2.first = PhiGapNumberForWheel(gap2.first);

	if(gap1 != gap2) gaperr = -1;

	//JT.>>> 
	SetFoldArea(P1.z());                 // set fold type
	SetHalfWave(P1.z());                 // set halfwave parameters for substep
	G4double p1[3];
	p1[0] = P1.x();
	p1[1] = P1.y();
	p1[2] = P1.z();

	//std::cout<<"===>>>GetCurrent1: before GetPhiGap** CosPhiGap,SinPhiGap="
	// <<CosPhiGap<<" "<<SinPhiGap<<std::endl;

	GetPhiGap(p1);  //set rotation angle from the wheel to the efield Map system; once for each G4 step

	//std::cout<<"===>>>GetCurrent1: after GetPhiGap** CosPhiGap,SinPhiGap="
	// <<CosPhiGap<<" "<<SinPhiGap<<std::endl;
	//<<<JT

	const G4double step = (P2 - P1).mag();
	const G4int nofstep = G4int(step / GridSize) + 1;  // step is divided to substeps
	G4double current = 0.;                          //current to be returned
	const G4double step_current = edep / nofstep / AverageCurrent; // base current for each step

//std::cout << "gap1: " << gap1.first << ", " << gap1.second
//          << "; gap2: " << gap2.first << ", " << gap2.second << std::endl;
//std::cout << "input: (" << P1.x() << ", " << P1.y() << ", " << P1.z() << "), ("
//          << P2.x() << ", " << P2.y() << ", " << P2.z() << "), " << edep << std::endl;


	for(G4int i = 0; i < nofstep; ++ i){
		G4double ds = (i + 0.5) / nofstep;
		G4ThreeVector Pe = P1 * (1. - ds) + P2 * ds;
		G4ThreeVector Pa = Pe;
//std::cout << "step " << i << std::endl;
//std::cout << "\tpoint (" << Pe.x() << ", " << Pe.y() << ", " << Pe.z() << ")" << std::endl;

		SetFoldArea(Pe.z());                 // set fold type
		SetHalfWave(Pe.z());                 // set halfwave parameters for substep
		G4double vstep[3], vmap[3];
		vstep[0] = Pe.x();
		vstep[1] = Pe.y();
		vstep[2] = Pe.z();
		TransformWheeltoFieldMap(vstep, vmap);  //get corresponding point in Map

		//JT.>>> 
		const G4double rvstep2=vstep[0]*vstep[0]+vstep[1]*vstep[1] ;
		const G4double rforalpha=sqrt( rvstep2-vmap[1]*vmap[1] );
		const G4double gap=HalfLArGapSize(rforalpha,rforalpha);// LAr gap at the straight section on the fieldmap
		const G4double yshift_on_map  =    rforalpha*M_PI/lwc()->GetNumberOfFans()-(FanAbsThickness+FanEleThickness)/2.;
		const G4double yshift_on_wheel=sqrt(rvstep2)*M_PI/lwc()->GetNumberOfFans()-(FanAbsThickness+FanEleThickness)/2.;
		const G4double cylgapcorr=yshift_on_wheel/yshift_on_map; // scale difference between plane and cylindrical surface
		/*
		std::cout<< " GetCurrent1**Nabs="<<lwc()->GetNumberOfFans()<<" absthick="<<FanAbsThickness<<" elethick="<<FanEleThickness
			 <<" cylgapcorr-1="<<cylgapcorr-1
			 <<" ZinHalfWave="<<ZinHalfWave<<" HalfWaveNumber="<<HalfWaveNumber
			 <<std::endl;
		*/
		//<<<JT

//std::cout << "\tvmap: (" << vmap[0] << ", " << vmap[1] << ", " << vmap[2] << ")" << std::endl;

		const G4double HV_value = get_HV_value(Pe, gap1);

		int Pe_fan = 0;
		const G4double dte = elc()->DistanceToTheNearestFan(Pe, Pe_fan);
		G4int side = signof(dte);
		int Pa_fan = 0;
		const G4double dta = lwc()->DistanceToTheNearestFan(Pa, Pa_fan);
//std::cout << "\tdte: " << dte << ", dta: " << dta << std::endl;
//std::cout << "\tPe: (" << Pe.x() << ", " << Pe.y() << ", " << Pe.z() << ")" << std::endl;
//std::cout << "\tPa: (" << Pa.x() << ", " << Pa.y() << ", " << Pa.z() << ")" << std::endl;

		const G4double suppression_range = ElectrodeFanHalfThickness + CHC_Esr;
		if(fabs(dte) < suppression_range){
//			std::cout << " S";
			continue; //skip point if too close to the electrode
		}

		const G4double agap = fabs(dte) - ElectrodeFanHalfThickness
		              + fabs(dta) - lwc()->GetFanHalfThickness();   //correction to electrode suppression not to
		G4double suppression = agap / (agap - CHC_Esr); // change av. signal in the gap
		if(suppression < 0.) suppression = 1.;
//std::cout << "\tagap: " << agap << ", suppression: " << suppression << std::endl;

		//search for radial layers the substep is in between;
/*		G4int irlayer = 0;
		const G4int &n_layers = ChCollWheelType->NumberOfRadialLayer;
		G4double R = Pe.r();
		for(G4int j = 1; j < n_layers; ++ j){
			if(R < ChCollWheelType->RadiusOfLayers[j]){
				irlayer = j - 1;
				break;
			}
		}*/
	// AMS
		G4int irlayer = 0;
		const G4int &n_layers = ChCollWheelType->NumberOfRadialLayer;


		//JT>>>>>G4double fold_angle = M_PI - 2. * parameterized_slant_angle(Pe.r());

		const G4double fold_angle = M_PI - 2. * lwc()->parameterized_slant_angle(rforalpha); //<<<JT
		const G4double fold_angle_deg=fold_angle/M_PI*180.;  //<<<JT

		for(G4int j = 1; j < n_layers; ++ j){
			if(fold_angle_deg > ChCollWheelType->FoldinAngleOfLayers[j]){
				irlayer = j - 1;
				break;
			}
		}
//std::cout << "\tirlayer: " << irlayer << ", R " << R << std::endl;
		assert(irlayer >= 0 && irlayer < n_layers - 1);

		G4double cur = 1.;

		if(!ChCollFoldType->FieldMapPrepared // fieldmap is not available for this fold
		|| (   vmap[2] < ChCollFoldType->MinZofLayer[irlayer] // or out of z range of both maps
		    && vmap[2] < ChCollFoldType->MinZofLayer[irlayer + 1])
		){
			// interpolation not needed
//std::cout << "\tinterpolation not needed" << std::endl;
		} else {        // interpolation needed
			// compute radial interpolation parameter
/*			G4double dr = (R - ChCollWheelType->RadiusOfLayers[irlayer])
				/ (ChCollWheelType->RadiusOfLayers[irlayer + 1] - ChCollWheelType->RadiusOfLayers[irlayer]);
*/
	// AMS  
		  // JT>>
			G4double dr = -(ChCollWheelType->FoldinAngleOfLayers[irlayer] - fold_angle_deg)  
				/ (ChCollWheelType->FoldinAngleOfLayers[irlayer + 1] - ChCollWheelType->FoldinAngleOfLayers[irlayer]);
		
			//std::cout<<"**** GetCurrent1 dr="<<dr<<" n_layers="<< n_layers<<" irlayer= "<<irlayer<<std::endl; 
		  //<<JT

			if(dr < 0.) dr = 0.;
			if(dr > 1.) dr = 1.;

//std::cout << "\tdr " << dr << " <- ";
//std::cout << R << " " << ChCollWheelType->RadiusOfLayers[irlayer]
//          << " " << ChCollWheelType->RadiusOfLayers[irlayer+1]
//		  << " " << ChCollWheelType->RadiusOfLayers[irlayer] << std::endl;

			const G4double a_e = elc()->AmplitudeOfSurface(Pe, side, Pe_fan);
			const G4double a_a = lwc()->AmplitudeOfSurface(Pa, signof(dta), Pa_fan);
//std::cout << "\ta_e: " << a_e << ", a_a: " << a_a << std::endl;

			const G4double x_e = fabs(Pe.x() - a_e);
			const G4double x_a = fabs(Pa.x() - a_a);

			//check geom. err condition if point is outside of LAr gap
			if(gaperr > -100 && (x_e < 0. || x_a < 0.)) gaperr -= 100;


			const G4double z1 = Pe.z() - lwc()->GetStraightStartSection();
			G4int nhwave = G4int(z1 / lwc()->GetHalfWaveLength());
			if(nhwave < 0) nhwave = 0;
			if(nhwave >= lwc()->GetNumberOfHalfWaves()) nhwave = lwc()->GetNumberOfHalfWaves() - 1;
			if((nhwave % 2) == 1) side = -side;

			// get relative y coordinate
			G4double yratio = ((side > 0)? x_e: x_a) / (x_e + x_a);
//std::cout << "\tyratio: " << yratio << std::endl;

			if(yratio <= 0.) yratio = 0.00001;      // pull the point into the gap if it wouldn't be there;
			else if(yratio >= 1.) yratio = 0.99999; // this may happen bec.G4 does not grantee that the full
                                                    // step-line is within the same volume

			//std::cout<<"**GetCurrent1: side="<<side<<" x_a,x_e="<<x_a<<" "<<x_e<<" yratio="<<
			//  yratio<<std::endl;

// get corresponding y coordinates on the radial layers
// where the weight is to be taken from;
			G4double shift = lwc()->GetStraightStartSection();
			if(PointFoldMapArea != 0) shift += WaveLength;
			SetHalfWave(shift + vmap[2]);
//std::cout << "\tPointFoldMapArea: " << PointFoldMapArea << ", SetHalfWave(" << shift + vmap[2] << ")" << std::endl;

			SetYlimitsofPhigapinFieldMap(irlayer);    //on the lower layer
//std::cout << "\tFieldmap limits: " << Ylimits[0] << " " << Ylimits[1] << " " << Ylimits[2] << " " << Ylimits[3] << std::endl;
			G4double pos_low = 0.;
			if(side < 0) pos_low = Ylimits[0] * (1. - yratio) + Ylimits[1] * yratio;
			else if(side > 0) pos_low = Ylimits[2] * (1. - yratio) + Ylimits[3] * yratio;

			SetYlimitsofPhigapinFieldMap(irlayer + 1);  //on the upper layer
			G4double pos_up = 0.;
			if(side < 0) pos_up = Ylimits[0] * (1. - yratio) + Ylimits[1] * yratio;
			else if(side > 0) pos_up = Ylimits[2] * (1. - yratio) + Ylimits[3] * yratio;

// get weights from the maps of lower and upper layer
			const G4double w_low = GetWeightfromFieldMap(irlayer, vmap[2], pos_low);
			const G4double w_up = GetWeightfromFieldMap(irlayer + 1, vmap[2], pos_up);

// get the interpolated normalized Efield alias cur
			cur = w_low * (1. - dr) + w_up * dr;
		}// end of interpolation

//(25-05-2005) new current calculation: edep*1/U*IonReco*E*v_drift
//   normalized so that signal=edep in the straight section of the same gap
//                             as above;
//       where: HV=1250V,E=9.615 [kv/cm] ,T=88.16 K,
//       vdrift=4.62 [mm/mikrosec],  Ionreco=96.4%
// ==> AverageCurrent=3.425/[mikrosec];

//		G4double gap = HalfLArGapSize(Pe.r(), Pe.r());       // gapsize in the straight section
// AMS: why "straight section"? I would better use agap...
// JT: because the fieldmap table is normalized to the efield in the straight section

		// JT>>const G4double &gap = agap;

		const G4double efield = cur * (HV_value * CLHEP::volt) / (gap * CLHEP::mm) / (CLHEP::kilovolt / CLHEP::cm) /cylgapcorr;
		//	efield = efield/cylgapcorr;  //<<JT // DM division by cylgapcorr combined to prev. line
		const G4double substep_current =  step_current * cur / gap * IonReco(efield) *
				   DriftVelo(LArTemperature_av, efield) * suppression / cylgapcorr;
		//substep_current = substep_current / cylgapcorr; //  <<JT // DM division by cylgapcorr combined to prev. line
 
//		if(birksLaw){
//			current += (*birksLaw)(substep_current, step/nofstep/cm, efield);
//		} else {
			current += substep_current;
//		}
	} // end of loop for substeps
//std::cout << std::endl;
	if(current < 0.){
		gaperr -= 1000;
		current = 0.;
	}

	if(gaperr != 0){
		CHCEbad = CHCEbad + edep;
		if(CHCIprint < CHCMaxPrint){
			++ CHCIprint;
			(*m_msg) << MSG::WARNING
			       << "GetCurrent has strange step, gaperr="
				   << gaperr << " correction still computed"
				   << " bad edep ratio=" << CHCEbad/CHCEtotal
				   << endreq;
		}
	}

	//std::cout<<"GetCurrent1::edep="<<edep<<" current="<<current <<" gaperr="<<gaperr<<std::endl;
	return current;
}


G4double EnergyCalculator::get_HV_value(
		       const G4ThreeVector& p, const std::pair<G4int, G4int> &gap) const
{
	const G4int atlas_side = (lwc()->GetAtlasZside() > 0) ? 0 : 1;

	G4ThreeVector p1 ( p );
	p1[2] += lwc()->GetElecFocaltoWRP() + lwc()->GetdWRPtoFrontFace();
	const G4double eta = p1.pseudoRapidity();
	G4int eta_section = -1;
	for(G4int i = 1; i <= NofEtaSection; ++ i){
	  if(eta <= HV_Etalim[i]){
	    eta_section = i - 1;
	    break;
	  }
	}
	if(!(eta_section>=0 && eta_section <=NofEtaSection-1)) throw std::runtime_error("Index out of range");	

	//assert(eta_section >= 0 && eta_section < NofEtaSection);

                  /*(right side of e large phi)*/   /*left side of electrode(small phi)*/
	const G4int e_side = (gap.second > 0) ?   1   :    0;              

	const G4int first_electrode = HV_Start_phi[atlas_side][eta_section][e_side];

	if(first_electrode < 0 || first_electrode >= lwc()->GetNumberOfFans()){
		(*m_msg) << MSG::FATAL
		       << " get_HV_value: first_electrode number is out of range"
			   << endreq;
		G4Exception("EnergyCalculator", "ElectrodeOutOfRange", FatalException, 
			    "get_HV_value: first_electrode number is out of range");
	}

	G4int e_index = lwc()->PhiGapNumberForWheel(gap.first) - first_electrode;
	if(e_index < 0) e_index += lwc()->GetNumberOfFans();
	const G4int &phi_section = e_index;
/*	G4int nofelectrodesinphisection = lwc()->GetNumberOfFans() / NofPhiSections;//24(8) for outer(inner) wheel
	G4int phi_section = e_index / nofelectrodesinphisection;

	if(phi_section < 0 || phi_section >= NofPhiSections){
		(*m_msg) << MSG::FATAL
		       << " get_HV_value: phi_section number is out of range"
			   << endreq;
		G4Exception("EnergyCalculator::get_HV_value");
	}*/
/*
printf("HV: %.2f at eta %.4f, phi %.4f",
       HV_Values[atlas_side][eta_section][e_side][phi_section], eta, p.phi());
std::cout << " " << gap.first << " | "
		  << atlas_side << " " << eta_section << " " << e_side << " " << phi_section
		  << std::endl;
*/
	return HV_Values[atlas_side][eta_section][e_side][phi_section];
}

/*
G4double EnergyCalculator::transform_z_to_fieldmap(const G4ThreeVector &p)
{
  // fieldmap is defined for
  //  --  'normal' folds  :
  //      in the coord syst. of a halfwave of pos.slope
  //      in the range of 0<z<QuarterWaveLength;
  //  --  0th and last fold in the range -13mm < z < QuarterWaveLength/2.
  // SetHalfWave() and GetPhiGap() should be called previously;

	G4double result = 0.;

	if(PointFoldMapArea == 0){
		if(HalfWaveNumber == 0) result = ZinHalfWave; // 0th fold
		else result = -ZinHalfWave; // last fold
	} else {    // normal folds
		result = fabs(ZinHalfWave);
	}
	return result;
}

*/
