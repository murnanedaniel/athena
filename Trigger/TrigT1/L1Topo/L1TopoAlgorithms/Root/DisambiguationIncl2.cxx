/*********************************
 * DisambiguationIncl2.cpp
 * Created by Joerg Stelzer / V Sorin on 2014.
 * Copyright (c) 2012 Joerg Stelzer. All rights reserved.
 *
 * @brief algorithm calculates the dR distance between objects in two lists, accept if at least a pair does not match in dR 
 *
 * @param NumberLeading
**********************************/

#include <cmath>

#include "L1TopoAlgorithms/DisambiguationIncl2.h"
#include "L1TopoCommon/Exception.h"
#include "L1TopoInterfaces/Decision.h"

REGISTER_ALG_TCS(DisambiguationIncl2)

using namespace std;

// not the best solution but we will move to athena where this comes for free
#define LOG cout << "TCS::DisambiguationIncl2:     "




namespace {
   unsigned int
   calcDeltaR2(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2) {
      double deta = ( tob1->etaDouble() - tob2->etaDouble() );
      double dphi = fabs( tob1->phiDouble() - tob2->phiDouble() );
      if(dphi>M_PI)
         dphi = 2*M_PI - dphi;

      return round ( 100 * ((dphi)*(dphi) + (deta)*(deta) )) ;

   }

   unsigned int
   calcDeltaR2BW(const TCS::GenericTOB* tob1, const TCS::GenericTOB* tob2) {

      int detaB = abs( tob1->eta() - tob2->eta() );
      int dphiB = abs( tob1->phi() - tob2->phi() );
      if(dphiB>32)
         dphiB = 64 - dphiB;

      unsigned int bit_dr2 = dphiB*dphiB + detaB*detaB;
      return bit_dr2;

   }
  
}


TCS::DisambiguationIncl2::DisambiguationIncl2(const std::string & name) : DecisionAlg(name)
{
   defineParameter("InputWidth1", 9);
   defineParameter("InputWidth2", 9);
   defineParameter("MaxTob1", 0); 
   defineParameter("MaxTob2", 0); 
   defineParameter("NumResultBits", 2);
   defineParameter("ClusterOnly",0);
   defineParameter("ApplyDR",0);
   defineParameter("MinET1",1,0);
   defineParameter("MinET2",1,0);
   defineParameter("DisambDRSqrMin",0,0);
   defineParameter("MinET1",1,1);
   defineParameter("MinET2",1,1);
   defineParameter("DisambDRSqrMin",0,1);
   setNumberOutputBits(2);
}

TCS::DisambiguationIncl2::~DisambiguationIncl2(){}


TCS::StatusCode
TCS::DisambiguationIncl2::initialize() {
   p_NumberLeading1 = parameter("InputWidth1").value();
   p_NumberLeading2 = parameter("InputWidth2").value();
   if(parameter("MaxTob1").value() > 0) p_NumberLeading1 = parameter("MaxTob1").value();
   if(parameter("MaxTob2").value() > 0) p_NumberLeading2 = parameter("MaxTob2").value();

   for(unsigned int i=0; i<numberOutputBits(); ++i) {
      p_MinET1[i] = parameter("MinET1",i).value();
      p_MinET2[i] = parameter("MinET2",i).value();

      p_DisambDR[i] = parameter("DisambDRSqrMin", i).value();
   
      TRG_MSG_INFO("MinET1          " << i << " : "<< p_MinET1[i]);
      TRG_MSG_INFO("MinET2          " << i << " : "<< p_MinET2[i]);
      TRG_MSG_INFO("DisambDR         " << i << " : "<< p_DisambDR[i]);


   }


   TRG_MSG_INFO("number output : " << numberOutputBits());


   return StatusCode::SUCCESS;
}



TCS::StatusCode
TCS::DisambiguationIncl2::processBitCorrect( const std::vector<TCS::TOBArray const *> & input,
                             const std::vector<TCS::TOBArray *> & output,
                             Decision & decision )
{

      
   if( input.size() == 2) {

      
      for( TOBArray::const_iterator tob1 = input[0]->begin(); 
           tob1 != input[0]->end() && distance(input[0]->begin(), tob1) < p_NumberLeading1;
           ++tob1)
         {


            for( TCS::TOBArray::const_iterator tob2 = input[1]->begin(); 
                 tob2 != input[1]->end() && distance(input[1]->begin(), tob2) < p_NumberLeading2;
                 ++tob2) {


               // test DeltaR2Min, DeltaR2Max
               unsigned int deltaR2 = calcDeltaR2BW( *tob1, *tob2 );
               
               for(unsigned int i=0; i<numberOutputBits(); ++i) {
                   bool accept = false;
                   if( parType_t((*tob1)->Et()) <= p_MinET1[i]) continue; // ET cut
                   if( parType_t((*tob2)->Et()) <= p_MinET2[i]) continue; // ET cut
                   accept = deltaR2 > p_DisambDR[i] ;
                   if( accept ) {
                       decision.setBit(i, true);
                       output[i]->push_back(TCS::CompositeTOB(*tob1, *tob2));
                   }
                   TRG_MSG_DEBUG("Decision " << i << ": " << (accept?"pass":"fail") << " deltaR2 = " << deltaR2);
               }
            }
         }
   } else {

      TCS_EXCEPTION("DisambiguationIncl2 alg must have  2 inputs, but got " << input.size());

   }
   return TCS::StatusCode::SUCCESS;

}

TCS::StatusCode
TCS::DisambiguationIncl2::process( const std::vector<TCS::TOBArray const *> & input,
                             const std::vector<TCS::TOBArray *> & output,
                             Decision & decision )
{

      
   if( input.size() == 2) {

      
      for( TOBArray::const_iterator tob1 = input[0]->begin(); 
           tob1 != input[0]->end() && distance(input[0]->begin(), tob1) < p_NumberLeading1;
           ++tob1)
         {


            for( TCS::TOBArray::const_iterator tob2 = input[1]->begin(); 
                 tob2 != input[1]->end() && distance(input[1]->begin(), tob2) < p_NumberLeading2;
                 ++tob2) {


               // test DeltaR2Min, DeltaR2Max
               unsigned int deltaR2 = calcDeltaR2( *tob1, *tob2 );
               
               for(unsigned int i=0; i<numberOutputBits(); ++i) {
                   bool accept = false;
                   if( parType_t((*tob1)->Et()) <= p_MinET1[i]) continue; // ET cut
                   if( parType_t((*tob2)->Et()) <= p_MinET2[i]) continue; // ET cut
                   accept = deltaR2 > p_DisambDR[i] ;
                   if( accept ) {
                       decision.setBit(i, true);
                       output[i]->push_back(TCS::CompositeTOB(*tob1, *tob2));
                   }
                   TRG_MSG_DEBUG("Decision " << i << ": " << (accept?"pass":"fail") << " deltaR2 = " << deltaR2);
               }
            }
         }
   } else {

      TCS_EXCEPTION("DisambiguationIncl2 alg must have  2 inputs, but got " << input.size());

   }
   return TCS::StatusCode::SUCCESS;
}
