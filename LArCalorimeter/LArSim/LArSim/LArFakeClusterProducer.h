/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/* Generated by Together */

#ifndef LARFAKECLUSTERPRODUCER_H
#define LARFAKECLUSTERPRODUCER_H
#include "LArSim/LArHitProducer.h"

/* author : Johann Collot */

/* date of creation : 31/01/2001 */
/* date of last modification : 13/09/2002 -  create fake clusters in EMB */

/**
 * Example of creation of LArHits through a LArHitProducer . <br><br> 
 * The only method which has to be implement is hitconstruction() .<br>
 * @author Johann Collot
 * @version 00-00-03
 * 
 */
class LArFakeClusterProducer : public LArHitProducer {

public:

    /** constructor */
    LArFakeClusterProducer(const std::string & name, ISvcLocator * pSvcLocator) ;

    /** destructor */
    virtual ~LArFakeClusterProducer() { }

    /** method in which the specific code <br> to produce the hits is placed */
    virtual StatusCode hitConstruction();
    
private:

    /** Total Energy of created Hits <br> 
      * This is a property of this algorithm <br>
      * Default value : 100 GeV 
      */
    double m_Etot ;
    
    /** Number of generated Clusters in eta <br> 
      * Calculated from ClusterSpacing and EtaStart 
      */
    int m_NumberOfClusters ;    

    /** number of cells in eta between 2 clusters <br> 
      * This is a property of this algorithm <br>
      * Default value : 4
      */
    int m_ClusterSpacing ;    

    /** Eta cell of first cluster <br> 
      * This is a property of this algorithm <br>
      * Default value : 3
      */
    int m_EtaStart ;

    /**
     * Sampling Fraction of EMB presampler <br>
     * Property of this algorithm <br>
     * Default value : 0.2 
     */
    double m_SampFracEMBPS;

    /**
     * Sampling Fraction of EMB first compartment  <br>
     * Property of this algorithm <br>
     * Default value : 0.2
     */
    double m_SampFracEMB1;

    /**
     * Sampling Fraction of EMB middle compartment <br>
     * Property of this algorithm <br>
     * Default value : 0.2
     */
    double m_SampFracEMB2;

    /**
     * Sampling Fraction of EMB back compartment <br>
     * Property of this algorithm <br>
     * Default value : 0.2
     */
    double m_SampFracEMB3;
};
#endif  //LARFAKECLUSTERPRODUCER_H
