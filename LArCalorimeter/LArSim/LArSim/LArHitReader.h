/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/* Generated by Together */

/* author : Johann Collot */

/* date of creation : 31/01/2001 */
/* modifications :                                                                              */
/*  10/06/2001 :                                                                                */
/*  13/09/2002 : corrections of comments                                                        */
/*  04/10/2002 : corrections of comments                                                        */

#ifndef LARHITREADER_H
#define LARHITREADER_H

#include <string>

#include "GaudiKernel/Algorithm.h"
#include "LArSimEvent/LArHitContainer.h"

class LArHitReader : public Algorithm {

public:

    /** usual ATHENA constructor of an algorithm*/
    LArHitReader(const std::string & name, ISvcLocator * pSvcLocator);

    /** Destructor */
    virtual ~LArHitReader() { }

    /** Initialize method , executed once at the beginning <br> of the job execution by the control framework */
    virtual StatusCode initialize();

    /** execute() method executed for each event <br> by the LArHitMaker instance which owns <br> this producer */
    virtual StatusCode execute();

    /** finalize() method executed at the end <br> of the job execution by the control framework */
    virtual StatusCode finalize();

private:


};
#endif    //LARHITPRODUCER_H


