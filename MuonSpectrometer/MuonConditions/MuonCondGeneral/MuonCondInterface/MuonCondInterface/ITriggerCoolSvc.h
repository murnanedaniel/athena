/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONCONDSVC_ITRIGGERCOOLSVC_H
#define MUONCONDSVC_ITRIGGERCOOLSVC_H

#include <string>

#include "CxxUtils/checker_macros.h"
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Service.h"

namespace MuonCalib {

    class ATLAS_NOT_THREAD_SAFE ITriggerCoolSvc : virtual public IInterface {
    public:
        static const InterfaceID& interfaceID();

        // put data in file into COOL

        virtual StatusCode writeToDBEta(const std::string& etaFolder, const std::string& filename_CM, const std::string& filename_Th0,
                                        const int chan_Eta, const std::string& sequence, const std::string& infoEta) const = 0;
        virtual StatusCode writeToDBPhi(const std::string& phiFolder, const std::string& filename_CM_Phi,
                                        const std::string& filename_Th0_Phi, const int chan_Phi, const std::string& infoPhi) const = 0;

        virtual StatusCode getData(const std::string& etaFolder, const int chan) const = 0;
    };

    inline const InterfaceID& ITriggerCoolSvc::interfaceID() {
        static const InterfaceID IID_ITriggerCoolSvc("ITriggerCoolSvc", 1, 0);
        return IID_ITriggerCoolSvc;
    }
}  // namespace MuonCalib
#endif  // MUONCONDSVC_ITRIGGERCOOLSVC_H
