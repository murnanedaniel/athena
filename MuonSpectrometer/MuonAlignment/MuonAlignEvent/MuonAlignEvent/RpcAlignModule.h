/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONALIGNGENALGS_RPCALIGNMODULE_H
#define MUONALIGNGENALGS_RPCALIGNMODULE_H

#include <iostream>
#include <vector>

#include "MuonAlignEvent/CombinedMuonAlignModule.h"
#include "TrkAlignEvent/AlignModule.h"
#include "TrkAlignEvent/AlignPar.h"

/**
   @file RpcAlignModule.h
   @class RpcAlignModule

   @brief RpcAlignModule is a grouping of RpcReadoutElements, grouped
   according to the type of alignment.

   @author Robert Harrington <roberth@bu.edu>
   @date 06/09/09
*/

class TFile;
class TTree;

namespace MuonGM {
    class RpcReadoutElement;
}

namespace Muon {

    class RpcAlignModule : public CombinedMuonAlignModule {
    public:
        /** Constructor using AlgTool gives a MsgStream with RpcAlignModule for name and the the same output level as the AlgTool. */
        RpcAlignModule(const AlgTool* algtool, const Amg::Transform3D& transform = Amg::Transform3D::Identity());

        virtual ~RpcAlignModule();

        // void   shiftSurface(Trk::TrkDetElementBase* det, Identifier id) const;
        // void   restoreSurfaces(Trk::TrkDetElementBase* tre) const;

    private:
        MsgStream* m_log;
        RpcAlignModule& operator=(const RpcAlignModule& right);
        RpcAlignModule(const RpcAlignModule&);

    };  // end class

}  // namespace Muon

#endif  // MUONALIGNGENALGS_RPCALIGNMODULE_H
