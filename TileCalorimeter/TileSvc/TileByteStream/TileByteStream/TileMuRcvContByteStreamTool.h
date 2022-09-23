/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//****************************************************************************
// Filename : TileMuRcvContByteStreamTool.h
// Author   : Joao Gentil Saraiva (jmendes@cern.ch)
// Created  : October 2014
//
// DESCRIPTION
//    AlgTool class to provide conversion from TileMuRcvContainer to ByteStream
//    and fill it in RawEvent
//
// BUGS:
//
// History:
//
//****************************************************************************

#ifndef TILEBYTESTREAM_TILEMURCVCONTBYTESTREAMTOOL_H
#define TILEBYTESTREAM_TILEMURCVCONTBYTESTREAMTOOL_H

#include "ByteStreamCnvSvcBase/FullEventAssembler.h"
#include "StoreGate/ReadCondHandleKey.h"

#include "AthenaBaseComps/AthAlgTool.h"
#include "TileEvent/TileMuonReceiverObj.h"
#include "TileEvent/TileMuonReceiverContainer.h"
#include "TileByteStream/TileHid2RESrcID.h"

class TileHWID;

#include <string>

/**
 * @class TileMuRcvContByteStreamTool
 * @brief This AlgTool class provides conversion from TileMuonReceiverContainer to ByteStream and fill it in RawEvent
 * @author Joao Gentil Saraiva
 **/

class TileMuRcvContByteStreamTool: public AthAlgTool {

 public:

  /** Constructor */
  TileMuRcvContByteStreamTool( const std::string& type, const std::string& name, const IInterface* parent );

  /** Destructor */
  virtual ~TileMuRcvContByteStreamTool();

  /** AlgTool InterfaceID */
  static const InterfaceID& interfaceID( );

  virtual StatusCode initialize() override;
  virtual StatusCode finalize() override;

  /** Provides conversion from TileMuRcvContainer to bytestream */
  StatusCode convert(TileMuonReceiverContainer* cont, FullEventAssembler<TileHid2RESrcID> *fea) const;

 private:

  Gaudi::Property<bool> m_initializeForWriting{this, "InitializeForWriting", false, "Initialize for writing"};

  SG::ReadCondHandleKey<TileHid2RESrcID> m_hid2RESrcIDKey{this,
     "TileHid2RESrcID", "TileHid2RESrcIDHLT", "TileHid2RESrcID key"};

  const TileHWID* m_tileHWID;

  int m_runPeriod;
};

#endif
