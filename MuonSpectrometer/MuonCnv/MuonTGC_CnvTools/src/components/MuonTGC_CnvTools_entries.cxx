#include "../TgcRdoToPrepDataToolMT.h"
#include "../TgcRdoContByteStreamTool.h"
#include "../TGC_RodDecoderReadout.h"
#include "../TGC_RodDecoderRawdata.h"
#include "../TGC_RawDataProviderTool.h"
#include "../TGC_RawDataProviderToolMT.h"
#include "../TgcRDO_Decoder.h"
#include "../TgcPrepDataReplicationAlg.h"
#include "../TgcPrepDataReplicationTool3BCtoAllBC.h"
#include "../TgcPrepDataReplicationToolAllBCto3BC.h"


DECLARE_COMPONENT( Muon::TgcRdoToPrepDataToolMT )
DECLARE_COMPONENT( Muon::TgcRdoContByteStreamTool )
DECLARE_COMPONENT( Muon::TgcPrepDataReplicationTool3BCtoAllBC )
DECLARE_COMPONENT( Muon::TgcPrepDataReplicationToolAllBCto3BC )
DECLARE_COMPONENT( Muon::TGC_RodDecoderReadout )
DECLARE_COMPONENT( Muon::TGC_RodDecoderRawdata )
DECLARE_COMPONENT( Muon::TGC_RawDataProviderTool )
DECLARE_COMPONENT( Muon::TGC_RawDataProviderToolMT )
DECLARE_COMPONENT( Muon::TgcRDO_Decoder )
DECLARE_COMPONENT( Muon::TgcPrepDataReplicationAlg )

