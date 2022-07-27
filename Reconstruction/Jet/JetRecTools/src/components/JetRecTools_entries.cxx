#include "JetRecTools/JetTrackSelectionTool.h"
#include "JetRecTools/JetTrackSelectionTool2.h"
#include "JetRecTools/SimpleJetTrackSelectionTool.h"
#include "JetRecTools/TrackVertexAssociationTool.h"
#include "JetRecTools/JetConstituentModSequence.h"
#include "JetRecTools/JetConstituentModifierBase.h"
#include "JetRecTools/CaloClusterConstituentsOrigin.h"
#include "JetRecTools/SoftKillerWeightTool.h"
#include "JetRecTools/VoronoiWeightTool.h"
#include "JetRecTools/ClusterAtEMScaleTool.h"
#include "JetRecTools/ConstitTimeCutTool.h"
#include "JetRecTools/ConstituentSubtractorTool.h"
#include "JetRecTools/JetInputElRemovalTool.h"
#include "JetRecTools/CorrectPFOTool.h"
#include "JetRecTools/ChargedHadronSubtractionTool.h"
#include "JetRecTools/PuppiWeightTool.h"
#include "JetRecTools/JetTrackSelectionAlg.h"
#include "JetRecTools/JetPFlowSelectionAlg.h"
#include "JetRecTools/JetTrackVtxAssoAlg.h"

#if !defined(XAOD_ANALYSIS)
#include "JetRecTools/JetUsedInFitTrackDecoratorTool.h"
#endif


DECLARE_COMPONENT( JetTrackSelectionTool )
DECLARE_COMPONENT( JetTrackSelectionTool2 )
DECLARE_COMPONENT( SimpleJetTrackSelectionTool )
DECLARE_COMPONENT( TrackVertexAssociationTool )
DECLARE_COMPONENT( JetConstituentModSequence )
DECLARE_COMPONENT( JetConstituentModifierBase )
DECLARE_COMPONENT( CaloClusterConstituentsOrigin )
DECLARE_COMPONENT( SoftKillerWeightTool )
DECLARE_COMPONENT( VoronoiWeightTool )
DECLARE_COMPONENT( ClusterAtEMScaleTool )
DECLARE_COMPONENT( ConstitTimeCutTool )
DECLARE_COMPONENT( ConstituentSubtractorTool )
DECLARE_COMPONENT( JetInputElRemovalTool )
DECLARE_COMPONENT( CorrectPFOTool )
DECLARE_COMPONENT( ChargedHadronSubtractionTool )
DECLARE_COMPONENT( PuppiWeightTool )
DECLARE_COMPONENT( JetTrackSelectionAlg )
DECLARE_COMPONENT( JetPFlowSelectionAlg )
DECLARE_COMPONENT( JetTrackVtxAssoAlg )

#if !defined(XAOD_ANALYSIS)
DECLARE_COMPONENT( JetUsedInFitTrackDecoratorTool )
#endif
