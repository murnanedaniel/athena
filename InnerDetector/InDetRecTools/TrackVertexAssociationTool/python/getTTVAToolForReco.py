# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# These are *constants* (OK, the user actually has the ability to change the track/vertex containers in the function calls,
# but the default container names are constants)
_DEFAULT_TRACK_CONT  = 'InDetTrackParticles'
_DEFAULT_VERTEX_CONT = 'PrimaryVertices'
_DECO_TOOL_NAME      = 'InDetUsedInFitDecoratorTool'
_DECO_ALG_NAME       = 'InDetUsedInFitDecorator'
_VERTEX_DECO         = 'TTVA_AMVFVertices_forReco'
_WEIGHT_DECO         = 'TTVA_AMVFWeights_forReco'

def _assertDecoValue(kwargs, key, default):
    if key in kwargs.keys():
        if kwargs[key] != default:
            raise RuntimeError('Property \'%s\' must be \'%s\': %s = %s' % (key, default, key, kwargs[key]))
    else:
        kwargs[key] = default


def addUsedInFitDecoratorForReco(TrackContName = _DEFAULT_TRACK_CONT, VertexContName = _DEFAULT_VERTEX_CONT, add2Seq = None):
    '''Function for adding a used-in-fit decorator alg to a sequence'''
    from AthenaCommon.AlgSequence import AlgSequence
    from InDetUsedInVertexFitTrackDecorator.UsedInVertexFitTrackDecoratorCfg import idForTrackVtxContainers

    topSequence = AlgSequence()
    # If not sequence is provided, add the decoration alg to the top sequence
    if add2Seq is None:
        add2Seq = topSequence

    # the alg name must be unique for a given trk,vtx and deco combination. Otherwise multiple algs
    # decorating identically may be scheduled (this function is called from many places) and results in errors (ATLASRECTS-6925)
    # the function idForTrackVtxContainers( ) returns a unique name for such a combination.
    basename = idForTrackVtxContainers(TrackContName, VertexContName, _VERTEX_DECO, _WEIGHT_DECO)
    if not hasattr(add2Seq, basename+"Alg"):
        import AthenaCommon.CfgMgr as CfgMgr
        from AthenaCommon.AppMgr import ToolSvc
        ToolSvc += CfgMgr.InDet__InDetUsedInFitTrackDecoratorTool(  name                    = basename+"Tool",
                                                                    AMVFVerticesDecoName    = _VERTEX_DECO,
                                                                    AMVFWeightsDecoName     = _WEIGHT_DECO,
                                                                    TrackContainer          = TrackContName,
                                                                    VertexContainer         = VertexContName )
        add2Seq += CfgMgr.InDet__InDetUsedInVertexFitTrackDecorator(name                    = basename+"Alg",
                                                                    UsedInFitDecoratorTool  = getattr(ToolSvc, basename+"Tool") )
    return

def getTTVAToolForReco(name = None, returnCompFactory = False, addDecoAlg = True, add2Seq = None, **kwargs):
    '''Function for retrieving a TTVA tool instance and, if it doesn't exist,
    for adding a used-in-fit track decorator to a sequence'''

    # First some argument checking
    _assertDecoValue(kwargs, 'AMVFVerticesDeco', _VERTEX_DECO)
    _assertDecoValue(kwargs, 'AMVFWeightsDeco',  _WEIGHT_DECO)

    # If a name has been provided for the TTVA tool, update our kwargs
    if name is not None:
        kwargs['name'] = name

    # Fetch the track and vertex containers for use in the used-in-fit deco tool
    TrackContName, VertexContName = kwargs.get('TrackContName', _DEFAULT_TRACK_CONT), kwargs.pop('VertexContName', _DEFAULT_VERTEX_CONT)

    # Instantiate our TTVA tool
    if returnCompFactory:
        from AthenaConfiguration.ComponentFactory import CompFactory
        tool = CompFactory.getComp("CP::TrackVertexAssociationTool")(**kwargs)
    else:
        import AthenaCommon.CfgMgr as CfgMgr
        tool = CfgMgr.CP__TrackVertexAssociationTool(**kwargs)

    # Add the decorator to the top sequence (if requested)
    if addDecoAlg:
        addUsedInFitDecoratorForReco(TrackContName, VertexContName, add2Seq)

    return tool
