from InDetRecExample.InDetJobProperties import InDetFlags
from InDetRecExample.InDetKeys import InDetKeys
from InDetRecExample.TrackingCommon import getInDetxAODParticleCreatorTool, makePublicTool, makeName

def getCollectionNameIfInFile(coll_type,coll_name) :
    from RecExConfig.AutoConfiguration import IsInInputFile
    if not IsInInputFile(coll_type,coll_name) :
        printfunc ('DEBUG getRecTrackParticleNameIfInFile set %s' % coll_name)
        return coll_name
    else :
        return ""

def passCollectionName(coll_name,condition) :
    return coll_name if condition else ""

def getRecTrackParticleNameIfInFile(coll_name) :
    return getCollectionNameIfInFile('Rec::TrackParticle',coll_name)

def getRecVertexNameIfInFile(coll_name) :
    return getCollectionNameIfInFile('Rec::Vertex',coll_name)


from AthenaCommon.GlobalFlags import globalflags
is_mc = (globalflags.DataSource == 'geant4')

doCreation = ( InDetFlags.doNewTracking() or InDetFlags.doPseudoTracking() or InDetFlags.doLargeD0() or InDetFlags.runLRTReco() \
                   or InDetFlags.doLowPtLargeD0() )  and InDetFlags.doParticleCreation()
doConversion = not InDetFlags.doNewTracking()  and not InDetFlags.doPseudoTracking() and not InDetFlags.doLargeD0() and not InDetFlags.runLRTReco()\
                    and not InDetFlags.doLowPtLargeD0() and InDetFlags.doParticleConversion()

if doCreation:
    printfunc ("Creating xAOD::TrackParticles from Trk::Tracks")
if doConversion:
    printfunc ("Converting Rec::TrackParticles to xAOD::TrackParticles")



# Run the xAOD truth builder for PU if needed
if InDetFlags.doSplitReco()  and is_mc:
    from xAODTruthCnv.xAODTruthCnvConf import xAODMaker__xAODTruthCnvAlg
    xAODTruthCnvPU = xAODMaker__xAODTruthCnvAlg("xAODTruthCnvPU")
    xAODTruthCnvPU.WriteInTimePileUpTruth = False
    xAODTruthCnvPU.WriteAllPileUpTruth = True
    xAODTruthCnvPU.AODContainerName = "GEN_EVENT_PU"
    xAODTruthCnvPU.xAODTruthEventContainerName = "TruthEvents_PU" #output
    xAODTruthCnvPU.xAODTruthPileupEventContainerName = "TruthPileupEvents_PU" #output
    xAODTruthCnvPU.xAODTruthParticleContainerName = "TruthParticles_PU" #output
    xAODTruthCnvPU.xAODTruthVertexContainerName = "TruthVertices_PU" #output
    xAODTruthCnvPU.TruthLinks = "xAODTruthLinks_PU" #output/intermediate
    xAODTruthCnvPU.MetaObjectName = "TruthMetaData_PU" #output
    topSequence += xAODTruthCnvPU

def getTrackCollectionCnvTool(prd_to_track_map=None, suffix="",trt_pid_tool=True, pixel_dedx=True) :
    from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__TrackCollectionCnvTool
    return xAODMaker__TrackCollectionCnvTool("TrackCollectionCnvTool"+suffix,
                                             TrackParticleCreator = getInDetxAODParticleCreatorTool(prd_to_track_map,
                                                                                                    suffix,
                                                                                                    trt_pid_tool=trt_pid_tool,
                                                                                                    pixel_dedx=pixel_dedx))

def getRecTrackParticleContainerCnvTool(prd_to_track_map=None, suffix="") :
    from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__RecTrackParticleContainerCnvTool
    return xAODMaker__RecTrackParticleContainerCnvTool("RecTrackParticleContainerCnvTool"+suffix,
                                                       TrackParticleCreator = getInDetxAODParticleCreatorTool(prd_to_track_map, suffix))

def isValid(name) :
    return name is not None and name != ""

def createTrackParticles(track_in, track_particle_truth_in,track_particle_out, topSequence, prd_to_track_map=None, suffix="",trt_pid_tool=True, pixel_dedx=True) :
    '''
    create algorithm to convert the input tracks into track xAOD track particles.
    @param track_in the name of the input TrackCollection
    @param track_particle_truth_in optional truth track collection to link to
    @param track_particle_out the name of the output xAOD track particle collection
    @param topSequence the sequence to which the algorithm is added
    @param prd_to_track_map None or if shared hits are to be recomputed a PRDtoTrackMap filled by a preceding
           algorithms e.g. a TrackCollectionMerger.
    @param suffix which makes the names of the particle creator tool and sub-tools unique in case a prd_to_track_map
           is provided.
    @param trt_pid_tool if true run TRT_ElectronPidTool for each track and decorate track particle with the results
    @param pixel_dedx   if true run PixelToTPIDTool for each track and decorate track particle with the results
    '''
    if not trt_pid_tool and not pixel_dedx :
        suffix += "NoPID"
    elif not trt_pid_tool :
        suffix += "NoTRTPID"
    elif not pixel_dedx :
        suffix += "NoPixPID"

    if isValid(track_in) and isValid(track_particle_out) :
        from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__TrackParticleCnvAlg
        xAODTrackParticleCnvAlg = xAODMaker__TrackParticleCnvAlg(track_particle_out)
        xAODTrackParticleCnvAlg.xAODContainerName = ""
        xAODTrackParticleCnvAlg.xAODTrackParticlesFromTracksContainerName = track_particle_out
        xAODTrackParticleCnvAlg.TrackContainerName = track_in
        xAODTrackParticleCnvAlg.TrackParticleCreator = getInDetxAODParticleCreatorTool(prd_to_track_map, suffix,trt_pid_tool=trt_pid_tool, pixel_dedx=pixel_dedx)
        xAODTrackParticleCnvAlg.TrackCollectionCnvTool = getTrackCollectionCnvTool(prd_to_track_map, suffix, trt_pid_tool=trt_pid_tool, pixel_dedx=pixel_dedx)
        xAODTrackParticleCnvAlg.AODContainerName = ""
        xAODTrackParticleCnvAlg.AODTruthContainerName = ""
        xAODTrackParticleCnvAlg.ConvertTrackParticles = False
        xAODTrackParticleCnvAlg.ConvertTracks = True
        xAODTrackParticleCnvAlg.AddTruthLink = InDetFlags.doTruth() and is_mc and isValid(track_particle_truth_in)
        xAODTrackParticleCnvAlg.xAODTruthLinkVector =  passCollectionName( 'xAODTruthLinks', InDetFlags.doTruth() and is_mc and isValid(track_particle_truth_in) )
        xAODTrackParticleCnvAlg.TrackTruthContainerName = passCollectionName(track_particle_truth_in,(is_mc and InDetFlags.doTruth()))
        from MCTruthClassifier.MCTruthClassifierBase import MCTruthClassifier
        xAODTrackParticleCnvAlg.MCTruthClassifier = MCTruthClassifier
        topSequence += xAODTrackParticleCnvAlg

def convertTrackParticles(aod_track_particles_in, track_particle_truth_in,track_particle_out, topSequence) :
    if isValid(aod_track_particles_in) and isValid(track_particle_out) :
        from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__TrackParticleCnvAlg
        xAODTrackParticleCnvAlg = xAODMaker__TrackParticleCnvAlg(track_particle_out)
        xAODTrackParticleCnvAlg.xAODContainerName = track_particle_out
        xAODTrackParticleCnvAlg.xAODTrackParticlesFromTracksContainerName = ""
        xAODTrackParticleCnvAlg.TrackContainerName =  ""
        xAODTrackParticleCnvAlg.TrackParticleCreator = getInDetxAODParticleCreatorTool()
        xAODTrackParticleCnvAlg.RecTrackParticleContainerCnvTool = getRecTrackParticleContainerCnvTool()
        xAODTrackParticleCnvAlg.AODContainerName = aod_track_particles_in
        xAODTrackParticleCnvAlg.AODTruthContainerName = passCollectionName(track_particle_truth_in,(is_mc and InDetFlags.doTruth()) )
        xAODTrackParticleCnvAlg.ConvertTrackParticles = True
        xAODTrackParticleCnvAlg.ConvertTracks = False
        xAODTrackParticleCnvAlg.AddTruthLink = InDetFlags.doTruth() and is_mc and isValid(track_particle_truth_in)
        xAODTrackParticleCnvAlg.xAODTruthLinkVector =  passCollectionName( 'xAODTruthLinks', InDetFlags.doTruth() and is_mc and isValid(track_particle_truth_in) )
        xAODTrackParticleCnvAlg.TrackTruthContainerName = ""
        from MCTruthClassifier.MCTruthClassifierBase import MCTruthClassifier
        xAODTrackParticleCnvAlg.MCTruthClassifier = MCTruthClassifier
        topSequence += xAODTrackParticleCnvAlg

if (doCreation or doConversion):# or InDetFlags.useExistingTracksAsInput()) : <---- [XXX JDC Should we included this?
                                #                                                    problems appear when nothing should
                                #                                                    be done but
                                #                                                    useExistinTracksAsInput...
    # [XXX JDC: to deal with the MergedTracks case, the truth collections are
    #           defined in the InputTrackCollectionTruth variable. To be deprecated
    #           if finally there is no need of the special "MergedTrack" name
    if 'InputTrackCollectionTruth' not in dir():
        InputTrackCollectionTruth = InDetKeys.TracksTruth()
    if not InDetFlags.doDBMstandalone():
        if doCreation :
            createTrackParticles(InputTrackCollection, InputTrackCollectionTruth, InDetKeys.xAODTrackParticleContainer(),topSequence,
                                 trt_pid_tool=True, pixel_dedx=True)
            from  InDetPhysValMonitoring.InDetPhysValJobProperties import InDetPhysValFlags
            from  InDetPhysValMonitoring.ConfigUtils import extractCollectionPrefix
            for col in InDetPhysValFlags.validateExtraTrackCollections() :
                prefix=extractCollectionPrefix(col)
                createTrackParticles(col,"", prefix+"TrackParticles", topSequence,
                                     trt_pid_tool=False, pixel_dedx=False)
        if doConversion :
            convertTrackParticles(getRecTrackParticleNameIfInFile(InDetKeys.TrackParticles()),
                                  InDetKeys.TrackParticlesTruth() ,
                                  InDetKeys.xAODTrackParticleContainer(),
                                  topSequence)


    if (InDetFlags.doDBMstandalone() or InDetFlags.doDBM() ) and doCreation :
        # or instead of InDetKeys.DBMTracksTruth()  rather InDetKeys.DBMDetailedTracksTruth() ?
        createTrackParticles( InDetKeys.xAODDBMTrackParticleContainer(), InDetKeys.DBMTracksTruth(), InDetKeys.xAODDBMTrackParticleContainer(),topSequence,
                              trt_pid_tool=False, pixel_dedx=False)

if not InDetFlags.doVertexFinding():
    if (not InDetFlags.doDBMstandalone() and
        not IsInInputFile ('xAOD::VertexContainer', InDetKeys.xAODVertexContainer()) and
        IsInInputFile ('VxContainer', InDetKeys.PrimaryVertices())) :
        if len(getRecVertexNameIfInFile(InDetKeys.PrimaryVertices()))>0 :
            from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__VertexCnvAlg
            xAODVertexCnvAlg = xAODMaker__VertexCnvAlg("VertexCnvAlg")
            xAODVertexCnvAlg.xAODContainerName = InDetKeys.xAODVertexContainer()
            xAODVertexCnvAlg.AODContainerName = InDetKeys.PrimaryVertices()
            xAODVertexCnvAlg.TPContainerName = InDetKeys.xAODTrackParticleContainer()
            topSequence += xAODVertexCnvAlg

    if InDetFlags.doDBMstandalone() or InDetFlags.doDBM():
        if (IsInInputFile ('VxContainer', InDetKeys.PrimaryVertices())) :
            from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__VertexCnvAlg
            xAODVertexCnvAlgDBM = xAODMaker__VertexCnvAlg("VertexCnvAlgDBM")
            xAODVertexCnvAlgDBM.xAODContainerName = InDetKeys.xAODVertexContainer()
            xAODVertexCnvAlgDBM.AODContainerName = InDetKeys.PrimaryVertices()
            xAODVertexCnvAlgDBM.TPContainerName = InDetKeys.xAODDBMTrackParticleContainer()
            topSequence += xAODVertexCnvAlgDBM

#For forward tracks, no separate collection for ITK, since they are already merged
if (InDetFlags.doForwardTracks() and InDetFlags.doParticleCreation()) or doConversion:
    if doCreation :
        createTrackParticles(InDetKeys.ResolvedForwardTracks(), InDetKeys.ResolvedForwardTracksTruth(), InDetKeys.xAODForwardTrackParticleContainer(),topSequence,
                             trt_pid_tool=False, pixel_dedx=False)
    if doConversion :
        convertTrackParticles(getRecTrackParticleNameIfInFile(InDetKeys.ResolvedForwardTrackParticles()),
                              InDetKeys.ResolvedForwardTrackParticlesTruth(),
                              InDetKeys.xAODForwardTrackParticleContainer(),
                              topSequence)

if InDetFlags.doPseudoTracking():
    if doCreation :
        createTrackParticles(InDetKeys.PseudoTracks(), InDetKeys.PseudoTracksTruth(), InDetKeys.xAODPseudoTrackParticleContainer(),topSequence,
                             trt_pid_tool=False, pixel_dedx=False) # @TODO PID for pseudo particles ? 
    if doConversion :
        convertTrackParticles(getRecTrackParticleNameIfInFile(InDetKeys.PseudoTracks()),
                              InDetKeys.TrackParticlesTruth(),
                              InDetKeys.xAODPseudoTrackParticleContainer(),
                              topSequence)


if InDetFlags.runLRTReco() and InDetFlags.storeSeparateLargeD0Container():
    if doCreation :
        createTrackParticles(InDetKeys.ExtendedLargeD0Tracks(), InDetKeys.ExtendedLargeD0TracksTruth(), InDetKeys.xAODLargeD0TrackParticleContainer(),topSequence,
                             trt_pid_tool=False, pixel_dedx=False) # @TODO TRT PID for LRT ?


if InDetFlags.doTrackSegmentsPixel() and InDetFlags.doParticleCreation():
    if doCreation :
        createTrackParticles(InDetKeys.PixelTracks(), InDetKeys.PixelTracksTruth(), InDetKeys.xAODPixelTrackParticleContainer(),topSequence,
                             trt_pid_tool=False, pixel_dedx=False)
    if doConversion :
        convertTrackParticles(getRecTrackParticleNameIfInFile(InDetKeys.xAODPixelTrackParticleContainer()),
                              "",
                              InDetKeys.xAODPixelTrackParticleContainer(),
                              topSequence)


if InDetFlags.doTrackSegmentsDisappearing() and InDetFlags.doParticleCreation():
    if doCreation :
        createTrackParticles(InDetKeys.DisappearingTracks(),
                             InDetKeys.DisappearingTracksTruth(),
                             InDetKeys.xAODDisappearingTrackParticleContainer(),
                             topSequence,
                             trt_pid_tool=True, pixel_dedx=True) # @TODO TRT PID necessary?


if InDetFlags.doTrackSegmentsSCT() and InDetFlags.doParticleCreation():
    if doCreation :
        createTrackParticles(InDetKeys.SCTTracks(),
                             "",
                             InDetKeys.xAODSCTTrackParticleContainer(),
                             topSequence,
                             trt_pid_tool=False, pixel_dedx=False)


if InDetFlags.doTrackSegmentsTRT() and InDetFlags.doParticleCreation():
    if doCreation :
        if InDetFlags.doCosmics() :
            # to reproduce shared hits of TRTTrackParticles stored in the ESD of q220
            createTrackParticles(InDetKeys.TRTTracks(),
                                 "",
                                 InDetKeys.xAODTRTTrackParticleContainer(),
                                 topSequence,
                                 ("PRDtoTrackMap" + InDetKeys.UnslimmedTracks()),
                                 "TRTTracks")
        else :
            createTrackParticles(InDetKeys.TRTTracks(),
                                 "",
                                 InDetKeys.xAODTRTTrackParticleContainer(),
                                 topSequence,
                                 trt_pid_tool=False, pixel_dedx=False)

if InDetFlags.doStoreTrackSeeds() and InDetFlags.doParticleCreation():
    if doCreation :
        createTrackParticles(InDetKeys.SiSPSeedSegments(),
                             "",
                             InDetKeys.SiSPSeedSegments()+"TrackParticle",
                             topSequence,
                             trt_pid_tool=False, pixel_dedx=False)
