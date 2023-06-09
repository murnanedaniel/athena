################################################################################
# Package: InDetTrigRawDataProvider
################################################################################

# Declare the package name:
atlas_subdir( InDetTrigRawDataProvider )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaBaseComps
                          GaudiKernel
                          InnerDetector/InDetRawEvent/InDetRawData
                          InnerDetector/InDetTrigRecTools/InDetTrigToolInterfaces
                          PRIVATE
                          Control/AthenaKernel
                          Control/StoreGate
                          DetectorDescription/IRegionSelector
                          Event/ByteStreamCnvSvcBase
                          InnerDetector/InDetDetDescr/InDetIdentifier
                          InnerDetector/InDetDetDescr/PixelCabling
                          InnerDetector/InDetDetDescr/SCT_Cabling
                          InnerDetector/InDetDetDescr/TRT_Cabling
                          InnerDetector/InDetEventCnv/PixelRawDataByteStreamCnv
                          InnerDetector/InDetEventCnv/SCT_RawDataByteStreamCnv
                          InnerDetector/InDetEventCnv/TRT_RawDataByteStreamCnv
                          Trigger/TrigEvent/TrigSteeringEvent )

# Component(s) in the package:
atlas_add_component( InDetTrigRawDataProvider
                     src/*.cxx
                     src/components/*.cxx
                     LINK_LIBRARIES AthenaBaseComps GaudiKernel InDetRawData AthenaKernel StoreGateLib SGtests IRegionSelector ByteStreamCnvSvcBaseLib InDetIdentifier SCT_CablingLib TrigSteeringEvent )

# Install files from the package:
atlas_install_headers( InDetTrigRawDataProvider )

