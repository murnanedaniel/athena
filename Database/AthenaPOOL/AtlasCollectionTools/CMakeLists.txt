################################################################################
# Package: AtlasCollectionTools
################################################################################

# Declare the package name:
atlas_subdir( AtlasCollectionTools )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PRIVATE
                          DataQuality/GoodRunsLists
                          Database/APR/CollectionBase
                          Database/APR/CollectionUtilities
                          Database/APR/POOLCore
                          Database/APR/PersistencySvc
                          Database/APR/StorageSvc
                          TestPolicy )

# Component(s) in the package:
atlas_add_library( AtlasCollectionTools
                   src/LumiMetaHandler.cxx
                   src/GlobalUniqueKey.cxx
                   src/GlobalSum.cxx
                   NO_PUBLIC_HEADERS
                   LINK_LIBRARIES GoodRunsListsLib
                   PRIVATE_LINK_LIBRARIES CollectionBase CollectionUtilities POOLCore PersistencySvc StorageSvc )

atlas_add_executable( CollAppend
                      src/CollAppend.cxx
                      LINK_LIBRARIES GoodRunsListsLib CollectionBase CollectionUtilities POOLCore PersistencySvc StorageSvc AtlasCollectionTools )

atlas_add_executable( CollAppendXML
                      src/CollAppendXML.cxx
                      LINK_LIBRARIES GoodRunsListsLib CollectionBase CollectionUtilities POOLCore PersistencySvc StorageSvc AtlasCollectionTools )

# Install files from the package:
atlas_install_python_modules( python/countGuidsClient.py python/eventLookupClient.py )
atlas_install_scripts( python/listDatasets.py )
atlas_install_scripts( python/runEventLookup.py python/runGuidsCount.py python/tagExtract.py python/tagDSList.py )

