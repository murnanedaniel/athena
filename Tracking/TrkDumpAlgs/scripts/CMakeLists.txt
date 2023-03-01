cmake_minimum_required(VERSION 3.14)

# Project's name
message(STATUS "compiling converter ROOT 2 CSV")
project (Dumper
         VERSION 0.1
         LANGUAGES CXX)

# Set the output folder where your program will be created
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#cmake_path(GET ${CMAKE_CURRENT_SOURCE_DIR} PARENT_PATH PROJECT_DIR)
get_filename_component(PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)

# minimal dependency versions. they are defined here in a single place so
# they can be easily upgraded, although they might not be used if the
# dependency is included via `add_subdirectory(...)`.
set (boost_version 1.71)
set (root_version 6.20)

# required packages
# NOTE we use boost::filesystem instead of std::filesystem since the later
#   is not uniformly supported even on compilers that nominally support C++17 

find_package (Boost ${boost_version} REQUIRED CONFIG COMPONENTS filesystem program_options unit_test_framework graph)
include_directories (${Boost_INCLUDE_DIRS})

find_package(ROOT ${root_version} REQUIRED CONFIG COMPONENTS Geom)
find_package(ROOT ${root_version} REQUIRED CONFIG COMPONENTS Core Geom GenVector Hist Tree TreePlayer)
message(STATUS "include dir: ${PROJECT_DIR}")
include_directories(${PROJECT_DIR}/src)
include_directories(${PROJECT_DIR}/TrkDumpAlgs)
include(${ROOT_DIR}/RootMacros.cmake)

# =====================
# Graph Creator Library
# =====================

set(SOURCE_FILES {$PROJECT_DIR}/src/ROOT2CSVconverter.cpp)
add_library(
Dumper SHARED
${PROJECT_DIR}/TrkDumpAlgs/ROOT2CSVconverter.hpp
${PROJECT_DIR}/src/ROOT2CSVconverter.cpp
)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard used for compiling")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# removes leading or trailing whitespace
#string (STRIP ${SOURCE_FILES} SOURCE_FILES)
# defines how to build the executable from source files
#add_executable(ROOT2CSVconverter.exe ${PROJECT_DIR}/src/ROOT2CSVconverter.cpp)
add_executable(ROOT2CSVconverter.exe ${PROJECT_DIR}/src/ROOT2CSVconverter.cpp)

target_compile_features (Dumper INTERFACE cxx_std_17)

# ---Rootcint dictionary ---------------------------------------------------
# rootcint -f dic.cpp -c LinkDef.hpp
root_generate_dictionary(ROOT2CSVconverter.hpp MODULE Dumper LINKDEF ${PROJECT_DIR}/TrkDumpAlgs/LinkDef.hpp)
target_link_libraries(ROOT2CSVconverter.exe Dumper ${Boost_LIBRARIES} ${ROOT_LIBRARIES})# ${ROOT_EXE_LINKER_FLAGS})
target_include_directories(Dumper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
