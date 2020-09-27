cmake_minimum_required(VERSION 3.15)

option(BUILD_EXE "Switch whether build exe file or dll file" OFF)

if(WIN32)
    cmake_policy(SET CMP0091 NEW)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif ()

set(PROJECT unPoseSymmetricMaster)
project(${PROJECT})
# version
set(serial "1.0")
# basically, we use Release
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

########################
##### Dependencies #####
########################
###
# eigen3
###
find_package(Eigen3 3.3 REQUIRED)

###
# cpprestsdk (for json)
###
find_package(cpprestsdk CONFIG REQUIRED)

##############################
##### Dependencies (end) #####
##############################

##########################################
##### Include path / cpp files       #####
##########################################
# Include directories for this project
set(INCLUDE_DIR
    src
    ${EIGEN3_INCLUDE_DIR}
    submodule/FromZ/src
    submodule/libigl/include
)

if(BUILD_EXE)
    # Target for this build
    set(CPP_FILES
        main.cpp
        UnPoseSymmetricMaster.cpp
    )
else()
    # Create shared library
    set(CPP_FILES
        UnPoseSymmetricMaster.cpp
    )
endif()
##########################################
##### Include path / cpp files (end) #####
##########################################


###########################
##### Libraries       #####
###########################
if(WIN32)
    set(LINK_LIBS
        Eigen3::Eigen
        cpprestsdk::cpprest
        cpprestsdk::cpprestsdk_zlib_internal
        cpprestsdk::cpprestsdk_brotli_internal
    )
elseif (APPLE)
    set(LINK_LIBS
        Eigen3::Eigen
        cpprestsdk::cpprest
        cpprestsdk::cpprestsdk_zlib_internal
        cpprestsdk::cpprestsdk_boost_internal
        cpprestsdk::cpprestsdk_openssl_internal
    )
endif ()

###########################
##### Libraries (end) #####
###########################

if(BUILD_EXE)
    # Target for this build
    add_executable(${PROJECT} ${CPP_FILES})
else()
    # Create shared library
    add_library(${PROJECT} SHARED ${CPP_FILES})
endif()

# Add a include files
target_include_directories(${PROJECT} PUBLIC ${INCLUDE_DIR})
# Add external libraries
target_link_libraries(${PROJECT} PRIVATE ${LINK_LIBS})

if(WIN32)
set(
    MY_RELEASE_OPTIONS
    /O2
    /MT
)
set(
    MY_DEBUG_OPTIONS
    -bigobj
)
elseif (APPLE)
set(
    MY_RELEASE_OPTIONS
    -O2
)
set(
    MY_DEBUG_OPTIONS
    -g
)
endif ()

target_compile_options(${PROJECT} PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")
target_compile_options(${PROJECT} PUBLIC "$<$<CONFIG:RELEASE>:${MY_RELEASE_OPTIONS}>")
target_compile_features(${PROJECT} PUBLIC cxx_std_17)
target_compile_definitions(${PROJECT} PUBLIC _TURN_OFF_PLATFORM_STRING)
target_compile_definitions(${PROJECT} PUBLIC _USE_MATH_DEFINES)
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
endif ()
