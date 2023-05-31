# #######################################################################
# NOTICE
#
# This software (or technical data) was produced for the U.S. Government
# under contract, and is subject to the Rights in Data-General Clause
# 52.227-14, Alt. IV (DEC 2007).
#
# Copyright 2022 The MITRE Corporation. All Rights Reserved.
# ####################################################################### 

cmake_minimum_required(VERSION 3.1)
project(BIQTIris)
set(CMAKE_CXX_STANDARD 11)

OPTION(BUILD_SHARED_LIBS "Builds shared libraries for certain dependencies. Recommended: ON" ON)
OPTION(BUILD_STATIC_LIBS "Builds static libraries for certain dependencies. Recommended: OFF" OFF)

if(NOT WIN32)
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	set(CMAKE_CXX_FLAGS "-g -fPIC")	
endif()

if(DEFINED ENV{BIQT_HOME})
	if(BIQT_HOME)
		message(WARNING "Both -DBIQT_HOME and a BIQT_HOME environment variable are defined. Using -DBIQT_HOME (${BIQT_HOME}).")
	else()
    file(TO_CMAKE_PATH $ENV{BIQT_HOME} BIQT_HOME)
  endif()
endif()

if(NOT BIQT_HOME)
	message(FATAL_ERROR "BIQT_HOME must be defined as an environment variable or passed to CMake as a parameter using -DBIQT_HOME.")
else()
  file(TO_CMAKE_PATH ${BIQT_HOME} BIQT_HOME)
endif()

if(CMAKE_INSTALL_PREFIX)
	# Ignore CMAKE_INSTALL_PREFIX, but let the user know about it.
	message(WARNING "Overriding provided CMAKE_INSTALL_PREFIX with ${BIQT_HOME}/providers.")
endif()

set(CMAKE_INSTALL_PREFIX "${BIQT_HOME}/providers")

find_package(OpenCV REQUIRED core highgui)
find_library(JSONCPP NAMES libjsoncpp jsoncpp HINTS ${BIQT_HOME}/lib ${BIQT_HOME}/lib64 ${BIQT_HOME}/bin ${BIQT_HOME}/../../lib ${BIQT_HOME}/../../lib64 ${BIQT_HOME}/../../bin)
find_library(BIQTAPI NAMES libbiqtapi biqtapi HINTS ${BIQT_HOME}/lib ${BIQT_HOME}/lib64 ${BIQT_HOME}/bin ${BIQT_HOME}/../../lib ${BIQT_HOME}/../../lib64 ${BIQT_HOME}/../../bin)

if(NOT JSONCPP)
	message(FATAL_ERROR "Failed to find library: jsoncpp.")
endif()

if(NOT BIQTAPI)
	message(FATAL_ERROR "Failed to find library: biqtapi.")
endif()

# BUILD THE IRIS LIBRARY FILE #################################################

include_directories("include" ${OpenCV_INCLUDE_DIRS} ${BIQT_HOME}/include ${BIQT_HOME}/../../include)

file(GLOB BIQTIris_SRC "src/*.cpp")

add_library(BIQTIris SHARED ${BIQTIris_SRC})

target_link_libraries(BIQTIris ${JSONCPP} ${BIQTAPI} ${OpenCV_LIBS} ${CMAKE_DL_LIBS})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/descriptor.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# INSTALL THE IRIS LIBRARY FILE ###############################################
if (WIN32)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CopyWindowsLibraries.cmake)  
	install(TARGETS   BIQTIris           DESTINATION "${CMAKE_PROJECT_NAME}")
	install(FILES     descriptor.json    DESTINATION "${CMAKE_PROJECT_NAME}")
  install(FILES     ${DEPENDENCY_DLLS} DESTINATION "${CMAKE_PROJECT_NAME}")
else ()
	install(TARGETS   BIQTIris           DESTINATION "${CMAKE_PROJECT_NAME}")
	install(FILES     descriptor.json    DESTINATION "${CMAKE_PROJECT_NAME}")
endif()

set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "RPM")
set(CPACK_PACKAGE_NAME "biqt-iris")
set(CPACK_PACKAGE_VERSION_MAJOR 0)
set(CPACK_PACKAGE_VERSION_MINOR 4)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_VERSION 0.4)
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "MITRE Corporation")
set(CPACK_PACKAGE_VENDOR "MITRE Corporation")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_RPM_PACKAGE_REQUIRES "biqt >= 0.1 jsoncpp >= 1.8 opencv >= 2.4")
include(CPack)

