if(NOT EMBED_libfmftsaxs)
  cmake_minimum_required(VERSION 3.0)
endif()

set(libfmftsaxs_version_major 1)
set(libfmftsaxs_version_minor 0)
set(libfmftsaxs_version_patch 0)
set(libfmftsaxs_version
  "${libfmftsaxs_version_major}.${libfmftsaxs_version_minor}.${libfmftsaxs_version_patch}")

project(fmftsaxs
  LANGUAGES C
  VERSION ${libfmftsaxs_version})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

#################
# BUILD OPTIONS #
#################
option( SONAME             "Set the (SO)VERSION of the target"	   ON  )
option( BUILD_SHARED_LIBS  "Build Shared Library (OFF for Static)" OFF )
option( PROFILE            "Generate profiling information"        OFF )
option( BUILD_TESTS        "Build tests"                           ON  )
option( BUILD_DOCS         "Create and install the HTML based 
                            API documentation (requires Doxygen)"  ON )
option( USE_LTO            "Use Link Time Optimization"            ON  )
option( MPI                "Message passing interface"             OFF )
option( VERBOSE            "Verbosity"                             ON  )

# Default installation directory is set to $HOME
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH "Base directory for installation" FORCE)
endif()

set (CMAKE_PREFIX_PATH $ENV{HOME})

##################
# COMPILER FLAGS #
##################
if (BUILD_SHARED_LIBS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden -fPIC")
endif ()

if(PROFILE)
  set(CMAKE_C_FLAGS "-pg ${CMAKE_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if (MPI)
  add_definitions(-D_MPI_)
  set(CMAKE_C_COMPILER mpicc)
endif ()

if (VERBOSE)
  add_definitions(-D_SXS_VERBOSE_)
endif ()

set(CMAKE_C_FLAGS
  "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_C_FLAGS_DEBUG
  "${CMAKE_C_FLAGS_DEBUG} -Wall -Wshadow -Wpointer-arith -Wcast-qual -Winline -Werror -Wextra -Wfatal-errors -Wstrict-prototypes -g -pg")
set(CMAKE_C_FLAGS_RELEASE
  "${CMAKE_C_FLAGS_RELEASE}")

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND USE_LTO)
  include(CheckCCompilerFlag)
  check_c_compiler_flag(-flto HAS_LINK_TIME_OPTIMIZATION)
  #check whether ar supports plugins
  set(AR_SUPPORT 1)
  if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    execute_process(
      COMMAND ar --help
      COMMAND grep plugin
      OUTPUT_VARIABLE AR_SUPPORT)
  endif()
  if(HAS_LINK_TIME_OPTIMIZATION AND AR_SUPPORT)
    set(CMAKE_C_FLAGS
      "${CMAKE_C_FLAGS} -flto")
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
      set(CMAKE_AR "gcc-ar")
      set(CMAKE_RANLIB "gcc-ranlib")
    endif()
  endif()
endif()

############
# SOURCES  #
############
file(GLOB LIBLBFGSB_SOURCES 
	lbfgsb/src/lbfgsb.c
	lbfgsb/src/linesearch.c 
	lbfgsb/src/subalgorithms.c 
	lbfgsb/src/print.c
    lbfgsb/src/linpack.c
    lbfgsb/src/miniCBLAS.c
    lbfgsb/src/timer.c)
file(GLOB LIBLBFGSB_PUBLIC_HEADERS lbfgsb/src/lbfgsb.h)

file(GLOB libfmftsaxs_SOURCES src/*.c)
file(GLOB libfmftsaxs_PUBLIC_HEADERS src/*.h )

#################
# FIND PACKAGES #
#################
find_package(mol2 REQUIRED)

if (MPI)
  find_package(MPI REQUIRED)
endif ()

############
# LIBRARY  #
############
#set install subdirs
set(LIB_INSTALL_DIR     "lib"                        CACHE PATH "Where to install libraries to.")
set(INCLUDE_INSTALL_DIR "include"                    CACHE PATH "Where to install headers to.")
set(CONFIG_INSTALL_DIR  "lib/cmake/${PROJECT_NAME}"  CACHE PATH "Where to install configs to.")

add_library(lbfgsb OBJECT ${LIBLBFGSB_SOURCES})

add_library(fmftsaxs ${libfmftsaxs_SOURCES} $<TARGET_OBJECTS:lbfgsb>)
target_link_libraries(fmftsaxs fftw3 mol2 gmp check rt m)
target_include_directories(fmftsaxs
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
    "$<BUILD_INTERFACE:${CMAKE_HOME_DIRECTORY}/lbfgsb/src>"
    "$<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>")
    
set_target_properties(fmftsaxs PROPERTIES PUBLIC_HEADER "${libfmftsaxs_PUBLIC_HEADERS};${LIBLBFGSB_PUBLIC_HEADERS}")
set_target_properties(fmftsaxs PROPERTIES VERSION "${libfmftsaxs_version}")

if (MPI)
  target_include_directories(fmftsaxs 
				PUBLIC 
				"$<BUILD_INTERFACE:${MPI_INCLUDE_PATH}>")
  if(MPI_COMPILE_FLAGS)
    set_target_properties(fmftsaxs
				PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(fmftsaxs 
				PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif ()

add_subdirectory(tools)

############
# INSTALL  #
############
install(TARGETS               fmftsaxs
	EXPORT                    fmftsaxsTargets
	LIBRARY DESTINATION       "${LIB_INSTALL_DIR}"
	ARCHIVE DESTINATION       "${LIB_INSTALL_DIR}"
	PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}/fmftsaxs")

############
#  CONFIG  #
############
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/fmftsaxs/fmftsaxsConfigVersion.cmake
                                 VERSION       ${libfmftsaxs_version}
                                 COMPATIBILITY AnyNewerVersion)

configure_package_config_file(cmake/fmftsaxsConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/fmftsaxs/fmftsaxsConfig.cmake
                              INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}")

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/fmftsaxs/fmftsaxsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/fmftsaxs/fmftsaxsConfigVersion.cmake
  DESTINATION ${CONFIG_INSTALL_DIR})

install(EXPORT fmftsaxsTargets
  FILE fmftsaxsTargets.cmake
  DESTINATION ${CONFIG_INSTALL_DIR})

#################
# DOCUMENTATION #
#################
# add a target to generate API documentation with Doxygen
if(BUILD_DOCS)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
    add_custom_target(doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
      COMMENT "Generating API documentation with Doxygen" 
      VERBATIM)
    #set install dirs
    set(DOC_INSTALL_DIR	"doc" CACHE PATH "Where to install docs to.")
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc/" 
  	    DESTINATION "${DOC_INSTALL_DIR}/fmftsaxs")
  else()
    message("Could not find Doxygen.")
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCS)

###########
# TESTING #
###########
if(BUILD_TESTS)
  find_package(Check REQUIRED)
  if(CHECK_FOUND)
    enable_testing()
    add_subdirectory(tests)
  else()
    message("Could not find Check for testing")
  endif()
endif()
