#--------------------CMake file for ProtonDetector in Geant 4.9.6.p02--------------
#Author: D. Perez-Loureiro 12/02/2013
#
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ProtonDetector)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Cmake Path for installation
#
#CMAKE_INSTALL_PREFIX 
#----------------------------------------------------------------------------
# Load some basic macros which are needed later on
include(FindROOT.cmake)

#---------------------------------------------------------------------------
# Create the directory for the ROOT files
#
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files")


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
# Setup ROOT include directories

include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
#Setup ROOT linking directories
#
link_directories( ${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
#Setup of compiler flags
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
   EXEC_PROGRAM("sw_vers -productVersion | cut -d . -f 1-2" OUTPUT_VARIABLE MAC_OS_VERSION)
   MESSAGE("--- Found a Mac OS X System ${MAC_OS_VERSION}")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

      STRING(COMPARE EQUAL "10.5" "${MAC_OS_VERSION}" MAC_OS_10_5)
      IF(MAC_OS_10_5)
        SET(CMAKE_CXX_FLAGS "-m64")
        SET(CMAKE_Fortran_FLAGS "-m64")
      ENDIF(MAC_OS_10_5)

      SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
      SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")

      
      
      endif (CMAKE_COMPILER_IS_GNUCXX)

endif (CMAKE_SYSTEM_NAME MATCHES Darwin) 

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Set everything needed for the root dictonary and create the
# dictionary
set(ProtonDetector_HEADERS
		     ${PROJECT_SOURCE_DIR}/include/ProtonDetectorTrack.hh 
		     ${PROJECT_SOURCE_DIR}/include/ProtonDetectorSimpleTrack.hh
		     ${PROJECT_SOURCE_DIR}/include/ProtonDetectorBeamInfo.hh
		     ${PROJECT_SOURCE_DIR}/include/ProtonDetectorData.hh
		     )
set(ProtonDetector_LINKDEF ${PROJECT_SOURCE_DIR}/src/ProtonDetectorLinkDef.hh )
set(ProtonDetector_DICTIONARY ProtonDetectorDict.cc) 
ROOT_GENERATE_DICTIONARY("${ProtonDetector_HEADERS}" "${ProtonDetector_LINKDEF}" "${ProtonDetector_DICTIONARY}" "${ROOT_INCLUDE_DIR}")


#----------------------------------------------------------------------------
# Set the library version in the main CMakeLists.txt
SET(ProtonDetector_MAJOR_VERSION 0)
SET(ProtonDetector_MINOR_VERSION 1)
SET(ProtonDetector_PATCH_VERSION 0)
SET(ProtonDetector_VERSION "${ProtonDetector_MAJOR_VERSION}.${ProtonDetector_MINOR_VERSION}.${ProtonDetector_PATCH_VERSION}")
SET(ProtonDetector_LIBRARY_PROPERTIES ${ProtonDetector_LIBRARY_PROPERTIES}
    VERSION "${ProtonDetector_VERSION}"
    SOVERSION "${ProtonDetector_MAJOR_VERSION}"
    SUFFIX ".sl"
)

 



add_library(Proton SHARED ProtonDetector.cc ${ProtonDetector_DICTIONARY} ${sources} ${headers})
target_link_libraries(Proton ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
set_target_properties(Proton PROPERTIES ${ProtonDetector_LIBRARY_PROPERTIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ProtonDetector ProtonDetector.cc ${ProtonDetector_DICTIONARY} ${sources} ${headers})
target_link_libraries(ProtonDetector ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_dependencies(ProtonDetector Proton)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build the executable. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(PROTONDETECTOR_SCRIPTS
vis.mac
Input.mac
Zpos_al23.root
)

foreach(_script ${PROTONDETECTOR_SCRIPTS})
 configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(ProtonDetector DEPENDS ProtonDetector)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS actarsim DESTINATION bin)
#install(TARGETS actarsim DESTINATION lib)

#install(FILES ${headers}
#        DESTINATION "./include"
#        DESTINATION "${PROJECT_BINARY_DIR}/include"
#         )
