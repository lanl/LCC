cmake_minimum_required(VERSION 3.10.0)
project(lcc VERSION 1.0.0 LANGUAGES Fortran)

#include(GNUInstallDirs)

set(dir ${CMAKE_CURRENT_SOURCE_DIR}/build/)
set(CMAKE_BUILD_DIRECTORY ${dir})
#set(CMAKE_INSTALL_BINDIR ${dir})

include(FindPkgConfig)

find_package(PROGRESS CONFIG QUIET)
pkg_check_modules(PROGRESS REQUIRED progress)
message(STATUS "Found progress: ${PROGRESS_LDFLAGS}")
list(APPEND LINK_LIBRARIES ${PROGRESS_LDFLAGS})

find_package(BML CONFIG QUIET)
pkg_check_modules(BML REQUIRED bml)
list(APPEND LINK_LIBRARIES BML::bml)
list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})
message(STATUS "Found bml: ${BML_LDFLAGS}")

if(PROGRESS_MPI)
  message(STATUS "Will build with MPI")
  add_definitions(-DDO_MPI)
endif()

message(STATUS "Extra FC Flags ${EXTRA_FCFLAGS}")

if(DEFINED EXTRA_FCFLAGS)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FCFLAGS}")
endif()

list(APPEND LINK_LIBRARIES ${EXTRA_LINK_FLAGS})

message(STATUS "Project sources = " ${PROJECT_SOURCE_DIR} )
include_directories(${PROJECT_SOURCE_DIR}/)
include_directories(${CMAKE_BINARY_DIR}/)
include_directories(${BML_INCLUDEDIR})
include_directories(${PROGRESS_INCLUDEDIR})

function(progress_appendix myappendix main_and_srcs)
list(GET main_and_srcs 0 main)
include_directories(${PROGRESS_INCLUDEDIR})
add_executable(${myappendix} ${main})
target_sources(${myappendix} PRIVATE ${ARGN})
target_link_libraries(${myappendix} PUBLIC
    ${LINK_LIBRARIES})
  set_target_properties(${myappendix}
    PROPERTIES
    LINK_FLAGS "")
endfunction(progress_appendix)

progress_appendix(lcc_main lcc_main.F90
                            lcc_allocation_mod.F90
                            lcc_check_mod.F90 
                            lcc_message_mod.F90
                            lcc_parser_mod.F90
                            lcc_constants_mod.F90
                            lcc_structs_mod.F90 
                            lcc_lattice_mod.F90 
                            lcc_aux_mod.F90 
                            lcc_build_mod.F90 
                            lcc_regular_mod.F90 
                            lcc_message_mod.F90
                            lcc_radial_mod.F90
                            lcc_string_mod.F90
                            lcc_template_mod.F90
                            )

install(TARGETS lcc_main DESTINATION ${CMAKE_INSTALL_BINDIR})
