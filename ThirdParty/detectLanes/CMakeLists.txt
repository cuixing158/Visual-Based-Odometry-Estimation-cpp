######################################################################
# CMakeLists.txt generated for component detectLaneMarkerRidge
# Product type: static library
######################################################################
cmake_minimum_required(VERSION 3.12)
# project(detectLaneMarkerRidge)
if(DEFINED ENV{CMAKE_EXPORT_COMPILE_COMMANDS} AND NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_EXPORT_COMPILE_COMMANDS $ENV{CMAKE_EXPORT_COMPILE_COMMANDS})
endif()


######################################################################
## Path variables
######################################################################
get_filename_component(START_DIR ./ ABSOLUTE)
# set(MATLAB_ROOT /opt_disk2/rd22946/MATLAB/R2023a CACHE PATH "")
if(${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(BUILDING_IN_SOURCE_DIR TRUE)
    set(BINARY_START_DIR "${START_DIR}")
else()
    set(BUILDING_IN_SOURCE_DIR FALSE)
    set(BINARY_START_DIR "${CMAKE_BINARY_DIR}")
endif()
if(APPLE)
    if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL x86_64)
        list(APPEND MATLAB_ROOT_SYSLIB_PATHS ${MATLAB_ROOT}/bin/maci64)
        list(APPEND MATLAB_ROOT_SYSLIB_PATHS ${MATLAB_ROOT}/sys/os/maci64)
        list(APPEND MATLAB_ROOT_SYSINCLUDE_PATHS ${MATLAB_ROOT}/toolbox/eml/externalDependency/omp/maci64/include)
    elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL arm64)
        list(APPEND MATLAB_ROOT_SYSLIB_PATHS ${MATLAB_ROOT}/bin/maca64)
        list(APPEND MATLAB_ROOT_SYSINCLUDE_PATHS ${MATLAB_ROOT}/toolbox/eml/externalDependency/omp/maca64/include)
    endif()
elseif(UNIX AND ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL x86_64) AND ("${CMAKE_SYSTEM_NAME}" STREQUAL Linux))
    list(APPEND MATLAB_ROOT_SYSLIB_PATHS ${MATLAB_ROOT}/bin/glnxa64)
elseif(WIN32 AND ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL AMD64))
    list(APPEND MATLAB_ROOT_SYSLIB_PATHS ${MATLAB_ROOT}/bin/win64
        ${MATLAB_ROOT}/lib/win64)
endif()

if(DEFINED MATLAB_ROOT_SYSLIB_PATHS)
    list(APPEND CMAKE_LIBRARY_PATH ${MATLAB_ROOT_SYSLIB_PATHS})
    link_directories(${MATLAB_ROOT_SYSLIB_PATHS})
endif()
if(DEFINED MATLAB_ROOT_SYSINCLUDE_PATHS)
    list(APPEND CMAKE_INCLUDE_PATH ${MATLAB_ROOT_SYSINCLUDE_PATHS})
endif()

if(WIN32 AND NOT MINGW)
    set(SYSLIB_PREFIX "lib")
else()
    set(SYSLIB_PREFIX "")
endif()
find_package(OpenMP REQUIRED)

######################################################################
## System libraries
######################################################################
find_library(FOUND_LIBM m NO_SYSTEM_ENVIRONMENT_PATH PATHS ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
    ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
find_library(FOUND_LIBSTDCXX stdc++ NO_SYSTEM_ENVIRONMENT_PATH PATHS ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
    ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})

######################################################################
## Target definition and commands
######################################################################
add_library(detectLaneMarkerRidge SHARED ${START_DIR}/rt_nonfinite.cpp
    ${START_DIR}/rtGetNaN.cpp
    ${START_DIR}/rtGetInf.cpp
    ${START_DIR}/detectLaneMarkerRidge.cpp)
set_target_properties(detectLaneMarkerRidge PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE ON
    RUNTIME_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    LIBRARY_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    ARCHIVE_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>")
target_compile_features(detectLaneMarkerRidge PUBLIC cxx_std_11)
target_compile_definitions(detectLaneMarkerRidge PRIVATE -DMODEL=detectLaneMarkerRidge)
target_include_directories(detectLaneMarkerRidge PUBLIC $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${MATLAB_ROOT}/extern/include>)
target_link_libraries(detectLaneMarkerRidge INTERFACE $<$<BOOL:${FOUND_LIBM}>:m>
    $<$<BOOL:${FOUND_LIBSTDCXX}>:stdc++>
    OpenMP::OpenMP_C
    OpenMP::OpenMP_CXX)


######################################################################
## Build success message
######################################################################
add_custom_command(TARGET detectLaneMarkerRidge POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "\\#\\#\\# Created library: $<TARGET_FILE:detectLaneMarkerRidge>")

######################################################################
## Call toolchain hook function if defined
######################################################################
if(COMMAND toolchain_target_hook)
    toolchain_target_hook(detectLaneMarkerRidge)
endif()
