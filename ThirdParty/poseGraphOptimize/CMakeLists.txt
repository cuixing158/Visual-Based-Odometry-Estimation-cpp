######################################################################
# CMakeLists.txt generated for component poseGraphOptimize
# Product type: static library
######################################################################
# cmake_minimum_required(VERSION 3.12)
# project(poseGraphOptimize)
if(DEFINED ENV{CMAKE_EXPORT_COMPILE_COMMANDS} AND NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_EXPORT_COMPILE_COMMANDS $ENV{CMAKE_EXPORT_COMPILE_COMMANDS})
endif()


######################################################################
## Path variables
######################################################################
get_filename_component(START_DIR . ABSOLUTE)
if(${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(BUILDING_IN_SOURCE_DIR TRUE)
    set(BINARY_START_DIR "${START_DIR}")
else()
    set(BUILDING_IN_SOURCE_DIR FALSE)
    set(BINARY_START_DIR "${CMAKE_BINARY_DIR}")
endif()
if(WIN32 AND NOT MINGW)
    set(SYSLIB_PREFIX "lib")
else()
    set(SYSLIB_PREFIX "")
endif()

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
add_library(poseGraphOptimize SHARED ${START_DIR}/coder_posix_time.c
    ${START_DIR}/CXSparse/Source/cs_ri_bundle.cpp
    ${START_DIR}/CXSparse/Source/cs_ci_bundle.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/solve_from_lu.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/solve_from_qr.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/makeCXSparseMatrix.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/unpackCXStruct.cpp
    ${START_DIR}/rt_nonfinite.cpp
    ${START_DIR}/rtGetNaN.cpp
    ${START_DIR}/rtGetInf.cpp
    ${START_DIR}/poseGraphOptimize.cpp)
set_target_properties(poseGraphOptimize PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE ON
    RUNTIME_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    LIBRARY_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    ARCHIVE_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>")
target_compile_features(poseGraphOptimize PUBLIC cxx_std_11)
target_compile_definitions(poseGraphOptimize PRIVATE -DMODEL=poseGraphOptimize)
target_include_directories(poseGraphOptimize PUBLIC $<BUILD_INTERFACE:${START_DIR}/poseGraphOptimize>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/poseGraphOptimize>
    $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/Include>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/Include>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/SuiteSparse_config>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/SuiteSparse_config>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/CXSparseSupport>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/CXSparseSupport>)
target_link_libraries(poseGraphOptimize INTERFACE $<$<BOOL:${FOUND_LIBM}>:m>
    $<$<BOOL:${FOUND_LIBSTDCXX}>:stdc++>)


######################################################################
## Build success message
######################################################################
add_custom_command(TARGET poseGraphOptimize POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "\\#\\#\\# Created library: $<TARGET_FILE:poseGraphOptimize>")

######################################################################
## Call toolchain hook function if defined
######################################################################
if(COMMAND toolchain_target_hook)
    toolchain_target_hook(poseGraphOptimize)
endif()
