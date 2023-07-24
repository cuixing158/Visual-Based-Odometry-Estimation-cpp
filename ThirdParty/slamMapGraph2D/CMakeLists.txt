######################################################################
# CMakeLists.txt generated for component slamMapGraph2D
# Product type: static library
######################################################################
cmake_minimum_required(VERSION 3.12)
project(slamMapGraph2D)
if(DEFINED ENV{CMAKE_EXPORT_COMPILE_COMMANDS} AND NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_EXPORT_COMPILE_COMMANDS $ENV{CMAKE_EXPORT_COMPILE_COMMANDS})
endif()


######################################################################
## Path variables
######################################################################
get_filename_component(START_DIR ./ ABSOLUTE)
set(MATLAB_ROOT /opt/matlab/R2023a CACHE PATH "")
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
add_library(slamMapGraph2D SHARED ${START_DIR}/coder_posix_time.c
    ${START_DIR}/CXSparse/Source/cs_ri_bundle.cpp
    ${START_DIR}/CXSparse/Source/cs_ci_bundle.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/solve_from_lu.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/solve_from_qr.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/makeCXSparseMatrix.cpp
    ${START_DIR}/CXSparse/CXSparseSupport/unpackCXStruct.cpp
    ${START_DIR}/slamMapGraph2D_data.cpp
    ${START_DIR}/rt_nonfinite.cpp
    ${START_DIR}/rtGetNaN.cpp
    ${START_DIR}/rtGetInf.cpp
    ${START_DIR}/optimizePoseGraph.cpp
    ${START_DIR}/tic.cpp
    ${START_DIR}/svd.cpp
    ${START_DIR}/xnrm2.cpp
    ${START_DIR}/xdotc.cpp
    ${START_DIR}/mtimes.cpp
    ${START_DIR}/sparse.cpp
    ${START_DIR}/norm.cpp
    ${START_DIR}/toc.cpp
    ${START_DIR}/sqrt.cpp
    ${START_DIR}/xaxpy.cpp
    ${START_DIR}/xrotg.cpp
    ${START_DIR}/xrot.cpp
    ${START_DIR}/xswap.cpp
    ${START_DIR}/insertionsort.cpp
    ${START_DIR}/introsort.cpp
    ${START_DIR}/heapsort.cpp
    ${START_DIR}/TrustRegionIndefiniteDogLegInterface.cpp
    ${START_DIR}/BlockInserter2.cpp
    ${START_DIR}/TrustRegionIndefiniteDogLegSE2.cpp
    ${START_DIR}/BlockMatrix.cpp
    ${START_DIR}/poseGraph.cpp
    ${START_DIR}/slamPoseGraph.cpp
    ${START_DIR}/sparse1.cpp
    ${START_DIR}/PoseGraphHelpers.cpp
    ${START_DIR}/SEHelpers.cpp
    ${START_DIR}/Sim3Helpers.cpp
    ${START_DIR}/CXSparseAPI.cpp
    ${START_DIR}/PoseGraphOptimizer.cpp
    ${START_DIR}/myGraph.cpp
    ${START_DIR}/slamMapGraph2D_rtwutil.cpp)
set_target_properties(slamMapGraph2D PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE ON
    RUNTIME_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    LIBRARY_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    ARCHIVE_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>")
target_compile_features(slamMapGraph2D PUBLIC cxx_std_11)
target_compile_definitions(slamMapGraph2D PRIVATE -DMODEL=slamMapGraph2D)
target_include_directories(slamMapGraph2D PUBLIC $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/Include>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/Include>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/SuiteSparse_config>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/SuiteSparse_config>
    $<BUILD_INTERFACE:${START_DIR}/CXSparse/CXSparseSupport>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/CXSparse/CXSparseSupport>
    $<BUILD_INTERFACE:${MATLAB_ROOT}/extern/include>)
target_link_libraries(slamMapGraph2D INTERFACE $<$<BOOL:${FOUND_LIBM}>:m>
    $<$<BOOL:${FOUND_LIBSTDCXX}>:stdc++>)


######################################################################
## Build success message
######################################################################
add_custom_command(TARGET slamMapGraph2D POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "\\#\\#\\# Created library: $<TARGET_FILE:slamMapGraph2D>")

######################################################################
## Call toolchain hook function if defined
######################################################################
if(COMMAND toolchain_target_hook)
    toolchain_target_hook(slamMapGraph2D)
endif()
