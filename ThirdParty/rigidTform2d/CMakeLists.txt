######################################################################
# CMakeLists.txt generated for component estimateAffineRigid2D
# Product type: static library
######################################################################
cmake_minimum_required(VERSION 3.12)
project(estimateAffineRigid2D)
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

######################################################################
## Target definition and commands
######################################################################
add_library(estimateAffineRigid2D SHARED ${START_DIR}/rt_nonfinite.c
    ${START_DIR}/rtGetNaN.c
    ${START_DIR}/rtGetInf.c
    ${START_DIR}/estimateAffineRigid2D.c
    ${START_DIR}/estimateAffineRigid2D_emxutil.c
    ${START_DIR}/estimateAffineRigid2D_emxAPI.c)
set_target_properties(estimateAffineRigid2D PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE ON
    RUNTIME_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    LIBRARY_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>"
    ARCHIVE_OUTPUT_DIRECTORY "${BINARY_START_DIR}/$<0:>")
target_compile_features(estimateAffineRigid2D PUBLIC c_std_99)
target_compile_definitions(estimateAffineRigid2D PRIVATE -DMODEL=estimateAffineRigid2D)
target_include_directories(estimateAffineRigid2D PUBLIC $<BUILD_INTERFACE:${START_DIR}/codegen/lib/estimateAffineRigid2D>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/codegen/lib/estimateAffineRigid2D>
    $<BUILD_INTERFACE:${START_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>>
    $<BUILD_INTERFACE:${MATLAB_ROOT}/extern/include>)
target_link_libraries(estimateAffineRigid2D INTERFACE $<$<BOOL:${FOUND_LIBM}>:m>)


######################################################################
## Build success message
######################################################################
add_custom_command(TARGET estimateAffineRigid2D POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "\\#\\#\\# Created library: $<TARGET_FILE:estimateAffineRigid2D>")

######################################################################
## Call toolchain hook function if defined
######################################################################
if(COMMAND toolchain_target_hook)
    toolchain_target_hook(estimateAffineRigid2D)
endif()
