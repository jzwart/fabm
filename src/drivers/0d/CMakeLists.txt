cmake_minimum_required(VERSION 2.8.8)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(fabm0d Fortran)

# Find GOTM directory.
# This directory will also be configurable in the CMake GUI or interactive mode.
find_path(GOTM_BASE
          NAMES src/gotm/gotm.F90
          HINTS ENV GOTMDIR
          PATH_SUFFIXES GOTM/gotm-git gotm-git
          DOC "Directory with GOTM source code (www.gotm.net)")

# Make NetCDF usage a configurable option.
option(FABM_USE_NetCDF "Enable NetCDF output" ON)
if(FABM_USE_NetCDF)
  add_definitions(-DREAL_4B=real\(4\))
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${GOTM_BASE}/src/cmake/Modules/")
  find_package(NetCDF REQUIRED)
  if (NetCDF_STATIC_MSVC_BUILD)
    add_compile_options("/libs:static")
  endif()
  set (OUTPUT_MANAGER ${GOTM_BASE}/src/output/output_manager_core.F90 ${GOTM_BASE}/src/output/output_manager.F90 ${GOTM_BASE}/src/output/text_output.F90 ${GOTM_BASE}/src/output/netcdf_output.F90)
  set (FIELD_MANAGER ${GOTM_BASE}/src/util/field_manager.F90)
endif(FABM_USE_NetCDF)

# Specify additional compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  endif()
endif()

# Force FABM_HOST to 0d and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST 0d)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm)

# Build 0d executable.
add_executable (fabm0d
                main.F90
                fabm0d.F90
                shared.F90
                output.F90
                fabm_0d.h
                ${GOTM_BASE}/src/airsea/solar_zenith_angle.F90
                ${GOTM_BASE}/src/airsea/albedo_water.F90
                ${GOTM_BASE}/src/airsea/short_wave_radiation.F90
                ${GOTM_BASE}/src/input/input.F90
                ${FIELD_MANAGER}
                ${OUTPUT_MANAGER}
                ${GOTM_BASE}/src/util/eqstate.F90
                ${GOTM_BASE}/src/util/gridinterpol.F90
                ${GOTM_BASE}/src/util/time.F90
                )

# Group GOTM source files in a folder in IDEs that support this.
source_group(gotm REGULAR_EXPRESSION ${GOTM_BASE})

# Link with FABM library.
target_link_libraries(fabm0d fabm)

if(FABM_USE_NetCDF)
  target_link_libraries(fabm0d ${NetCDF_LIBRARIES})
  set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES ${NetCDF_INCLUDE_DIRS})
  set_property(TARGET fabm0d APPEND PROPERTY COMPILE_DEFINITIONS NETCDF4)
  set_property(SOURCE ${OUTPUT_MANAGER} APPEND PROPERTY COMPILE_DEFINITIONS NETCDF_FMT)
  if (NetCDF_STATIC_MSVC_BUILD)
    set_property(TARGET fabm0d PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmt\"")
  endif()
endif(FABM_USE_NetCDF)

# Set module directory (output of *.mod)
set_property(TARGET fabm0d PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Set include directories with modules (*.mod) - our own and FABM's.
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/fabm/modules")

# Set include directories with header files (*.h) for GOTM - version.h, cppdefs.h, ode_solver_template.F90.
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTM_BASE}/include")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTM_BASE}/src/util")

# Set include directories with header files (*.h) - our own and FABM's.
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/../../../include")

# Set default installation prefix.
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/fabm/0d")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/fabm/0d")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/fabm/0d")
endif()

# Set cache variable with installation prefix.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # None set - set our default prefix.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM's 0d driver in" FORCE)
else()
  # User specified a prefix - just update the variable description.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM's 0d driver in")
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

install(TARGETS fabm0d
  DESTINATION bin)
