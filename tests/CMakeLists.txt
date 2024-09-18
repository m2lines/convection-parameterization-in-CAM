cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
#policy CMP0076 - target_sources source files are relative to file where target_sources is run
cmake_policy (SET CMP0076 NEW)

set(PROJECT_NAME YOGTest)

project(${PROJECT_NAME} LANGUAGES Fortran)

if(NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none")
endif()


set(NN_Module_DIR "../NN_module/")
set(CAM_Interface_DIR "../CAM_interface/")



#Add cmake directory to the environment module variable
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find the NetCDF installations and set the relevant variables for compilation
# Then link to executables
# Requires more legwork as NetCDF not provided by default
find_package(PkgConfig)
pkg_search_module(NETCDF_FORTRAN netcdf-fortran)
if (NETCDF_FORTRAN_FOUND)
    set(NETCDF_LIBRARIES "${NETCDF_FORTRAN_LDFLAGS}")
    set(NETCDF_INCLUDES "${NETCDF_FORTRAN_INCLUDE_DIRS}")
else()
    set(NETCDF_F90 "YES")
    find_package(NetCDF REQUIRED)
endif()
pkg_search_module(NETCDF_C netcdf)
if (NETCDF_C_FOUND)
    list(APPEND NETCDF_LIBRARIES "${NETCDF_C_LDFLAGS}")
    list(APPEND NETCDF_INCLUDES "${NETCDF_C_INCLUDE_DIRS}")
endif()

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../)

# Generate the main executable
add_executable ( test_NN_module 
  test_NN_module/test.f90
  test_utils.f90
  ${NN_Module_DIR}/precision.f90
  ${NN_Module_DIR}/SAM_consts.f90
  ${NN_Module_DIR}/nn_cf_net.f90
  ${NN_Module_DIR}/nn_convection_flux.f90
)

target_link_libraries( test_NN_module PRIVATE ${NETCDF_LIBRARIES} )
target_include_directories( test_NN_module PRIVATE ${NETCDF_INCLUDES} )

# Generate the main executable
add_executable ( test_CAM_interface 
  test_CAM_interface/test_cam_interface.f90
  test_utils.f90
  ${NN_Module_DIR}/precision.f90
  ${NN_Module_DIR}/SAM_consts.f90
  ${NN_Module_DIR}/nn_cf_net.f90
  ${NN_Module_DIR}/nn_convection_flux.f90
  ${CAM_Interface_DIR}/nn_interface_CAM.f90
)

target_link_libraries( test_CAM_interface PRIVATE ${NETCDF_LIBRARIES} )
target_include_directories( test_CAM_interface PRIVATE ${NETCDF_INCLUDES} )
