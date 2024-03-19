cmake_minimum_required(VERSION 3.14)
project(CAM-ML Fortran)

if(NOT CMAKE_Fortran_FLAGS_PROFILE)
set(CMAKE_Fortran_FLAGS_PROFILE
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -pg -fno-omit-frame-pointer"
    CACHE STRING "" FORCE
)
endif()
if(NOT CMAKE_EXE_LINKER_FLAGS_PROFILE)
set(CMAKE_EXE_LINKER_FLAGS_PROFILE
    "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} -pg"
    CACHE STRING "" FORCE
)
endif()
if(NOT CMAKE_SHARED_LINKER_FLAGS_PROFILE)
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE
    "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO} -pg"
    CACHE STRING "" FORCE
)
endif()
if(NOT CMAKE_STATIC_LINKER_FLAGS_PROFILE)
set(CMAKE_STATIC_LINKER_FLAGS_PROFILE
    "${CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO}"
    CACHE STRING "" FORCE
)
endif()
if(NOT CMAKE_MODULE_LINKER_FLAGS_PROFILE)
set(CMAKE_MODULE_LINKER_FLAGS_PROFILE
    "${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO} -pg"
    CACHE STRING "" FORCE
)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcpp -Mfreeform -Wall")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fopenmp -ffree-line-length-none -Wall")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -qopenmp -warn all")
endif()

find_package(netCDF-Fortran REQUIRED)
find_package(netCDF REQUIRED)

# Testing NN
find_package(FTorch)
message(STATUS "Building with the Fortran PyTorch coupling")

add_executable(test_nn test.f90)
target_link_libraries(test_nn PRIVATE FTorch::ftorch netCDF::netcdf netCDF::netcdff)
target_sources(test_nn PRIVATE nn_cf_net.f90 nn_cf_net_torch.f90 nn_convection_flux.f90 SAM_consts.f90)
#target_sources(test_nn PUBLIC nn_cf_net.f90 nn_cf_net_torch.f90 nn_convection_flux.f90 nn_convection_flux_torch.f90)
#target_link_libraries(test_nn PUBLIC FTorch::ftorch netCDF::netcdf netCDF::netcdff)

