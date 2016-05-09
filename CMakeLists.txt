cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_C_COMPILER gcc)

project(couple C)

option(USE_MPI "Use MPI!" OFF)
if(USE_MPI)
    find_package(MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})       ## WARNING!!!! THIS NEEDS TO BE CHECKED!!!
    add_definitions(-D_MPI)
endif(USE_MPI)

# Include Directory, needed for main/main.c
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/structs")

# Specify binary and library output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

add_subdirectory(initio)
add_subdirectory(solver)
add_subdirectory(structs)
add_subdirectory(main)
