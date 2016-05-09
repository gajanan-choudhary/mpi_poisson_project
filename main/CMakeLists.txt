##############################################################################
# Build executable
list(APPEND main_sources
    main.c
    )
# Add main as build target
add_executable(main ${main_sources})

list(APPEND lib_depends structs initio solver)

if(USE_MPI)
    find_package(MPI)
    list(APPEND lib_depends ${MPI_C_LIBRARIES})
endif(USE_MPI)

# Specify library dependencies
target_link_libraries(main ${lib_depends}  )
