##############################################################################
# Build executable
list(APPEND tests_sources
    tests.c
    )
list(REMOVE_DUPLICATES tests_sources)

list(APPEND lib_depends structs )

if(USE_MPI)
#    find_package(MPI)
    list(APPEND lib_depends ${MPI_C_LIBRARIES} messaging)
endif(USE_MPI)

add_library("tests" STATIC EXCLUDE_FROM_ALL ${tests_sources})
# Specify library dependencies
target_link_libraries(tests ${lib_depends} )
