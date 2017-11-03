# General source files (always included)
list(APPEND solver_srcs
  # Insert new files below here, alphabetically
  tools.c
  solver.c
  )
list(REMOVE_DUPLICATES solver_srcs)

if(USE_MPI)
#    find_package(MPI)
    list(APPEND lib_depends ${MPI_C_LIBRARIES}  messaging)
endif(USE_MPI)

add_library("solver" STATIC EXCLUDE_FROM_ALL ${solver_srcs})
target_link_libraries(solver ${lib_depends}  )
