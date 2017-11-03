# General source files (always included)
list(APPEND messaging_srcs
  # Insert new files below here, alphabetically
  message_error.c
  message_barrier.c
  )

list(REMOVE_DUPLICATES messaging_srcs)

list(APPEND lib_depends ${MPI_C_LIBRARIES})

add_library("messaging" STATIC EXCLUDE_FROM_ALL ${messaging_srcs})
target_link_libraries(messaging ${lib_depends})
