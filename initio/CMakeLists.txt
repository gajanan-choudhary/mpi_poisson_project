# General source files (always included)
list(APPEND initio_srcs
  # Insert new files below here, alphabetically
  read_mesh.c
  )
list(REMOVE_DUPLICATES initio_srcs)

list(APPEND lib_depends structs)

add_library("initio" STATIC EXCLUDE_FROM_ALL ${initio_srcs})

target_link_libraries(initio ${lib_depends})
