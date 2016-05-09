# General source files (always included)
list(APPEND structs_srcs
  # Insert new files below here, alphabetically
  element.c
#  node.c
  model.c
  )

list(REMOVE_DUPLICATES structs_srcs)

list(APPEND lib_depends solver)

add_library("structs" STATIC EXCLUDE_FROM_ALL ${structs_srcs})
target_link_libraries(structs ${lib_depends})
