# General source files (always included)
list(APPEND solver_srcs
  # Insert new files below here, alphabetically
  tools.c
  solver.c
  )
list(REMOVE_DUPLICATES solver_srcs)

add_library("solver" STATIC EXCLUDE_FROM_ALL ${solver_srcs})
