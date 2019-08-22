if(MPI_C_INCLUDE_PATH)
  include_directories(${MPI_C_INCLUDE_PATH})
endif()

# Public headers.
set(HEADERS-C
  bml.h
  bml_add.h
  bml_adjungate_triangle.h
  bml_allocate.h
  bml_convert.h
  bml_copy.h
  bml_diagonalize.h
  bml_export.h
  bml_getters.h
  bml_import.h
  bml_init.h
  bml_introspection.h
  bml_inverse.h
  bml_logger.h
  bml_multiply.h
  bml_norm.h
  bml_normalize.h
  bml_parallel.h
  bml_scale.h
  bml_setters.h
  bml_shutdown.h
  bml_submatrix.h
  bml_threshold.h
  bml_trace.h
  bml_transpose.h
  bml_transpose_triangle.h
  bml_types.h
  bml_utilities.h)
install(FILES ${HEADERS-C} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Private headers.
set(HEADERS-C-PRIVATE
  blas.h
  typed.h)

set(SOURCES-C
  bml_add.c
  bml_adjungate_triangle.c
  bml_allocate.c
  bml_convert.c
  bml_copy.c
  bml_diagonalize.c
  bml_export.c
  bml_getters.c
  bml_import.c
  bml_init.c
  bml_introspection.c
  bml_inverse.c
  bml_logger.c
  bml_multiply.c
  bml_norm.c
  bml_normalize.c
  bml_parallel.c
  bml_scale.c
  bml_setters.c
  bml_shutdown.c
  bml_submatrix.c
  bml_threshold.c
  bml_trace.c
  bml_transpose.c
  bml_transpose_triangle.c
  bml_utilities.c)
add_library(bml-c OBJECT ${SOURCES-C})
set_target_properties(bml-c
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-c
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
if(MPI_C_FOUND AND MPI_C_COMPILE_FLAGS)
  set_target_properties(bml-c
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

add_subdirectory(dense)
add_subdirectory(ellpack)
add_subdirectory(ellblock)
add_subdirectory(ellsort)
