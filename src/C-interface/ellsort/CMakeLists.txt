include_directories(..)

set(HEADERS-ELLSORT
  bml_add_ellsort.h
  bml_adjungate_triangle_ellsort.h
  bml_allocate_ellsort.h
  bml_convert_ellsort.h
  bml_copy_ellsort.h
  bml_export_ellsort.h
  bml_getters_ellsort.h
  bml_import_ellsort.h
  bml_introspection_ellsort.h
  bml_multiply_ellsort.h
  bml_normalize_ellsort.h
  bml_norm_ellsort.h
  bml_parallel_ellsort.h
  bml_scale_ellsort.h
  bml_setters_ellsort.h
  bml_submatrix_ellsort.h
  bml_threshold_ellsort.h
  bml_trace_ellsort.h
  bml_transpose_ellsort.h
  bml_types_ellsort.h
  bml_utilities_ellsort.h)

set(SOURCES-ELLSORT
  bml_add_ellsort.c
  bml_adjungate_triangle_ellsort.c
  bml_allocate_ellsort.c
  bml_convert_ellsort.c
  bml_copy_ellsort.c
  bml_export_ellsort.c
  bml_getters_ellsort.c
  bml_import_ellsort.c
  bml_introspection_ellsort.c
  bml_multiply_ellsort.c
  bml_normalize_ellsort.c
  bml_norm_ellsort.c
  bml_parallel_ellsort.c
  bml_scale_ellsort.c
  bml_setters_ellsort.c
  bml_submatrix_ellsort.c
  bml_threshold_ellsort.c
  bml_trace_ellsort.c
  bml_transpose_ellsort.c
  bml_utilities_ellsort.c)

add_library(bml-ellsort OBJECT ${SOURCES-ELLSORT})
set_target_properties(bml-ellsort
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-ellsort
    PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(SOURCES-ELLSORT-TYPED
  bml_add_ellsort_typed.c
  bml_adjungate_triangle_ellsort_typed.c
  bml_allocate_ellsort_typed.c
  bml_convert_ellsort_typed.c
  bml_copy_ellsort_typed.c
  bml_export_ellsort_typed.c
  bml_getters_ellsort_typed.c
  bml_import_ellsort_typed.c
  bml_introspection_ellsort_typed.c
  bml_multiply_ellsort_typed.c
  bml_normalize_ellsort_typed.c
  bml_norm_ellsort_typed.c
  bml_parallel_ellsort_typed.c
  bml_scale_ellsort_typed.c
  bml_setters_ellsort_typed.c
  bml_submatrix_ellsort_typed.c
  bml_threshold_ellsort_typed.c
  bml_trace_ellsort_typed.c
  bml_transpose_ellsort_typed.c
  bml_utilities_ellsort_typed.c)

include(${CMAKE_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-ellsort single_real "${SOURCES-ELLSORT-TYPED}")
bml_add_typed_library(bml-ellsort double_real "${SOURCES-ELLSORT-TYPED}")
bml_add_typed_library(bml-ellsort single_complex "${SOURCES-ELLSORT-TYPED}")
bml_add_typed_library(bml-ellsort double_complex "${SOURCES-ELLSORT-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-ellsort-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellsort-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellsort-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellsort-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
