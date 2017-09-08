include_directories(..)

set(HEADERS-ELLPACK
  bml_add_ellpack.h
  bml_adjungate_triangle_ellpack.h
  bml_allocate_ellpack.h
  bml_convert_ellpack.h
  bml_copy_ellpack.h
  bml_diagonalize_ellpack.h
  bml_export_ellpack.h
  bml_getters_ellpack.h
  bml_import_ellpack.h
  bml_introspection_ellpack.h
  bml_inverse_ellpack.h
  bml_multiply_ellpack.h
  bml_normalize_ellpack.h
  bml_norm_ellpack.h
  bml_parallel_ellpack.h
  bml_scale_ellpack.h
  bml_setters_ellpack.h
  bml_submatrix_ellpack.h
  bml_threshold_ellpack.h
  bml_trace_ellpack.h
  bml_transpose_ellpack.h
  bml_types_ellpack.h
  bml_utilities_ellpack.h)

set(SOURCES-ELLPACK
  bml_add_ellpack.c
  bml_adjungate_triangle_ellpack.c
  bml_allocate_ellpack.c
  bml_convert_ellpack.c
  bml_copy_ellpack.c
  bml_diagonalize_ellpack.c
  bml_export_ellpack.c
  bml_getters_ellpack.c
  bml_import_ellpack.c
  bml_introspection_ellpack.c
  bml_inverse_ellpack.c
  bml_multiply_ellpack.c
  bml_normalize_ellpack.c
  bml_norm_ellpack.c
  bml_parallel_ellpack.c
  bml_scale_ellpack.c
  bml_setters_ellpack.c
  bml_submatrix_ellpack.c
  bml_threshold_ellpack.c
  bml_trace_ellpack.c
  bml_transpose_ellpack.c
  bml_utilities_ellpack.c)

add_library(bml-ellpack OBJECT ${SOURCES-ELLPACK})
set_target_properties(bml-ellpack
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-ellpack
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

set(SOURCES-ELLPACK-TYPED
  bml_add_ellpack_typed.c
  bml_adjungate_triangle_ellpack_typed.c
  bml_allocate_ellpack_typed.c
  bml_convert_ellpack_typed.c
  bml_copy_ellpack_typed.c
  bml_diagonalize_ellpack_typed.c
  bml_export_ellpack_typed.c
  bml_getters_ellpack_typed.c
  bml_import_ellpack_typed.c
  bml_introspection_ellpack_typed.c
  bml_inverse_ellpack_typed.c
  bml_multiply_ellpack_typed.c
  bml_normalize_ellpack_typed.c
  bml_norm_ellpack_typed.c
  bml_parallel_ellpack_typed.c
  bml_scale_ellpack_typed.c
  bml_setters_ellpack_typed.c
  bml_submatrix_ellpack_typed.c
  bml_threshold_ellpack_typed.c
  bml_trace_ellpack_typed.c
  bml_transpose_ellpack_typed.c
  bml_utilities_ellpack_typed.c)

include(${CMAKE_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-ellpack single_real "${SOURCES-ELLPACK-TYPED}")
bml_add_typed_library(bml-ellpack double_real "${SOURCES-ELLPACK-TYPED}")
bml_add_typed_library(bml-ellpack single_complex "${SOURCES-ELLPACK-TYPED}")
bml_add_typed_library(bml-ellpack double_complex "${SOURCES-ELLPACK-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-ellpack-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellpack-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellpack-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellpack-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
