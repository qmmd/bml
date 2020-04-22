set(HEADERS-DENSE
  bml_add_dense.h
  bml_adjungate_triangle_dense.h
  bml_allocate_dense.h
  bml_convert_dense.h
  bml_copy_dense.h
  bml_diagonalize_dense.h
  bml_export_dense.h
  bml_getters_dense.h
  bml_import_dense.h
  bml_introspection_dense.h
  bml_inverse_dense.h
  bml_multiply_dense.h
  bml_normalize_dense.h
  bml_norm_dense.h
  bml_parallel_dense.h
  bml_scale_dense.h
  bml_setters_dense.h
  bml_threshold_dense.h
  bml_trace_dense.h
  bml_transpose_dense.h
  bml_transpose_triangle_dense.h
  bml_types_dense.h
  bml_utilities_dense.h)

set(SOURCES-DENSE
  bml_add_dense.c
  bml_adjungate_triangle_dense.c
  bml_allocate_dense.c
  bml_convert_dense.c
  bml_copy_dense.c
  bml_diagonalize_dense.c
  bml_export_dense.c
  bml_getters_dense.c
  bml_import_dense.c
  bml_introspection_dense.c
  bml_inverse_dense.c
  bml_multiply_dense.c
  bml_normalize_dense.c
  bml_norm_dense.c
  bml_parallel_dense.c
  bml_scale_dense.c
  bml_setters_dense.c
  bml_threshold_dense.c
  bml_trace_dense.c
  bml_transpose_dense.c
  bml_transpose_triangle_dense.c
  bml_utilities_dense.c)

add_library(bml-dense OBJECT ${SOURCES-DENSE})
set_target_properties(bml-dense
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-dense
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

set(SOURCES-DENSE-TYPED
  bml_add_dense_typed.c
  bml_adjungate_triangle_dense_typed.c
  bml_allocate_dense_typed.c
  bml_convert_dense_typed.c
  bml_copy_dense_typed.c
  bml_export_dense_typed.c
  bml_getters_dense_typed.c
  bml_import_dense_typed.c
  bml_introspection_dense_typed.c
  bml_inverse_dense_typed.c
  bml_multiply_dense_typed.c
  bml_normalize_dense_typed.c
  bml_norm_dense_typed.c
  bml_parallel_dense_typed.c
  bml_scale_dense_typed.c
  bml_setters_dense_typed.c
  bml_threshold_dense_typed.c
  bml_trace_dense_typed.c
  bml_transpose_dense_typed.c
  bml_transpose_triangle_dense_typed.c
  bml_utilities_dense_typed.c)

include(${PROJECT_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-dense single_real "${SOURCES-DENSE-TYPED}")
bml_add_typed_library(bml-dense double_real "${SOURCES-DENSE-TYPED}")
bml_add_typed_library(bml-dense single_complex "${SOURCES-DENSE-TYPED}")
bml_add_typed_library(bml-dense double_complex "${SOURCES-DENSE-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-dense-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-dense-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-dense-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-dense-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
