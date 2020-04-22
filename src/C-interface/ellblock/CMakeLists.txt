set(HEADERS-ELLBLOCK
  bml_add_ellblock.h
  bml_allocate_ellblock.h
  bml_convert_ellblock.h
  bml_copy_ellblock.h
  bml_diagonalize_ellblock.h
  bml_export_ellblock.h
  bml_getters_ellblock.h
  bml_import_ellblock.h
  bml_introspection_ellblock.h
  bml_inverse_ellblock.h
  bml_multiply_ellblock.h
  bml_normalize_ellblock.h
  bml_norm_ellblock.h
  bml_scale_ellblock.h
  bml_setters_ellblock.h
  bml_threshold_ellblock.h
  bml_trace_ellblock.h
  bml_transpose_ellblock.h
  bml_types_ellblock.h
  bml_utilities_ellblock.h)

set(SOURCES-ELLBLOCK
  bml_add_ellblock.c
  bml_allocate_ellblock.c
  bml_convert_ellblock.c
  bml_copy_ellblock.c
  bml_diagonalize_ellblock.c
  bml_export_ellblock.c
  bml_getters_ellblock.c
  bml_import_ellblock.c
  bml_introspection_ellblock.c
  bml_inverse_ellblock.c
  bml_multiply_ellblock.c
  bml_normalize_ellblock.c
  bml_norm_ellblock.c
  bml_scale_ellblock.c
  bml_setters_ellblock.c
  bml_threshold_ellblock.c
  bml_trace_ellblock.c
  bml_transpose_ellblock.c
  bml_utilities_ellblock.c
)

add_library(bml-ellblock OBJECT ${SOURCES-ELLBLOCK})
set_target_properties(bml-ellblock
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-ellblock
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

set(SOURCES-ELLBLOCK-TYPED
  bml_add_ellblock_typed.c
  bml_allocate_ellblock_typed.c
  bml_convert_ellblock_typed.c
  bml_copy_ellblock_typed.c
  bml_diagonalize_ellblock_typed.c
  bml_export_ellblock_typed.c
  bml_getters_ellblock_typed.c
  bml_import_ellblock_typed.c
  bml_introspection_ellblock_typed.c
  bml_inverse_ellblock_typed.c
  bml_multiply_ellblock_typed.c
  bml_normalize_ellblock_typed.c
  bml_norm_ellblock_typed.c
  bml_scale_ellblock_typed.c
  bml_setters_ellblock_typed.c
  bml_threshold_ellblock_typed.c
  bml_trace_ellblock_typed.c
  bml_transpose_ellblock_typed.c
  bml_utilities_ellblock_typed.c
)

include(${PROJECT_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-ellblock single_real "${SOURCES-ELLBLOCK-TYPED}")
bml_add_typed_library(bml-ellblock double_real "${SOURCES-ELLBLOCK-TYPED}")
bml_add_typed_library(bml-ellblock single_complex "${SOURCES-ELLBLOCK-TYPED}")
bml_add_typed_library(bml-ellblock double_complex "${SOURCES-ELLBLOCK-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-ellblock-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellblock-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellblock-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-ellblock-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
