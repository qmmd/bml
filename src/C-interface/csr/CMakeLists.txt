set(HEADERS-CSR
  bml_types_csr.h
  bml_allocate_csr.h
  bml_setters_csr.h
  bml_getters_csr.h
  bml_copy_csr.h
  bml_trace_csr.h
  bml_scale_csr.h
  bml_utilities_csr.h
  bml_norm_csr.h
  bml_convert_csr.h
  bml_export_csr.h
  bml_import_csr.h
  bml_introspection_csr.h
  bml_transpose_csr.h
  bml_threshold_csr.h
  bml_diagonalize_csr.h
  bml_inverse_csr.h
  bml_add_csr.h
  bml_normalize_csr.h
  bml_multiply_csr.h
)


set(SOURCES-CSR
  bml_allocate_csr.c
  bml_setters_csr.c
  bml_getters_csr.c
  bml_copy_csr.c
  bml_trace_csr.c
  bml_scale_csr.c
  bml_utilities_csr.c
  bml_norm_csr.c
  bml_convert_csr.c
  bml_export_csr.c
  bml_import_csr.c
  bml_introspection_csr.c
  bml_transpose_csr.c  
  bml_threshold_csr.c
  bml_diagonalize_csr.c
  bml_inverse_csr.c
  bml_add_csr.c
  bml_normalize_csr.c
  bml_multiply_csr.c  
)

add_library(bml-csr OBJECT ${SOURCES-CSR})
set_target_properties(bml-csr
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-csr
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

set(SOURCES-CSR-TYPED
  bml_allocate_csr_typed.c
  bml_setters_csr_typed.c
  bml_getters_csr_typed.c
  bml_copy_csr_typed.c
  bml_trace_csr_typed.c
  bml_scale_csr_typed.c
  bml_utilities_csr_typed.c
  bml_norm_csr_typed.c
  bml_convert_csr_typed.c
  bml_export_csr_typed.c
  bml_import_csr_typed.c
  bml_introspection_csr_typed.c
  bml_transpose_csr_typed.c  
  bml_threshold_csr_typed.c
  bml_diagonalize_csr_typed.c
  bml_inverse_csr_typed.c
  bml_add_csr_typed.c
  bml_normalize_csr_typed.c
  bml_multiply_csr_typed.c
  )
  
include(${CMAKE_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-csr single_real "${SOURCES-CSR-TYPED}")
bml_add_typed_library(bml-csr double_real "${SOURCES-CSR-TYPED}")
bml_add_typed_library(bml-csr single_complex "${SOURCES-CSR-TYPED}")
bml_add_typed_library(bml-csr double_complex "${SOURCES-CSR-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-csr-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-csr-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-csr-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-csr-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
