# Public headers.
set(HEADERS-INTERNAL-BLAS
  bml_gemm.h)
install(FILES ${HEADERS-C} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(SOURCES-INTERNAL-BLAS
  bml_gemm.c)

include(${CMAKE_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-internal-blas single_real "${SOURCES-INTERNAL-BLAS}")
bml_add_typed_library(bml-internal-blas double_real "${SOURCES-INTERNAL-BLAS}")
bml_add_typed_library(bml-internal-blas single_complex "${SOURCES-INTERNAL-BLAS}")
bml_add_typed_library(bml-internal-blas double_complex "${SOURCES-INTERNAL-BLAS}")

if(OPENMP_FOUND)
  set_target_properties(bml-internal-blas-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-internal-blas-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-internal-blas-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-internal-blas-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
