set(HEADERS-DISTRIBUTED
  bml_allocate_distributed2d.h
  bml_export_distributed2d.h
  bml_import_distributed2d.h
  bml_multiply_distributed2d.h )

set(SOURCES-DISTRIBUTED
  bml_allocate_distributed2d.c
  bml_export_distributed2d.c
  bml_import_distributed2d.c
  bml_multiply_distributed2d.c )

add_library(bml-distributed2d OBJECT ${SOURCES-DISTRIBUTED})
set_target_properties(bml-distributed2d
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)
if(OPENMP_FOUND)
  set_target_properties(bml-distributed2d
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

set(SOURCES-DISTRIBUTED-TYPED
  bml_export_distributed2d_typed.c
  bml_import_distributed2d_typed.c
  bml_multiply_distributed2d_typed.c
  )

include(${PROJECT_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bml-distributed2d single_real "${SOURCES-DISTRIBUTED-TYPED}")
bml_add_typed_library(bml-distributed2d double_real "${SOURCES-DISTRIBUTED-TYPED}")
bml_add_typed_library(bml-distributed2d single_complex "${SOURCES-DISTRIBUTED-TYPED}")
bml_add_typed_library(bml-distributed2d double_complex "${SOURCES-DISTRIBUTED-TYPED}")
if(OPENMP_FOUND)
  set_target_properties(bml-distributed2d-single_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-distributed2d-double_real
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-distributed2d-single_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
  set_target_properties(bml-distributed2d-double_complex
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()
