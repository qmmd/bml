include_directories(${CMAKE_SOURCE_DIR}/src/C-interface)

set(SOURCES_TYPED
  add_matrix_typed.c
  adjacency_matrix_typed.c
  adjungate_triangle_matrix_typed.c
  allocate_matrix_typed.c
  copy_matrix_typed.c
  diagonalize_matrix_typed.c
  get_element_typed.c
  get_set_diagonal_typed.c
  get_sparsity_typed.c
  import_export_matrix_typed.c
  inverse_matrix_typed.c
  multiply_matrix_typed.c
  normalize_matrix_typed.c
  norm_matrix_typed.c
  scale_matrix_typed.c
  set_row_typed.c
  submatrix_matrix_typed.c
  test_bml_gemm_typed.c
  threshold_matrix_typed.c
  trace_matrix_typed.c
  transpose_matrix_typed.c)

include(${CMAKE_SOURCE_DIR}/cmake/bmlAddTypedLibrary.cmake)
bml_add_typed_library(bmltests single_real "${SOURCES_TYPED}")
bml_add_typed_library(bmltests double_real "${SOURCES_TYPED}")
bml_add_typed_library(bmltests single_complex "${SOURCES_TYPED}")
bml_add_typed_library(bmltests double_complex "${SOURCES_TYPED}")

add_library(bmltests
  $<TARGET_OBJECTS:bmltests-single_real>
  $<TARGET_OBJECTS:bmltests-double_real>
  $<TARGET_OBJECTS:bmltests-single_complex>
  $<TARGET_OBJECTS:bmltests-double_complex>)
set_target_properties(bmltests
  PROPERTIES
  POSITION_INDEPENDENT_CODE yes)

add_executable(bml-test
  add_matrix.c
  adjacency_matrix.c
  adjungate_triangle_matrix.c
  allocate_matrix.c
  bml_test.c
  copy_matrix.c
  diagonalize_matrix.c
  get_element.c
  get_set_diagonal.c
  get_sparsity.c
  import_export_matrix.c
  inverse_matrix.c
  multiply_matrix.c
  normalize_matrix.c
  norm_matrix.c
  scale_matrix.c
  set_row.c
  submatrix_matrix.c
  test_bml_gemm.c
  threshold_matrix.c
  trace_matrix.c
  transpose_matrix.c)

target_link_libraries(bml-test bmltests bml ${LINK_LIBRARIES})
set_target_properties(bml-test
  PROPERTIES
  LINK_FLAGS "--coverage")
if(OPENMP_FOUND)
  set_target_properties(bml-test
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    LINK_FLAGS ${OpenMP_C_FLAGS})
endif()
if(MPI_C_FOUND AND MPI_C_COMPILE_FLAGS)
  set_target_properties(bml-test
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_C_LINK_FLAGS})
endif()

foreach(N
    add
    adjacency
    adjungate_triangle
    allocate
    bml_gemm
    copy
    diagonalize
    get_element
    get_set_diagonal
    get_sparsity
    import_export
    inverse
    multiply
    norm
    normalize
    scale
    set_row
    submatrix
    threshold
    trace
    transpose)
  foreach(T dense ellpack ellsort)
    foreach(P single_real double_real single_complex double_complex)
      #if(BML_MPI)
        #add_test(${N}-${T}-${P} mpirun -np 1 bml-test -n ${N} -t ${T} -p ${P})
      #else()
        add_test(${N}-${T}-${P} bml-test -n ${N} -t ${T} -p ${P})
      #endif()
      if(NOT BML_MPI AND NOT BML_OPENMP AND VALGRIND)
        add_test(${N}-${T}-${P}-valgrind ${VALGRIND} --error-exitcode=1
          --leak-check=full --read-var-info=yes --track-origins=yes
          ${CMAKE_CURRENT_BINARY_DIR}/bml-test -n ${N} -t ${T} -p ${P})
      endif()
    endforeach()
  endforeach()
endforeach()
