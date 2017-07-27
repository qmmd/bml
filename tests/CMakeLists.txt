include_directories(
  ${CMAKE_SOURCE_DIR}/src/C-interface
  ${CMAKE_BINARY_DIR}/src/Fortran-interface)

set(SOURCES_TYPED
  add_matrix_typed.c
  adjacency_matrix_typed.c
  adjungate_triangle_matrix_typed.c
  allocate_matrix_typed.c
  convert_matrix_typed.c
  copy_matrix_typed.c
  diagonalize_matrix_typed.c
  get_element_typed.c
  get_sparsity_typed.c
  inverse_matrix_typed.c
  multiply_matrix_typed.c
  normalize_matrix_typed.c
  norm_matrix_typed.c
  scale_matrix_typed.c
  get_set_diagonal_typed.c
  set_row_typed.c
  submatrix_matrix_typed.c
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
  convert_matrix.c
  copy_matrix.c
  diagonalize_matrix.c
  get_element.c
  get_sparsity.c
  inverse_matrix.c
  multiply_matrix.c
  normalize_matrix.c
  norm_matrix.c
  scale_matrix.c
  get_set_diagonal.c
  set_row.c
  submatrix_matrix.c
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

add_executable(bml-test-fortran
  test_m.F90
  allocate_matrix.F90
  bml_test.F90)
target_link_libraries(bml-test-fortran bml ${LINK_LIBRARIES})
set_target_properties(bml-test-fortran
  PROPERTIES
  LINK_FLAGS "--coverage")
if(OPENMP_FOUND)
  set_target_properties(bml-test-fortran
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()
if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(bml-test-fortran
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
endif()

foreach(N add adjacency adjungate_triangle allocate convert copy diagonalize
    get_element get_sparsity inverse multiply norm  normalize scale
    get_set_diagonal set_row submatrix threshold trace transpose)
  foreach(T dense ellpack ellsort)
    foreach(P single_real double_real single_complex double_complex)
      #if(BML_MPI)
        #add_test(${N}-${T}-${P} mpirun -np 1 bml-test -n ${N} -t ${T} -p ${P})
      #else()
        add_test(${N}-${T}-${P} bml-test -n ${N} -t ${T} -p ${P})
        add_test(${N}-${T}-${P}-Fortran bml-test-fortran ${N} ${T} ${P})
      #endif()
      if(NOT BML_MPI AND NOT BML_OPENMP AND VALGRIND)
        add_test(${N}-${T}-${P}-valgrind ${VALGRIND} --error-exitcode=1
          --leak-check=full --read-var-info=yes --track-origins=yes
          ${CMAKE_CURRENT_BINARY_DIR}/bml-test -n ${N} -t ${T} -p ${P})
      endif()
    endforeach()
  endforeach()
endforeach()
