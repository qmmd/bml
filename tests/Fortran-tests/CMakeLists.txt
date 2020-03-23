include_directories(${BML_INCLUDEDIR})
include_directories(${CMAKE_BINARY_DIR}/src/Fortran-interface)
include(${CMAKE_SOURCE_DIR}/cmake/bmlGetTypes.cmake)

# Adding the low-level source for each test
set(FORTRAN-SOURCES-TYPED
  add_matrix_typed.F90
  adjungate_triangle_typed.F90
  allocate_matrix_typed.F90
  convert_matrix_typed.F90
  copy_matrix_typed.F90
  get_bandwidth_typed.F90
  get_diagonal_typed.F90
  get_row_typed.F90
  diagonalize_matrix_typed.F90
  get_element_typed.F90
  io_matrix_typed.F90
  inverse_matrix_typed.F90
  multiply_matrix_typed.F90
  normalize_matrix_typed.F90
  scale_matrix_typed.F90
  set_diagonal_typed.F90
  set_element_typed.F90
  set_row_typed.F90
  threshold_matrix_typed.F90
  trace_matrix_typed.F90
  transpose_matrix_typed.F90
  )

# Setting all the relevant element/numeric types
if(BML_COMPLEX)
  set(TYPES
    single_real
    double_real
    single_complex
    double_complex
    )
else()
  set(TYPES
    single_real
    double_real
    )
endif()


# Preprocessing the low-level sources
foreach(T ${TYPES})
  foreach(S ${FORTRAN-SOURCES-TYPED})
    string(REPLACE "typed" ${T} S_TYPED ${S})
    get_filename_component(S_TYPED ${S_TYPED} NAME)

    if(${T} MATCHES "single_real")
      set(FORMATED_KIND "real")
      set(FORMATED_PREC "sp")
      set(NUMBER_PREC "")
    elseif(${T} MATCHES "double_real")
      set(FORMATED_KIND "real")
      set(FORMATED_PREC "dp")
      set(NUMBER_PREC "_dp")
    elseif(${T} MATCHES "single_complex")
      set(FORMATED_KIND "complex")
      set(FORMATED_PREC "sp")
      set(NUMBER_PREC "")
    elseif(${T} MATCHES "double_complex")
      set(FORMATED_KIND "complex")
      set(FORMATED_PREC "dp")
      set(NUMBER_PREC "_dp")
    endif()

    add_custom_command(OUTPUT ${S_TYPED}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${S}
      COMMAND ${CMAKE_BINARY_DIR}/convert-template
      -DDUMMY_KIND=${FORMATED_KIND}
      -DDUMMY_PREC=${FORMATED_PREC}
      -D_MP=${NUMBER_PREC}
      -D_typed="_"${T}
      -o ${S_TYPED}
      ${CMAKE_CURRENT_SOURCE_DIR}/${S})
    list(APPEND FORTRAN-SOURCES-PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${S_TYPED})
  endforeach()
endforeach()

# Setting link flags
#set(FORTRAN_TEST_LINK_FLAGS
#-L${CMAKE_SOURCE_DIR}/build/src -lbml_fortran -lbml
#${OpenMP_Fortran_FLAGS}
#${BLAS_LIBRARIES}
#${LINK_LIBRARIES}
#)

# Building the general executable
add_executable(bml-testf
  testf.F90
  add_matrix.F90
  adjungate_triangle.F90
  allocate_matrix.F90
  convert_matrix.F90
  copy_matrix.F90
  diagonalize_matrix.F90
  get_bandwidth.F90
  get_diagonal.F90
  get_element.F90
  get_row.F90
  inverse_matrix.F90
  io_matrix.F90
  multiply_matrix.F90
  normalize_matrix.F90
  scale_matrix.F90
  set_diagonal.F90
  set_element.F90
  set_row.F90
  threshold_matrix.F90
  trace_matrix.F90
  transpose_matrix.F90
  prec.F90
  ${FORTRAN-SOURCES-PUBLIC}
  ${FORTRAN-SOURCES-PRIVATE}
  )

target_link_libraries(bml-testf bml_fortran bml ${LINK_LIBRARIES})

# More flags
if(OPENMP_FOUND)
  set_target_properties(bml-testf
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS} )
endif()

# Function that adds the test
function(test_formats ${formats})
  set(language fortran)
  foreach(T ${formats} )
    set(precisions single_real double_real)
    if(BML_COMPLEX)
      list(APPEND precisions single_complex double_complex)
    endif()
    foreach(P ${precisions})
      set(testname ${language}-${N}-${T}-${P})
      add_test(${testname} ${CMAKE_CURRENT_BINARY_DIR}/bml-testf -n ${N} -t ${T} -p ${P})
      if(NOT BML_MPI AND NOT BML_OPENMP AND VALGRIND AND BML_VALGRIND)
        add_test(${testname}-valgrind ${VALGRIND} ${VALGRIND_COMMON_ARGS}
          ${CMAKE_CURRENT_BINARY_DIR}/bml-testf -n ${N} -t ${T} -p ${P})
        # [nicolasbock] Ignore errors in valgrind tests for now.
        # Please revert this change once all memory issues have been
        # fixed.
        set_tests_properties(${testname}-valgrind PROPERTIES
          SKIP_RETURN_CODE 1)
      endif()
    endforeach()
  endforeach()
endfunction(test_formats)

# Adding tests
foreach(N
    add
    #adjungate
    allocate
    multiply
    convert
    copy
    #diagonalize
    get_bandwidth
    get_diagonal
    #get_element
    get_row
    inverse
    io
    multiply
    #normalize
    #scale
    set_diagonal
    set_row
    threshold
    trace
    transpose
    )
  set(formats dense ellpack ellsort ellblock)
  # Setting specific formats for a particular bml function
  if(${N} MATCHES get_bandwidth)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES inverse)
    set(formats dense ellpack ellblock)
  elseif(${N} MATCHES threshold)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES multiply)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES get_diagonal)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES get_row)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES io)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES set_diagonal)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES set_row)
    set(formats dense ellpack ellsort)
  elseif(${N} MATCHES transpose)
    set(formats dense ellpack ellsort)
  endif()
  test_formats(${formats})
endforeach()
