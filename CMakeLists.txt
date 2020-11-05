cmake_minimum_required(VERSION 3.10)

message(STATUS "CMake version ${CMAKE_VERSION}")

set(LANGUAGES C CXX)
set(BML_BUILD_FORTRAN_INTERFACE TRUE CACHE BOOL
  "Build the Fortran API (requires a Fortran compiler)")
if(BML_BUILD_FORTRAN_INTERFACE)
  list(APPEND LANGUAGES Fortran)
endif()

project(bml ${LANGUAGES}) # don't move this line as it changes PROJECT_SOURCE_DIR

# The library version is versioned off the major version. If the API
# changes, the library version should be bumped.
set(PROJECT_VERSION_MAJOR "1")
set(PROJECT_VERSION_MINOR "3")
set(PROJECT_VERSION_PATCH "1")

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(PROJECT_DESCRIPTION
  "Basic Matrix Library (bml): "
  "A matrix library for linear algebra operations, "
  "supporting several data structures.")
string(REPLACE ";" " " PROJECT_DESCRIPTION "${PROJECT_DESCRIPTION}")

include(GNUInstallDirs)
set(CMAKE_INSTALL_PKG_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "No user specified build type, using default...")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
if(CMAKE_BUILD_TYPE_UPPER STREQUAL DEBUG)
  add_definitions(-DGLOBAL_DEBUG=BML_LOG_DEBUG)
else()
  add_definitions(-DGLOBAL_DEBUG=BML_LOG_INFO)
endif()

set(GNU_C_FLAGS_DEBUG -O0 -g -Wall -Wimplicit -Werror=implicit
  -Wunused-variable -Wno-unknown-pragmas -save-temps -std=c99)
set(GNU_C_FLAGS_RELEASE -O3 -g -std=c99 -DNDEBUG)
set(GNU_C_FLAGS_RELWITHDEBINFO -O3 -g -std=c99 -DNDEBUG)
set(GNU_Fortran_FLAGS_DEBUG -O0 -g --coverage -fcheck=all
  -ffree-line-length-none)
set(GNU_Fortran_FLAGS_RELEASE -O3 -g -ffree-line-length-none)
set(GNU_Fortran_FLAGS_RELWITHDEBINFO -O3 -g -ffree-line-length-none -DNDEBUG)

set(Intel_C_FLAGS_DEBUG -O0 -g -std=c99 -check=conversions,stack,uninit
  -traceback -qopt-report=5 -DINTEL_SDK)
set(Intel_C_FLAGS_RELEASE -O2 -g -std=c99 -qopt-report=5 -DNDEBUG -DINTEL_SDK)
set(Intel_C_FLAGS_RELWITHDEBINFO -O2 -g -std=c99 -qopt-report=5 -DNDEBUG -DINTEL_SDK)
set(Intel_Fortran_FLAGS_DEBUG -O0 -g -check all -assume realloc_lhs
  -traceback -qopt-report=5 -DINTEL_SDK)
set(Intel_Fortran_FLAGS_RELEASE -O2 -g -assume realloc_lhs -qopt-report=5 -DINTEL_SDK)
set(Intel_Fortran_FLAGS_RELWITHDEBINFO -O2 -g -assume realloc_lhs
  -qopt-report=5 -DNDEBUG -DINTEL_SDK)

set(XL_C_FLAGS_DEBUG -O0 -g)
set(XL_C_FLAGS_RELEASE -O3 -g -std=c99 -DNDEBUG)
set(XL_C_FLAGS_RELWITHDEBINFO -O2 -g -std=c99 -DNDEBUG)
set(XL_Fortran_FLAGS_DEBUG -O0 -g)
set(XL_Fortran_FLAGS_RELEASE -O3 -g)
set(XL_Fortran_FLAGS_RELWITHDEBINFO -O3 -g -DNDEBUG)

set(Clang_C_FLAGS_DEBUG -O0 -g -save-temps -std=c99)
set(Clang_C_FLAGS_RELEASE -O2 -g -std=c99 -DNDEBUG)
set(Clang_C_FLAGS_RELWITHDEBINFO -O2 -g -std=c99 -DNDEBUG)

set(MALLOC_ALIGNMENT 64 CACHE INT "Alignment boundary for memory allocations")
add_definitions(-DMALLOC_ALIGNMENT=${MALLOC_ALIGNMENT})

option(DONT_TOUCH_MY_FLAGS "Don't touch compiler flags" OFF)
if(CMAKE_BUILD_TYPE_UPPER AND NOT DONT_TOUCH_MY_FLAGS)
  if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG"
      OR CMAKE_BUILD_TYPE_UPPER STREQUAL "RELEASE"
      OR CMAKE_BUILD_TYPE_UPPER STREQUAL "RELWITHDEBINFO")
    if(CMAKE_C_FLAGS STREQUAL "")
      if(CMAKE_C_COMPILER_ID STREQUAL "GNU"
          OR CMAKE_C_COMPILER_ID STREQUAL "Intel"
          OR CMAKE_C_COMPILER_ID STREQUAL "Clang"
          OR CMAKE_C_COMPILER_ID STREQUAL "XL")
        set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${${CMAKE_C_COMPILER_ID}_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
        string(REPLACE ";" " "
          CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
        if(DEFINED EXTRA_CFLAGS)
          set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
            "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${EXTRA_CFLAGS}")
        endif()
        message(STATUS "Setting C compiler flags to "
          "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
      else()
        message(STATUS "Unknown C compiler ${CMAKE_C_COMPILER_ID}")
      endif()
    else()
      message(STATUS "Using user supplied C compiler flags: ${CMAKE_C_FLAGS}")
      set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
        ${CMAKE_C_FLAGS})
    endif()

    if(CMAKE_Fortran_FLAGS STREQUAL "")
      if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
          OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"
          OR CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
        set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${${CMAKE_Fortran_COMPILER_ID}_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
        string(REPLACE ";" " "
          CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
        if(DEFINED EXTRA_FFLAGS)
          set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
            "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${EXTRA_FFLAGS}")
        endif()
        message(STATUS "Setting Fortran compiler flags to "
          "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
      else()
        message(STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
      endif()
    else()
      message(STATUS "Using user supplied Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
      set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
      ${CMAKE_Fortran_FLAGS})
    endif()
  endif()
endif()

set(BML_LINK_FLAGS "" CACHE STRING "Additional linker flags")
list(APPEND LINK_LIBRARIES ${BML_LINK_FLAGS})

include(FindThreads)
list(APPEND LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

set(BML_OPENMP TRUE CACHE BOOL "Whether to compile with OpenMP support")
if(BML_OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    if(CMAKE_VERSION VERSION_LESS 3.1)
      # There is not Fortran support in FindOpenMP before version 3.1. We will
      # have to hack this.
      message(WARNING "Versions of <cmake-3.1 can not detect Fortran "
        "compiler OpenMP flags. We will assume that the flags found for "
        "the C compiler also work for the Fortran compiler. If this is "
        "not the case, please send email to <nbock@lanl.gov>.")
      set(OpenMP_Fortran_FLAGS ${OpenMP_C_FLAGS})
    endif()
  else()
    message(WARNING "Could not get the compilers to use OpenMP. "
      "Will pretend that this never happened and compile the library "
      "without OpenMP.")
    set(OpenMP_C_FLAGS "")
    set(OpenMP_CXX_FLAGS "")
    set(OpenMP_Fortran_FLAGS "")
  endif()
else()
  message(STATUS "Will not build with OpenMP")
endif()

set(BML_COVERAGE FALSE CACHE BOOL "Whether to generate coverage code")
if(BML_COVERAGE)
  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} --coverage)
endif()

set(BML_MPI FALSE CACHE BOOL "Compile with MPI support")
if(BML_MPI)
  message(STATUS "Will build with MPI")
  find_package(MPI REQUIRED)
  add_definitions(-DDO_MPI)
endif()

option(BML_OMP_OFFLOAD "Compile with OpenMP GPU Offload support" FALSE)
if(BML_OMP_OFFLOAD)
  message(STATUS "Will build with OpenMP GPU Offload support")
  add_definitions(-DUSE_OMP_OFFLOAD)

  find_package(OpenMP 4.5 REQUIRED)

  include(CheckCXXCompilerFlag)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OpenMP_CXX_OFFLOAD_FLAG "-foffload='-lm'")
  elseif(MAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(OpenMP_CXX_OFFLOAD_FLAG "-foffload='-lm'")
  elseif(MAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OpenMP_CXX_OFFLOAD_FLAG "-fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_TOOLKIT_ROOT_DIR}")
  elseif(MAKE_CXX_COMPILER_ID STREQUAL "XL")
    set(OpenMP_CXX_OFFLOAD_FLAG "-qsmp=omp -qoffload")
  endif()
  check_cxx_compiler_flag("${OpenMP_CXX_OFFLOAD_FLAG}" OFFLOAD_CXX_FOUND)
  if(OFFLOAD_CXX_FOUND AND NOT DONT_TOUCH_MY_FLAGS)
    message(STATUS "Found GPU offload flags")
    set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${OpenMP_CXX_OFFLOAD_FLAG}")
    set(OpenMP_C_FLAGS "${OpenMP_C_FLAGS} ${OpenMP_CXX_OFFLOAD_FLAG}")
    set(OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} ${OpenMP_CXX_OFFLOAD_FLAG}")
  endif()

  # list(APPEND LANGUAGES CUDA)

  set(CUDA_FIND_REQUIRED TRUE)
  include(FindCUDA)

  set(GPU_ARCH "sm_60" CACHE STRING "LAMMPS GPU CUDA SM architecture")
  set_property(CACHE GPU_ARCH PROPERTY STRINGS sm_30 sm_50 sm_60 sm_70)

  message("CUDA libraries: ${CUDA_LIBRARIES}")
  message("CUDA CUBLAS libraries: ${CUDA_CUBLAS_LIBRARIES}")
endif()

set(BML_CUDA FALSE CACHE BOOL "Whether to compile with CUDA support")
if(BML_CUDA OR BML_CUSOLVER)
  include(FindCUDA)
endif()
if (CUDA_FOUND)
   message("CUDA CUBLAS libraries: ${CUDA_CUBLAS_LIBRARIES}")
   include_directories(${CUDA_INCLUDE_DIRS})
endif()

set(BML_OPENCL FALSE CACHE BOOL "Whether to compiler with OpenCL support")
if(BML_OPENCL)
  include(FindOpenCL)
endif()

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(BML_MAGMA FALSE CACHE BOOL "Whether to use MAGMA library")
if(BML_MAGMA)
  FIND_PACKAGE(MAGMA)
endif()
if (MAGMA_FOUND)
  add_definitions(-DBML_USE_MAGMA)
  message(STATUS
    "MAGMA was found:\n"
    "  MAGMA_INCLUDE_DIRS: ${MAGMA_INCLUDE_DIRS}\n"
    "  MAGMA_LIBRARY_DIRS: ${MAGMA_LIBRARY_DIRS}\n"
    "  MAGMA_LIBRARIES: ${MAGMA_LIBRARIES}"
  )
  include_directories(${MAGMA_INCLUDE_DIRS})
else()
  message(STATUS
    "MAGMA not found:\n"
  )
endif()

if (MAGMA_FOUND)
  if(BML_CUSOLVER)
    message(STATUS "Will use cuSOLVER")
    add_definitions(-DBML_USE_CUSOLVER)
  endif()
endif()

set(BML_INTERNAL_BLAS FALSE
  CACHE BOOL "Use internal BLAS functions")

set(BML_INTERNAL_GEMM ${BML_INTERNAL_BLAS}
  CACHE BOOL "Use internal version of ?gemm BLAS functions")
if(BML_INTERNAL_GEMM)
  message(STATUS "Will use internal gemm")
  add_definitions(-DBML_INTERNAL_GEMM)
endif()

set(BLAS_VENDOR ""
  CACHE STRING "If set, the preferred BLAS/LAPACK vendor. Possible choices: {Intel,MKL,ACML,GNU,IBM}")

if(BLAS_VENDOR STREQUAL "Intel" OR BLAS_VENDOR STREQUAL "MKL")
  set(MKLROOT $ENV{MKLROOT})
  message(STATUS "Attempting to use Intel's BLAS/LAPACK (MKL)")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Look for MKL in ${MKLROOT}")
    if(BML_OPENMP)
      set(BLAS_LIBRARIES "-L${MKLROOT}/lib/intel64"
        -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
        -lpthread -lm)
    else()
      set(BLAS_LIBRARIES "-L${MKLROOT}/lib/intel64"
        -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm)
    endif()
  else()
    if(BML_OPENMP)
      set(BLAS_LIBRARIES "-L${MKLROOT}/lib/intel64"
        -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp
        -lpthread -lm -ldl)
    else()
      set(BLAS_LIBRARIES -Wl,--no-as-needed "-L${MKLROOT}/lib/intel64"
        -lmkl_lapack95_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
        -lmkl_sequential -lmkl_core -lpthread -lm)
    endif()
  endif()
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
elseif(BLAS_VENDOR STREQUAL "ACML")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if(BML_OPENMP)
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/ifort64_mp/lib" -lacml_mp)
    else()
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/ifort64/lib" -lacml)
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
      OR CMAKE_Fortran_COMPILER_ID STREQUAL "PathScale")
    if(BML_OPENMP)
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/gfortran64_mp/lib" -lacml_mp)
    else()
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/gfortran64/lib" -lacml)
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    if(OPENMP)
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/pgi64_mp/lib" -lacml_mp)
    else()
      set(BLAS_LIBRARIES "-L$ENV{ACML_DIR}/pgi64/lib" -lacml)
    endif()
  endif()
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
elseif(BLAS_VENDOR STREQUAL "GNU")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(BLAS_LIBRARIES  -llapack -lblas)
  endif()
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
elseif(BLAS_VENDOR STREQUAL "IBM")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
    set(BLAS_LIBRARIES "-L$ENV{ESSL_DIR}/lib64" -lesslsmp)
  endif()
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
elseif(BLAS_VENDOR STREQUAL "None")
  set(NOBLAS TRUE)
  add_definitions(-DNOBLAS)
elseif(BLAS_VENDOR STREQUAL "Auto")
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
  message(STATUS "I will skip BLAS/LAPACK detection")
elseif(BLAS_VENDOR STREQUAL "")
  message(STATUS "I will rely on the default paths or environment variables to find BLAS/LAPACK")
else()
  message(FATAL_ERROR "Unknown BLAS_VENDOR ${BLAS_VENDOR}")
endif()
message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
if(BLAS_LIBRARIES)
  set(BLAS_FOUND TRUE)
endif()
if(LAPACK_LIBRARIES)
  set(LAPACK_FOUND TRUE)
endif()

if(NOT (BLAS_FOUND OR NOBLAS))
  include(FindBLAS)
  message(STATUS "BLA_VENDOR=${BLA_VENDOR}")
  message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
  if(NOT BLAS_FOUND)
    message(FATAL_ERROR "Can not find suitable BLAS library")
  endif()
endif()

if(NOT (LAPACK_FOUND OR NOBLAS OR MAGMA_FOUND))
  include(FindLAPACK)
  if(NOT LAPACK_FOUND)
    message(FATAL_ERROR "Can not find suitable LAPACK library")
  endif()
endif()

include(CheckFunctionExists)
include(CheckFortranFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES ${LINK_LIBRARIES} -lm)
check_function_exists(fabs HAVE_FABS)
if(NOT HAVE_FABS)
  message(FATAL_ERROR "Could not find the fabs() function")
endif()
list(APPEND LINK_LIBRARIES -lm)

if(NOT (BLAS_FOUND OR NOBLAS))
  message(FATAL_ERROR "Could not find BLAS library.")
endif()

if(NOT (LAPACK_FOUND OR NOBLAS))
  if(NOT (MAGMA_FOUND) )
    message(FATAL_ERROR "Could not find LAPACK library.")
  endif()
endif()

include(${PROJECT_SOURCE_DIR}/cmake/bmlCheckCFortranFunctionExists.cmake)
if(MAGMA_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES
    ${MAGMA_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${MAGMA_LIBRARIES})
  if(BML_CUSOLVER)
    message(STATUS "cuSOLVER: ${CUDA_cusolver_LIBRARY}")
    list(APPEND LINK_LIBRARIES ${CUDA_cusolver_LIBRARY})
  endif()
endif()
if(BLAS_FOUND)
  add_definitions(-DHAVE_BLAS)

  set(CMAKE_REQUIRED_LIBRARIES
    ${LINK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenMP_Fortran_FLAGS})
  check_fortran_function_exists(sgemm HAVE_SGEMM)
  check_fortran_function_exists(dgemm HAVE_DGEMM)
  check_fortran_function_exists(cgemm HAVE_CGEMM)
  check_fortran_function_exists(zgemm HAVE_ZGEMM)
  if(NOT (HAVE_SGEMM AND HAVE_DGEMM AND HAVE_CGEMM AND HAVE_ZGEMM))
    message(FATAL_ERROR "Could not find {s,d,c,z}gemm")
  endif()
  add_definitions(-DHAVE_SGEMM -DHAVE_DGEMM)

  set(CMAKE_REQUIRED_LIBRARIES
    ${LINK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenMP_C_FLAGS})
  bml_check_C_Fortran_function_exists(dgemm C_DGEMM REQUIRED)
  bml_check_C_Fortran_function_exists(sgemm C_SGEMM REQUIRED)
  bml_check_C_Fortran_function_exists(cgemm C_CGEMM REQUIRED)
  bml_check_C_Fortran_function_exists(zgemm C_ZGEMM REQUIRED)
  bml_check_C_Fortran_function_exists(sscal C_SSCAL REQUIRED)
  bml_check_C_Fortran_function_exists(dscal C_DSCAL REQUIRED)
  bml_check_C_Fortran_function_exists(cscal C_CSCAL REQUIRED)
  bml_check_C_Fortran_function_exists(zscal C_ZSCAL REQUIRED)
  bml_check_C_Fortran_function_exists(saxpy C_SAXPY REQUIRED)
  bml_check_C_Fortran_function_exists(daxpy C_DAXPY REQUIRED)
  bml_check_C_Fortran_function_exists(caxpy C_CAXPY REQUIRED)
  bml_check_C_Fortran_function_exists(zaxpy C_ZAXPY REQUIRED)
  add_definitions(
    -DC_SGEMM=${C_SGEMM}
    -DC_DGEMM=${C_DGEMM}
    -DC_CGEMM=${C_CGEMM}
    -DC_ZGEMM=${C_ZGEMM}
    -DC_SSCAL=${C_SSCAL}
    -DC_DSCAL=${C_DSCAL}
    -DC_CSCAL=${C_CSCAL}
    -DC_ZSCAL=${C_ZSCAL}
    -DC_SAXPY=${C_SAXPY}
    -DC_DAXPY=${C_DAXPY}
    -DC_CAXPY=${C_CAXPY}
    -DC_ZAXPY=${C_ZAXPY})

  list(APPEND LINK_LIBRARIES ${BLAS_LIBRARIES})
endif()

if(LAPACK_FOUND)
  add_definitions(-DHAVE_LAPACK)

  set(CMAKE_REQUIRED_LIBRARIES
    ${LINK_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_Fortran_FLAGS})
  check_fortran_function_exists(ssyev HAVE_SSYEV)
  check_fortran_function_exists(dsyev HAVE_DSYEV)
  check_fortran_function_exists(cheevr HAVE_CHEEVR)
  check_fortran_function_exists(zheevr HAVE_ZHEEVR)
  check_fortran_function_exists(sgetrf HAVE_SGETRF)
  check_fortran_function_exists(dgetrf HAVE_DGETRF)
  check_fortran_function_exists(cgetri HAVE_CGETRI)
  check_fortran_function_exists(zgetri HAVE_ZGETRI)
  if(NOT (HAVE_SSYEV AND HAVE_DSYEV AND HAVE_CHEEVR AND HAVE_ZHEEVR))
    message(FATAL_ERROR "Could not find {s,d}syev or {c,z}heevr")
  endif()
  add_definitions(-DHAVE_SSYEV -DHAVE_DSYEV -DHAVE_CSYEVR -DHAVE_ZSYEVR)

  set(CMAKE_REQUIRED_LIBRARIES
    ${LINK_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_C_FLAGS})
  bml_check_C_Fortran_function_exists(ssyev C_SSYEV REQUIRED)
  bml_check_C_Fortran_function_exists(dsyev C_DSYEV REQUIRED)
  bml_check_C_Fortran_function_exists(cheevr C_CHEEVR REQUIRED)
  bml_check_C_Fortran_function_exists(zheevr C_ZHEEVR REQUIRED)
  bml_check_C_Fortran_function_exists(sgetrf C_SGETRF REQUIRED)
  bml_check_C_Fortran_function_exists(dgetrf C_DGETRF REQUIRED)
  bml_check_C_Fortran_function_exists(cgetrf C_CGETRF REQUIRED)
  bml_check_C_Fortran_function_exists(zgetrf C_ZGETRF REQUIRED)
  bml_check_C_Fortran_function_exists(sgetri C_SGETRI REQUIRED)
  bml_check_C_Fortran_function_exists(dgetri C_DGETRI REQUIRED)
  bml_check_C_Fortran_function_exists(cgetri C_CGETRI REQUIRED)
  bml_check_C_Fortran_function_exists(zgetri C_ZGETRI REQUIRED)
  add_definitions(
    -DC_SSYEV=${C_SSYEV}
    -DC_DSYEV=${C_DSYEV}
    -DC_CHEEVR=${C_CHEEVR}
    -DC_ZHEEVR=${C_ZHEEVR}
    -DC_SGETRF=${C_SGETRF}
    -DC_DGETRF=${C_DGETRF}
    -DC_CGETRF=${C_CGETRF}
    -DC_ZGETRF=${C_ZGETRF}
    -DC_SGETRI=${C_SGETRI}
    -DC_DGETRI=${C_DGETRI}
    -DC_CGETRI=${C_CGETRI}
    -DC_ZGETRI=${C_ZGETRI})

  list(APPEND LINK_LIBRARIES ${LAPACK_LIBRARIES})
endif()

# Check whether the compiler supports complex types. The Nvidia
# compiler is known not to support them for example.

set(BML_COMPLEX TRUE
  CACHE BOOL "Whether to build complex types")

if(BML_COMPLEX)
  try_compile(SUPPORTS_COMPLEX ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/cmake/test_complex.c
    OUTPUT_VARIABLE COMPLEX_OUTPUT
    )
  message(STATUS "Compile output: ${COMPLEX_OUTPUT}")
  if(SUPPORTS_COMPLEX)
    message(STATUS "The compiler supports complex types")
    add_definitions(-DBML_COMPLEX)
  else()
    message(STATUS "The compiler does not support complex types. Will skip them.")
    set(BML_COMPLEX FALSE)
  endif()
else()
  message(STATUS "Complex types disabled by user request")
endif()

set(BML_XSMM FALSE CACHE BOOL "Whether to use XSMM library")
if(NOT (NOBLAS))
  if(BML_XSMM)
    find_package(XSMM)
  endif()
  if(XSMM_FOUND)
    message(STATUS "Use XSMM")
    add_definitions(-DBML_USE_XSMM)
    include_directories(${XSMM_INCLUDE_DIRS})
    list(APPEND LINK_LIBRARIES ${XSMM_LIBRARY_DIRS}/libxsmm.so
                               ${BLAS_LIBRARIES})
  endif()
endif()

if(BML_ELLBLOCK_MEMPOOL)
  message(STATUS "Use memory pool for elllock")
  add_definitions(-DBML_ELLBLOCK_USE_MEMPOOL)
endif()

add_definitions(-D_POSIX_C_SOURCE=200112L)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if(HAVE_POSIX_MEMALIGN)
  add_definitions(-DHAVE_POSIX_MEMALIGN)
endif()

set(INTEL_OPT FALSE
  CACHE BOOL "Whether to use Intel compiler specific optimizations")

if(INTEL_OPT)
  if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    add_definitions(-DINTEL_OPT)
  else()
    message(FATAL_ERROR "Compiler ID must be Intel when INTEL_OPT=yes")
  endif()
endif()

find_program(ADDR2LINE addr2line)
if(ADDR2LINE)
  message(STATUS "Found addr2line, backtrace will resolve line numbers")
  add_definitions(-DHAVE_ADDR2LINE -DADDR2LINE="${ADDR2LINE}")
else()
  message(STATUS "Did not find addr2line")
endif()

add_definitions(-DPROJECT_VERSION="${PROJECT_VERSION}")

set(DIST_FILES CMakeLists.txt)

set(BML_OLD FALSE
  CACHE BOOL "Whether to build the old version of the library")

if(BML_OLD)
  message(STATUS "Building old version of library (unmaintained)")
  add_subdirectory(src-old)
else()
  add_subdirectory(src)
endif()

include(FindDoxygen)

if(DOXYGEN_FOUND)
  if(BML_OLD)
    set(DOXYGEN_INPUT
      ${PROJECT_SOURCE_DIR}/src-old
      ${CMAKE_BINARY_DIR}/src-old)
  else()
    set(DOXYGEN_INPUT
      ${PROJECT_SOURCE_DIR}/src/C-interface
      ${PROJECT_SOURCE_DIR}/src/Fortran-interface)
  endif()
  string(REPLACE ";" " " DOXYGEN_INPUT "${DOXYGEN_INPUT}")
  configure_file(documentation/Doxyfile.in Doxyfile)
  add_custom_target(docs
    COMMAND rm -rvf ${PROJECT_SOURCE_DIR}/docs/API/*
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
else()
  add_custom_target(docs)
endif()

set(BML_TESTING FALSE
  CACHE BOOL "Whether to build the test suite.")
set(BML_VALGRIND FALSE
  CACHE BOOL "Whether to test for memory leaks with valgrind")
set(VALGRIND_COMMON_ARGS --fullpath-after=src/ --fullpath-after=tests/
  --error-exitcode=1 --leak-check=full --show-leak-kinds=all
  --read-var-info=yes --track-origins=yes)

if(BML_TESTING)
  message(STATUS "Setting up test suite")
  find_program(VALGRIND valgrind
    DOC "valgrind - a suite of tools for debugging and profiling programs")
  if(VALGRIND AND BML_VALGRIND)
    message(STATUS "Will test for memory leaks")
  else()
    message(STATUS "Will _not_ test for memory leaks")
  endif()
  enable_testing()
  if(BML_OLD)
    add_subdirectory(tests-old)
  else()
    add_subdirectory(tests)
  endif()
endif()

find_program(TAR tar)
find_program(GIT git)
if(TAR AND GIT)
  execute_process(COMMAND ${GIT} status
    RESULT_VARIABLE IS_GIT_REPOSITORY
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_QUIET
    ERROR_QUIET)
  if(NOT IS_GIT_REPOSITORY EQUAL 0)
    message(WARNING "Not in git repository, disabling dist target")
  else()
    message(STATUS "Creating dist target")
    add_custom_target(dist
      COMMAND ${GIT} archive
      --format=tar.gz
      --prefix=bml-${PROJECT_VERSION}/
      --output=${CMAKE_BINARY_DIR}/bml-${PROJECT_VERSION}.tar.gz
      master
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endif()
else()
  add_custom_target(dist)
endif()

include(FindPythonInterp)
if(NOT PYTHONINTERP_FOUND)
  message(FATAL "Could not find Python interpreter")
endif()
configure_file(${PROJECT_SOURCE_DIR}/scripts/convert-template.in
  ${CMAKE_BINARY_DIR}/convert-template)

message(STATUS "When building executables with the bml library, please use the following link flags")

if(OPENMP_FOUND)
  message(STATUS "Additional Fortran compiler link flags: '${OpenMP_Fortran_FLAGS}'")
  message(STATUS "Additional C compiler link flags: '${OpenMP_C_FLAGS}'")
endif()

message(STATUS "LINK_LIBRARIES=${LINK_LIBRARIES}")
string(REPLACE ";" " " LINKLINE "${LINK_LIBRARIES}")
string(STRIP "${LINKLINE}" LINKLINE)
message(STATUS "Link flags: '${LINKLINE}'")

configure_file(src/bml.pc.in bml.pc)
install(FILES ${CMAKE_BINARY_DIR}/bml.pc
  DESTINATION ${CMAKE_INSTALL_PKG_CONFIG_DIR})
