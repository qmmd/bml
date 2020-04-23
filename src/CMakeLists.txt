add_subdirectory(Fortran-interface)
add_subdirectory(C-interface)
add_subdirectory(internal-blas)

add_library(bml
  $<TARGET_OBJECTS:bml-c>
  $<TARGET_OBJECTS:bml-internal-blas-double_complex>
  $<TARGET_OBJECTS:bml-internal-blas-double_real>
  $<TARGET_OBJECTS:bml-internal-blas-single_complex>
  $<TARGET_OBJECTS:bml-internal-blas-single_real>
  $<TARGET_OBJECTS:bml-dense-double_complex>
  $<TARGET_OBJECTS:bml-dense-double_real>
  $<TARGET_OBJECTS:bml-dense-single_complex>
  $<TARGET_OBJECTS:bml-dense-single_real>
  $<TARGET_OBJECTS:bml-dense>
  $<TARGET_OBJECTS:bml-ellpack-double_complex>
  $<TARGET_OBJECTS:bml-ellpack-double_real>
  $<TARGET_OBJECTS:bml-ellpack-single_complex>
  $<TARGET_OBJECTS:bml-ellpack-single_real>
  $<TARGET_OBJECTS:bml-ellpack>
  $<TARGET_OBJECTS:bml-ellblock-double_complex>
  $<TARGET_OBJECTS:bml-ellblock-double_real>
  $<TARGET_OBJECTS:bml-ellblock-single_complex>
  $<TARGET_OBJECTS:bml-ellblock-single_real>
  $<TARGET_OBJECTS:bml-ellblock>
  $<TARGET_OBJECTS:bml-ellsort-double_complex>
  $<TARGET_OBJECTS:bml-ellsort-double_real>
  $<TARGET_OBJECTS:bml-ellsort-single_complex>
  $<TARGET_OBJECTS:bml-ellsort-single_real>
  $<TARGET_OBJECTS:bml-ellsort>)
set_target_properties(bml
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  POSITION_INDEPENDENT_CODE yes)
target_include_directories(bml PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(OPENMP_FOUND)
  target_link_libraries(bml PUBLIC OpenMP::OpenMP_C)
endif()
target_link_libraries(bml PRIVATE ${LINK_LIBRARIES})

add_library(bml_fortran
  $<TARGET_OBJECTS:bml-fortran>)
set_target_properties(bml_fortran
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  POSITION_INDEPENDENT_CODE yes)
target_include_directories(bml_fortran PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/Fortran-interface>)
target_link_libraries(bml_fortran PUBLIC bml)
if(OPENMP_FOUND)
  target_link_libraries(bml_fortran PUBLIC OpenMP::OpenMP_Fortran)
endif()

install(TARGETS bml bml_fortran EXPORT BML_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT BML_Targets FILE BML_Targets.cmake NAMESPACE BML:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BML)
include(CMakePackageConfigHelpers)
configure_file(BMLConfig.cmakein ${CMAKE_CURRENT_BINARY_DIR}/BMLConfig.cmake @ONLY)
write_basic_package_version_file("BMLConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BMLConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/BMLConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BML)
add_library(BML::bml ALIAS bml)
add_library(BML::bml_fortran ALIAS bml_fortran)
