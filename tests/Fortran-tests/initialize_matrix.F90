program test

  use bml

  implicit none

  integer, parameter :: N = 12, M = 8

  type(bml_matrix_t) :: A

  real, allocatable :: A_dense_single(:, :)
  double precision, allocatable :: A_dense_double(:, :)

  call random_matrix(BML_MATRIX_DENSE, N, A, BML_PRECISION_DOUBLE, M)
  call convert_to_dense(A, A_dense_double)

  call random_matrix(BML_MATRIX_DENSE, N, A, BML_PRECISION_SINGLE, M)
  call convert_to_dense(A, A_dense_single)

  call deallocate_matrix(A)

end program test
