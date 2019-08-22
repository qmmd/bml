module scale_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: scale_matrix_t
  contains
    procedure, nopass :: test_function
  end type scale_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    double precision, parameter :: alpha = 1.2

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
    real :: abs_tol = 1e-6
#else
    double precision :: abs_tol = 1d-12
#endif

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: c

    REAL_TYPE, allocatable :: a_dense(:, :)
    REAL_TYPE, allocatable :: c_dense(:, :)

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, c)
    call bml_scale(alpha, a, c)

    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(c, c_dense)

    if(maxval(abs(alpha * a_dense - c_dense)) > abs_tol) then
      test_result = .false.
      call bml_print_matrix("A", alpha * a_dense, 1, n, 1, n)
      call bml_print_matrix("C", c_dense, 1, n, 1, n)
      print *, "maxval abs difference = ", maxval(abs(alpha * a_dense - c_dense))
      print *, "abs_tol = ", abs_tol
      print *, "matrix element mismatch"
    else
      test_result = .true.
    endif

    call bml_deallocate(a)
    call bml_deallocate(c)

    deallocate(a_dense)
    deallocate(c_dense)

  end function test_function

end module scale_matrix_m
