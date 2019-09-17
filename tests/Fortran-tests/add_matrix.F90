module add_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: add_matrix_t
  contains
    procedure, nopass :: test_function
  end type add_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    double precision, parameter :: alpha = 1.2
    double precision, parameter :: beta = 0.8
    double precision, parameter :: threshold = 0.0

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
    double precision :: rel_tol = 1e-6
#else
    double precision :: rel_tol = 1d-12
#endif

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b
    type(bml_matrix_t) :: c
    type(bml_matrix_t) :: d

    REAL_TYPE, allocatable :: a_dense(:, :)
    REAL_TYPE, allocatable :: b_dense(:, :)
    REAL_TYPE, allocatable :: c_dense(:, :)
    REAL_TYPE, allocatable :: d_dense(:, :)

    double precision :: expected, rel_diff
    integer :: i, j

    test_result = .true.

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_copy(a, b)
    call bml_copy(a, d)
    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & c)

    call bml_add(alpha, b, beta, c, threshold)
    call bml_add_identity(d, alpha, threshold)

    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(b, b_dense)
    call bml_export_to_dense(c, c_dense)
    call bml_export_to_dense(d, d_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_print_matrix("B", b_dense, 1, n, 1, n)
    call bml_print_matrix("C", c_dense, 1, n, 1, n)
    call bml_print_matrix("D", d_dense, 1, n, 1, n)
    do i = 1, n
      do j = 1, n
        expected = alpha * a_dense(i, j) + beta * c_dense(i, j)
        rel_diff = abs((expected - b_dense(i, j)) / expected)
        if(rel_diff > rel_tol) then
          print *, "rel. diff = ", rel_diff
          call bml_error(__FILE__, __LINE__, "add() matrices are not identical")
        end if

        if(i == j) then
          expected = a_dense(i, j) + alpha
        else
          expected = a_dense(i, j)
        end if
        rel_diff = abs((expected - d_dense(i, j)) / expected)
        if(rel_diff > rel_tol) then
          print *, "rel. diff = ", rel_diff
          call bml_error(__FILE__, __LINE__, "add_identity() matrices are not identical")
        end if
      end do
    end do

    call bml_deallocate(a)
    call bml_deallocate(b)
    call bml_deallocate(c)
    call bml_deallocate(d)

    deallocate(a_dense)
    deallocate(b_dense)
    deallocate(c_dense)
    deallocate(d_dense)

  end function test_function

end module add_matrix_m
