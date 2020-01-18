module add_matrix_typed

  use bml
  use prec
  use bml_add_m

  implicit none

  public :: test_add_matrix_typed

contains

  function test_add_matrix_typed(matrix_type, element_kind, element_precision,&
       & n, m) result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    !! \todo Fixme: For some reason the test fails when single prec is set
    real(dp), parameter :: alpha = 1.2_dp
    real(dp), parameter :: beta = 0.8_dp
    real(dp), parameter :: threshold = 0.0_dp
    real(dp) :: trace

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b
    type(bml_matrix_t) :: c
    type(bml_matrix_t) :: d

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    DUMMY_KIND(DUMMY_PREC), allocatable :: b_dense(:, :)
    DUMMY_KIND(DUMMY_PREC), allocatable :: c_dense(:, :)
    DUMMY_KIND(DUMMY_PREC), allocatable :: d_dense(:, :)

    real(DUMMY_PREC) :: expected, rel_diff, tol

    integer :: i, j

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    test_result = .true.

    !c = a + b
    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    call bml_copy_new(a, b)
    call bml_copy_new(a, d)
    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & c)

    call bml_add(b, c, alpha, beta, threshold)

    call bml_add_identity(d, alpha, threshold)


    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(b, b_dense)
    call bml_export_to_dense(c, c_dense)
    call bml_export_to_dense(d, d_dense)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_print_matrix("B", b, 1, n, 1, n)
    call bml_print_matrix("C", c, 1, n, 1, n)
    call bml_print_matrix("D", d, 1, n, 1, n)
    do i = 1, n
      do j = 1, n
        expected = alpha * a_dense(i, j) + beta * c_dense(i, j)
        rel_diff = abs((expected - b_dense(i, j)) / expected)
        if(rel_diff > tol) then
          print *, "rel. diff = ", rel_diff
          call bml_error(__FILE__, __LINE__, "add() matrices are not identical")
        end if

        if(i == j) then
          expected = a_dense(i, j) + alpha
        else
          expected = a_dense(i, j)
        end if
        rel_diff = abs((expected - d_dense(i, j)) / expected)
        if(rel_diff > tol) then
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

  end function test_add_matrix_typed

end module add_matrix_typed
