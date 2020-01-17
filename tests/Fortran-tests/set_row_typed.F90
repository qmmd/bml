module set_row_typed

  use bml
  use prec
  use bml_setters_m

  implicit none

  public :: test_set_row_typed

contains

  function test_set_row_typed(matrix_type, element_kind, element_precision,&
       & n, m) result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    integer :: i, j
    logical :: test_result
    type(bml_matrix_t) :: a
    DUMMY_KIND(DUMMY_PREC), allocatable :: row(:), a_dense(:,:)
    real(dp), parameter :: threshold = 0.0_dp
    real(DUMMY_PREC) :: diff, tol

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    test_result = .true.

    allocate(row(n))
    allocate(a_dense(n,n))

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)

    do i=1,n
      do j=1,n
        row(j) = real(i,DUMMY_PREC)*real(j,DUMMY_PREC)
      end do
      call bml_set_row(a, i, row, threshold)
    end do

    call bml_export_to_dense(a,a_dense)

    do i=1,n
      do j=1,n
        diff = abs(a_dense(i,j) - real(i,DUMMY_PREC)*real(j,DUMMY_PREC))
        if(diff > tol) then
          print *, "diff = ", diff
          call bml_error(__FILE__, __LINE__, "set_row failing")
          test_result = .false.
        end if
      end do
    end do

    deallocate(a_dense)
    deallocate(row)
    call bml_deallocate(a)

  end function test_set_row_typed

end module set_row_typed
