module io_matrix_typed

#ifdef INTEL_SDK
  use ifport
#endif

  use bml
  use prec
  use bml_utilities_m

  implicit none

  public :: test_io_matrix_typed

#ifdef CRAY_SDK
  interface
    integer(c_int) function getpid() bind(c,name="getpid")
      use iso_c_binding
    end function getpid
  end interface
#endif

contains

  function test_io_matrix_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    real(DUMMY_PREC) :: tol
    integer :: pid
    character(20) :: pid_char
    character(100) :: fname


    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    DUMMY_KIND(DUMMY_PREC), allocatable :: b_dense(:, :)

    pid = getpid()
    write(pid_char,*)pid
    write(fname,*)"ftest_matrix",trim(adjustl(pid_char)),".mtx"

    allocate(a_dense(n,n))
    allocate(b_dense(n,n))

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    call bml_write_matrix(a, fname)
    call bml_zero_matrix(matrix_type, element_kind, element_precision, n, m, b)
    call bml_read_matrix(b, fname)

    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(b, b_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_print_matrix("B", b_dense, 1, n, 1, n)
    if(maxval(abs(a_dense - b_dense)) > tol) then
      test_result = .false.
      print *, "matrices are not identical"
    else
      test_result = .true.
    end if

  end function test_io_matrix_typed

end module io_matrix_typed
