
!> General test driver.
!! \brief This program will execute a bml fortran test.
!! \usage testf -n <testName> -t <matrixFormat> -p <precision>
!! testName: Name of the test to be executed.
!! matrixFormat: bml format of the matrix (dense, ellpack, etc.)
!! precision: The element kind and precision (single_real, double_real, etc.)
program testf

  use bml
  use add_matrix
  use adjungate_triangle
  use allocate_matrix
  use convert_matrix
  use copy_matrix
  use diagonalize_matrix
  use get_bandwidth
  use get_diagonal
  use get_element
  use get_row
  use inverse_matrix
  use io_matrix
  use multiply_matrix
  use normalize_matrix
  use scale_matrix
  use set_diagonal
  use set_element
  use set_row
  use threshold_matrix
  use trace_matrix
  use transpose_matrix

  use prec

  implicit none
  character(20) :: args(6)
  integer, parameter :: N = 7, M = 7
  integer :: i, narg
  logical :: missingarg = .false.
  logical :: test_result = .false.
  character(20) :: test_name, matrix_type, element_type

  call get_arguments(test_name, matrix_type, element_type)

  select case(test_name)
  case("add")
    test_result = test_add_matrix(matrix_type, element_type, N, M)
  case("adjungate")
    test_result = test_adjungate_triangle(matrix_type, element_type, N, M)
  case("allocate")
    test_result = test_allocate_matrix(matrix_type, element_type, N, M)
  case("convert")
    test_result = test_convert_matrix(matrix_type, element_type, N, M)
  case("copy")
    test_result = test_copy_matrix(matrix_type, element_type, N, M)
  case("diagonalize")
    test_result = test_diagonalize_matrix(matrix_type, element_type, N, M)
  case("get_bandwidth")
    test_result = test_get_bandwidth(matrix_type, element_type, N, M)
  case("get_diagonal")
    test_result = test_get_diagonal(matrix_type, element_type, N, M)
  case("get_element")
    test_result = test_get_element(matrix_type, element_type, N, M)
  case("get_row")
    test_result = test_get_row(matrix_type, element_type, N, M)
  case("inverse")
    test_result = test_inverse_matrix(matrix_type, element_type, N, M)
  case("io")
    test_result = test_io_matrix(matrix_type, element_type, N, M)
  case("multiply")
    test_result = test_multiply_matrix(matrix_type, element_type, N, M)
  case("normalize")
    test_result = test_normalize_matrix(matrix_type, element_type, N, M)
  case("scale")
    test_result = test_scale_matrix(matrix_type, element_type, N, M)
  case("set_diagonal")
    test_result = test_set_diagonal(matrix_type, element_type, N, M)
  case("set_element")
    test_result = test_set_element(matrix_type, element_type, N, M)
  case("set_row")
    test_result = test_set_row(matrix_type, element_type, N, M)
  case("threshold")
    test_result = test_threshold_matrix(matrix_type, element_type, N, M)
  case("trace")
    test_result = test_trace_matrix(matrix_type, element_type, N, M)
  case("transpose")
    test_result = test_transpose_matrix(matrix_type, element_type, N, M)
  case default
    call error_usage("Invalid test name "//test_name)
  end select

  if(test_result)then
    write(*,*)"Test passed"
  else
    write(*,*)"Test failed"
    error stop
  endif

end program testf

!> Error messages for this program.
!! \param message Message to be passed.
subroutine error_usage(message)
  character(*) :: message

  write(*,*)""
  write(*,*)"ERROR:"
  write(*,*) message
  write(*,*)""
  write(*,*)"USAGE:"
  write(*,*)"testf -n <testName> -t <matrixFormat> -p <precision>"
  write(*,*)""
  error stop

end subroutine error_usage

!> Get the arguments from command line.
!! \brief Gets the arguments that are passed from command line.
!! test_name: Name of the test to be executed.
!! matrix_type: bml format of the matrix (dense, ellpack, etc.)
!! element_type: The element kind and precision (single_real, double_real, etc.)
subroutine get_arguments(test_name, matrix_type, element_type)
  character(20) :: args(6)
  integer :: i, narg
  logical :: missingarg = .false.
  character(20), intent(out) :: test_name, matrix_type, element_type

  do i = 1,6
    call getarg(i,args(i))
  end do

  write(*,*)args
  narg = 0
  do i = 1,6
    if (trim(adjustl(args(i))) == "") then
      call error_usage("No arguments passed")
    end if
    if (trim(adjustl(args(i))) == "-n") then
      if (trim(adjustl(args(i+1))) .NE. "") then
        test_name = args(i+1)
        narg = narg + 1
      end if
    end if

    if (trim(adjustl(args(i))) == "-t") then
      if (trim(adjustl(args(i+1))) .NE. "") then
        matrix_type = args(i+1)
        narg = narg + 1
      end if
    end if

    if (trim(adjustl(args(i))) == "-p") then
      if (trim(adjustl(args(i+1))) .NE. "") then
        element_type = args(i+1)
        narg = narg + 1
      end if
    end if
  end do

  if  (narg < 3) then
    write(*,*)"Number of arguments",narg
    call error_usage("Missing an argument")
  endif
  write(*,*)"End parsing"

end subroutine get_arguments
