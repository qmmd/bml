program bml_test

  use bml_types_m

  integer :: N = 11
  integer :: M = -1
  character(1000) :: test = ""
  character(1000) :: matrix_type = BML_MATRIX_DENSE
  character(1000) :: matrix_precision = BML_ELEMENT_REAL

  integer :: n_args

  ! Arguments are interpreted as
  !
  ! testname
  ! matrix_type
  ! precision

  n_args = command_argument_count()
  if(n_args /= 3) then
    print *, "incorrect number of command line arguments"
    error stop
  end if

  call get_command_argument(1, test)
  call get_command_argument(2, matrix_type)
  call get_command_argument(3, matrix_precision)

  select case(test)
  case default
    print *, "unknown test name "//trim(test)
    error stop
  end select

end program bml_test
