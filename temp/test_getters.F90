!#define SINGLE_REAL
!#define DOUBLE_REAL
!#define SINGLE_COMPLEX
#define DOUBLE_COMPLEX

#ifdef SINGLE_REAL
#define REAL_T real(kind(0e0))
#define BML_T single_real
#elif defined(DOUBLE_REAL)
#define REAL_T real(kind(0d0))
#define BML_T double_real
#elif defined(SINGLE_COMPLEX)
#define REAL_T complex(kind(0e0))
#define BML_T single_complex
#elif defined(DOUBLE_COMPLEX)
#define REAL_T complex(kind(0d0))
#define BML_T double_complex
#endif

program test_getters

  use bml

  double precision, parameter :: REL_TOL = 1d-12
  integer, parameter :: N = 4

  type(bml_matrix_t) :: A
  real(kind(0d0)) :: A_real(N, N)
  REAL_T :: A_dense(N, N)
  REAL_T :: row(N)

  integer :: i, j

  call random_number(A_real)
  A_dense = A_real
  write(*, "(A)") "A_dense ="
  call print_dense_matrix(A_dense)

  call bml_convert_from_dense(BML_MATRIX_DENSE, A_dense, A)
  call bml_print_matrix("A", A, 1, N, 1, N)

  do i = 1, N
    write (*, *) "Getting row", i
    call bml_get_row(A, i, row)
    call print_dense_vector(row)
    do j = 1, N
      rel_diff = abs((A_dense(i, j) - row(j)) / A_dense(i, j))
      if (rel_diff > REL_TOL) then
        write(*, "(A,I2,',',I2,A)") "matrices are not identical at A(", i, j, ")"
        write(*, "(A,I2)") "getting row ", i
        call print_dense_vector(row)
        error stop
      end if
    end do
  end do
  write(*, *) "test passed"

contains

  subroutine print_dense_vector(x)

    REAL_T, intent(in) :: x(:)

    integer :: i
    character(len=1000) :: format_string

    format_string = "("
    do i = 1, size(x, 1)
#if defined(SINGLE_REAL) || defined(DOUBLE_REAL)
      format_string = trim(format_string)//"f7.3"
#else
      format_string = trim(format_string)//"f7.3,f7.3"
#endif
      if (i < size(x, 1)) then
        format_string = trim(format_string)//","
      end if
    end do
    format_string = trim(format_string)//")"
    write(*, format_string) (x(i), i = 1, size(x, 1))

  end subroutine print_dense_vector

  subroutine print_dense_matrix(A)

    REAL_T, intent(in) :: A(:, :)

    integer :: i
    character(len=1000) :: format_string

    format_string = "("
    do i = 1, size(A, 2)
#if defined(SINGLE_REAL) || defined(DOUBLE_REAL)
      format_string = trim(format_string)//"f7.3"
#else
      format_string = trim(format_string)//"f7.3,f7.3"
#endif
      if (i < size(A, 2)) then
        format_string = trim(format_string)//","
      end if
    end do
    format_string = trim(format_string)//")"
    do i = 1, size(A, 1)
      write(*, format_string) (A(i, j), j = 1, size(A, 2))
    end do

  end subroutine print_dense_matrix

end program test_getters
