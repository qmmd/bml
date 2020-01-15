module adjungate_triangle_typed

  use bml
  use prec
  use bml_adjungate_triangle_m

  implicit none

  public :: test_adjungate_triangle_typed

contains

  function test_adjungate_triangle_typed(matrix_type, element_kind, element_precision,&
       & n, m) result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    integer :: i,j
    logical :: test_result
    real(DUMMY_PREC) :: aux, tol
    complex(DUMMY_PREC) :: aji

    type(bml_matrix_t) :: a

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)

    test_result = .true.

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)

    call bml_adjungate_triangle(a, "u")

    call bml_export_to_dense(a, a_dense)

    call bml_deallocate(a)

    do i=1,n
      do j=i+1,n
        aji = a_dense(j,i)
        aux = abs(a_dense(j,i) - conjg(aji))
        if(aux > tol)then
          test_result = .false.
        end if
      end do
    end do

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)

    call bml_adjungate_triangle(a, "l")

    call bml_export_to_dense(a, a_dense)

    call bml_deallocate(a)

    do i=1,n
      do j=i+1,n
        aji = a_dense(j,i)
        aux = abs(a_dense(j,i) - conjg(aji))
        if(aux > tol)then
          test_result = .false.
          write(*,*)"aij is not equal to conjugat(aji)",a_dense(i,j),aji
        end if
      end do
    end do

    deallocate(a_dense)

  end function test_adjungate_triangle_typed

end module adjungate_triangle_typed
