
!> Routine to set elements, rows and diagonals.
!!
module bml_setters_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private


  !> Routine to set elements to a new matrix
  !! This routine is just faster than set_element
  !! We just need to make sure that the matrix is empty and
  !! no elements are replaced.
  !! The dense version remains the same as set_element.
  interface bml_set_element_new
    module procedure bml_set_element_new_single_real
    module procedure bml_set_element_new_double_real
    module procedure bml_set_element_new_single_complex
    module procedure bml_set_element_new_double_complex
  end interface bml_set_element_new

  !> Routine to set elements.
  !!
  interface bml_set_element
    module procedure bml_set_element_single_real
    module procedure bml_set_element_double_real
    module procedure bml_set_element_single_complex
    module procedure bml_set_element_double_complex
  end interface bml_set_element

  !> Routine to set rows.
  !!
  interface bml_set_row
    module procedure bml_set_row_single_real
    module procedure bml_set_row_double_real
    module procedure bml_set_row_single_complex
    module procedure bml_set_row_double_complex
  end interface bml_set_row

  !> Routine to set diagonal.
  !!
  interface bml_set_diagonal
    module procedure bml_set_diagonal_single_real
    module procedure bml_set_diagonal_double_real
    module procedure bml_set_diagonal_single_complex
    module procedure bml_set_diagonal_double_complex
  end interface bml_set_diagonal

  public :: bml_set_row
  public :: bml_set_diagonal
  public :: bml_set_element
  public :: bml_set_element_new

contains

  !Setters for element new
  subroutine bml_set_element_new_single_real(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    real(C_FLOAT), target, intent(in) :: element

    call bml_set_element_new_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_new_single_real

  subroutine bml_set_element_new_double_real(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    real(C_DOUBLE), target, intent(in) :: element

    call bml_set_element_new_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_new_double_real

  subroutine bml_set_element_new_single_complex(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    complex(C_FLOAT_COMPLEX), target, intent(in) :: element

    call bml_set_element_new_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_new_single_complex

  subroutine bml_set_element_new_double_complex(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    complex(C_DOUBLE_COMPLEX), target, intent(in) :: element

    call bml_set_element_new_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_new_double_complex

  !Setters for element
  subroutine bml_set_element_single_real(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    real(C_FLOAT), target, intent(in) :: element

    call bml_set_element_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_single_real

  subroutine bml_set_element_double_real(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    real(C_DOUBLE), target, intent(in) :: element

    call bml_set_element_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_double_real

  subroutine bml_set_element_single_complex(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    complex(C_FLOAT_COMPLEX), target, intent(in) :: element

    call bml_set_element_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_single_complex

  subroutine bml_set_element_double_complex(a,i,j,element)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    integer(C_INT), intent(in) :: j
    complex(C_DOUBLE_COMPLEX), target, intent(in) :: element

    call bml_set_element_C(a%ptr,i-1,j-1,c_loc(element))

  end subroutine bml_set_element_double_complex


  !Setters for row
  subroutine bml_set_row_single_real(a, i, row, threshold)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    real(C_FLOAT), target, intent(in) :: row(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_row_C(a%ptr, i-1, c_loc(row), threshold_)

  end subroutine bml_set_row_single_real

  subroutine bml_set_row_double_real(a, i, row, threshold)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    real(C_DOUBLE), target, intent(in) :: row(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_row_C(a%ptr, i-1, c_loc(row), threshold_)

  end subroutine bml_set_row_double_real

  subroutine bml_set_row_single_complex(a, i, row, threshold)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    complex(C_FLOAT_COMPLEX), target, intent(in) :: row(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_row_C(a%ptr, i-1, c_loc(row), threshold_)

  end subroutine bml_set_row_single_complex

  subroutine bml_set_row_double_complex(a, i, row, threshold)

    type(bml_matrix_t), intent(inout) :: a
    integer(C_INT), intent(in) :: i
    complex(C_DOUBLE_COMPLEX), target, intent(in) :: row(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_row_C(a%ptr, i-1, c_loc(row), threshold_)

  end subroutine bml_set_row_double_complex

  !Setter for diagonal
  subroutine bml_set_diagonal_single_real(a, diagonal, threshold)

    type(bml_matrix_t), intent(inout) :: a
    real(C_FLOAT), target, intent(in) :: diagonal(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_diagonal_C(a%ptr, c_loc(diagonal), threshold_)

  end subroutine bml_set_diagonal_single_real

  subroutine bml_set_diagonal_double_real(a, diagonal, threshold)

    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), target, intent(in) :: diagonal(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_diagonal_C(a%ptr, c_loc(diagonal), threshold_)

  end subroutine bml_set_diagonal_double_real

  subroutine bml_set_diagonal_single_complex(a, diagonal, threshold)

    type(bml_matrix_t), intent(inout) :: a
    complex(C_FLOAT_COMPLEX), target, intent(in) :: diagonal(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_diagonal_C(a%ptr,c_loc(diagonal), threshold_)

  end subroutine bml_set_diagonal_single_complex

  subroutine bml_set_diagonal_double_complex(a, diagonal, threshold)

    type(bml_matrix_t), intent(inout) :: a
    complex(C_DOUBLE_COMPLEX), target, intent(in) :: diagonal(*)
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0d0
    end if

    call bml_set_diagonal_C(a%ptr,c_loc(diagonal), threshold_)

  end subroutine bml_set_diagonal_double_complex

end module bml_setters_m
