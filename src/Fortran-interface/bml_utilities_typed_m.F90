!> Utility matrix functions.
module bml_utilities_MATRIX_TYPE_m

  use bml_c_interface_m
  use bml_types_m
  use bml_interface_m

  implicit none
  private

  !> Print a dense matrix.
  interface bml_print_matrix
     module procedure bml_print_dense_matrix_MATRIX_TYPE
  end interface bml_print_matrix

  public :: bml_print_matrix

contains

  !> Print a dense matrix.
  !!
  !! \param tag A string to print before the matrix.
  !! \param a The matrix.
  !! \param i_l The lower row bound.
  !! \param i_u The upper row bound.
  !! \param j_l The lower column bound.
  !! \param j_u The upper column bound.
  subroutine bml_print_dense_matrix_MATRIX_TYPE(tag, a, i_l, i_u, j_l, j_u)

    character(len=*), intent(in) :: tag
    REAL_TYPE, target, intent(in) :: a(:, :)
    integer(C_INT), intent(in) :: i_l
    integer(C_INT), intent(in) :: i_u
    integer(C_INT), intent(in) :: j_l
    integer(C_INT), intent(in) :: j_u

    write(*, "(A)") tag
    associate(a_ptr => a(lbound(a, 1), lbound(a, 2)))
      ! Print bounds are inclusive here, i.e. [i_l, i_u], but are
      ! exclusive in the upper bound in the C code.
      call bml_print_dense_matrix_C(size(a, 1, kind=C_INT), &
           & get_element_id(REAL_NAME, REAL_KIND), &
           & BML_DENSE_COLUMN_MAJOR, &
           & c_loc(a_ptr), &
           & i_l-lbound(a, 1, kind=C_INT), i_u-lbound(a, 1, kind=C_INT)+1, &
           & j_l-lbound(a, 2, kind=C_INT), j_u-lbound(a, 2, kind=C_INT)+1)
    end associate

  end subroutine bml_print_dense_matrix_MATRIX_TYPE

end module bml_utilities_MATRIX_TYPE_m
