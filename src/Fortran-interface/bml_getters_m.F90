module bml_getters_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  interface bml_get_row
    module procedure bml_get_row_single_real
    module procedure bml_get_row_double_real
    module procedure bml_get_row_single_complex
    module procedure bml_get_row_double_complex
  end interface bml_get_row

  interface bml_get_diagonal
    module procedure bml_get_diagonal_single_real
    module procedure bml_get_diagonal_double_real
    module procedure bml_get_diagonal_single_complex
    module procedure bml_get_diagonal_double_complex
  end interface bml_get_diagonal

  public :: bml_get_row, bml_get_diagonal

contains

  !Getters for diagonal

  !> Get the diagonal i of matrix a
  !! \param a The matrix
  !! \param diagonal The diagonal that is extracted
  subroutine bml_get_diagonal_single_real(a, diagonal)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    real(C_FLOAT), allocatable, intent(inout) :: diagonal(:)
    real(C_FLOAT), pointer :: diagonal_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_diagonal_C(a%ptr)
    call c_f_pointer(ptr, diagonal_ptr, [bml_get_N(a)])
    diagonal = diagonal_ptr

  end subroutine bml_get_diagonal_single_real

  !> Get the diagonal i of matrix a
  !! \param a The matrix
  !! \param diagonal The diagonal that is extracted
  subroutine bml_get_diagonal_double_real(a, diagonal)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    real(C_DOUBLE), allocatable, intent(inout) :: diagonal(:)
    real(C_DOUBLE), pointer :: diagonal_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_diagonal_C(a%ptr)
    call c_f_pointer(ptr, diagonal_ptr, [bml_get_N(a)])
    diagonal = diagonal_ptr

  end subroutine bml_get_diagonal_double_real

  !> Get the diagonal i of matrix a
  !! \param a The matrix
  !! \param diagonal The diagonal that is extracted
  subroutine bml_get_diagonal_single_complex(a, diagonal)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    complex(C_FLOAT_COMPLEX), allocatable, intent(inout) :: diagonal(:)
    complex(C_FLOAT_COMPLEX), pointer :: diagonal_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_diagonal_C(a%ptr)
    call c_f_pointer(ptr, diagonal_ptr, [bml_get_N(a)])
    diagonal = diagonal_ptr

  end subroutine bml_get_diagonal_single_complex

  !> Get the diagonal i of matrix a
  !! \param a The matrix
  !! \param diagonal The diagonal that is extracted
  subroutine bml_get_diagonal_double_complex(a, diagonal)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    complex(C_DOUBLE_COMPLEX), allocatable, intent(inout) :: diagonal(:)
    complex(C_DOUBLE_COMPLEX), pointer :: diagonal_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_diagonal_C(a%ptr)
    call c_f_pointer(ptr, diagonal_ptr, [bml_get_N(a)])
    diagonal = diagonal_ptr

  end subroutine bml_get_diagonal_double_complex

  !Getter for row

  !> Get the row i of matrix a
  !! \param a The matrix
  !! \param i The row number
  !! \param row The row that is extracted
  subroutine bml_get_row_single_real(a, i, row)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), intent(in) :: i
    real(C_FLOAT), allocatable, intent(inout) :: row(:)
    real(C_FLOAT), pointer :: row_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_row_C(a%ptr, i - 1)
    call c_f_pointer(ptr, row_ptr, [bml_get_N(a)])
    row = row_ptr

  end subroutine bml_get_row_single_real

  !> Get the row i of matrix a
  !! \param a The matrix
  !! \param i The row number
  !! \param row The row that is extracted
  subroutine bml_get_row_double_real(a, i, row)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), intent(in) :: i
    real(C_DOUBLE), allocatable, intent(inout) :: row(:)
    real(C_DOUBLE), pointer :: row_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_row_C(a%ptr, i - 1)
    call c_f_pointer(ptr, row_ptr, [bml_get_N(a)])
    row = row_ptr

  end subroutine bml_get_row_double_real

  !> Get the row i of matrix a
  !! \param a The matrix
  !! \param i The row number
  !! \param row The row that is extracted
  subroutine bml_get_row_single_complex(a, i, row)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), intent(in) :: i
    complex(C_FLOAT_COMPLEX), allocatable, intent(inout) :: row(:)
    complex(C_FLOAT_COMPLEX), pointer :: row_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_row_C(a%ptr, i - 1)
    call c_f_pointer(ptr, row_ptr, [bml_get_N(a)])
    row = row_ptr

  end subroutine bml_get_row_single_complex

  !> Get the row i of matrix a
  !! \param a The matrix
  !! \param i The row number
  !! \param row The row that is extracted
  subroutine bml_get_row_double_complex(a, i, row)

    use bml_introspection_m

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), intent(in) :: i
    complex(C_DOUBLE_COMPLEX), allocatable, intent(inout) :: row(:)
    complex(C_DOUBLE_COMPLEX), pointer :: row_ptr(:)
    type(C_PTR) :: ptr

    ptr = bml_get_row_C(a%ptr, i - 1)
    call c_f_pointer(ptr, row_ptr, [bml_get_N(a)])
    row = row_ptr

  end subroutine bml_get_row_double_complex

end module bml_getters_m
