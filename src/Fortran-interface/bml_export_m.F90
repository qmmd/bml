module bml_export_m

  use bml_c_interface_m
  use bml_types_m
  use bml_interface_m
  use bml_introspection_m

  implicit none
  private

  interface bml_export_to_dense
    module procedure bml_export_to_dense_single
    module procedure bml_export_to_dense_double
    module procedure bml_export_to_dense_single_complex
    module procedure bml_export_to_dense_double_complex
  end interface bml_export_to_dense

  public :: bml_export_to_dense

contains

  !> Convert a matrix into a dense matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param a The bml matrix
  !! \param a_dense The dense matrix
  subroutine bml_export_to_dense_single(a, a_dense)

    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    real(C_FLOAT), allocatable, intent(inout) :: a_dense(:, :)

    type(C_PTR) :: a_ptr
    real(C_FLOAT), pointer :: a_dense_ptr(:, :)

    a_ptr = bml_export_to_dense_C(a%ptr, BML_DENSE_COLUMN_MAJOR)
    call c_f_pointer(a_ptr, a_dense_ptr, [bml_get_n(a), bml_get_n(a)])
    a_dense = a_dense_ptr
    call bml_free(a_ptr)

  end subroutine bml_export_to_dense_single

  !> Convert a matrix into a dense matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param a The bml matrix
  !! \param a_dense The dense matrix
  subroutine bml_export_to_dense_double(a, a_dense)

    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    real(C_DOUBLE), allocatable, intent(inout) :: a_dense(:, :)

    type(C_PTR) :: a_ptr
    real(C_DOUBLE), pointer :: a_dense_ptr(:, :)

    a_ptr = bml_export_to_dense_C(a%ptr, BML_DENSE_COLUMN_MAJOR)
    call c_f_pointer(a_ptr, a_dense_ptr, [bml_get_n(a), bml_get_n(a)])
    a_dense = a_dense_ptr
    call bml_free(a_ptr)

  end subroutine bml_export_to_dense_double

  !> Convert a matrix into a dense matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param a The bml matrix
  !! \param a_dense The dense matrix
  subroutine bml_export_to_dense_single_complex(a, a_dense)

    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    complex(C_FLOAT_COMPLEX), allocatable, intent(out) :: a_dense(:, :)

    type(C_PTR) :: a_ptr
    complex(C_FLOAT_COMPLEX), pointer :: a_dense_ptr(:, :)

    a_ptr = bml_export_to_dense_C(a%ptr, BML_DENSE_COLUMN_MAJOR)
    call c_f_pointer(a_ptr, a_dense_ptr, [bml_get_n(a), bml_get_n(a)])
    a_dense = a_dense_ptr
    call bml_free(a_ptr)

  end subroutine bml_export_to_dense_single_complex

  !> Convert a matrix into a dense matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param a The bml matrix
  !! \param a_dense The dense matrix
  subroutine bml_export_to_dense_double_complex(a, a_dense)

    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    complex(C_DOUBLE_COMPLEX), allocatable, intent(out) :: a_dense(:, :)

    type(C_PTR) :: a_ptr
    complex(C_DOUBLE_COMPLEX), pointer :: a_dense_ptr(:, :)

    a_ptr = bml_export_to_dense_C(a%ptr, BML_DENSE_COLUMN_MAJOR)
    call c_f_pointer(a_ptr, a_dense_ptr, [bml_get_n(a), bml_get_n(a)])
    a_dense = a_dense_ptr
    call bml_free(a_ptr)

  end subroutine bml_export_to_dense_double_complex

end module bml_export_m
