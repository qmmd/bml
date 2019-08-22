!> Matrix scaling for matrices.
module bml_scale_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  !> Scale a matrix.
  interface bml_scale
    module procedure scale_one_single_real
    module procedure scale_one_double_real
    module procedure scale_one_single_complex
    module procedure scale_one_double_complex
    module procedure scale_two_single_real
    module procedure scale_two_double_real
    module procedure scale_two_single_complex
    module procedure scale_two_double_complex
  end interface bml_scale

  public :: bml_scale

contains

  !> Scale a bml matrix.
  !!
  !! \f$ A \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  subroutine scale_one_single_real(alpha, a)

    real(C_FLOAT), target, intent(in) :: alpha
    type(bml_matrix_t), intent(inout) :: a

    call bml_scale_inplace_C(c_loc(alpha), a%ptr)

  end subroutine scale_one_single_real

  !> Scale a bml matrix.
  !!
  !! \f$ A \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  subroutine scale_one_double_real(alpha, a)

    real(C_DOUBLE), target, intent(in) :: alpha
    type(bml_matrix_t), intent(inout) :: a

    call bml_scale_inplace_C(c_loc(alpha), a%ptr)

  end subroutine scale_one_double_real

  !> Scale a bml matrix.
  !!
  !! \f$ A \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  subroutine scale_one_single_complex(alpha, a)

    complex(C_FLOAT_COMPLEX), target, intent(in) :: alpha
    type(bml_matrix_t), intent(inout) :: a

    call bml_scale_inplace_C(c_loc(alpha), a%ptr)

  end subroutine scale_one_single_complex

  !> Scale a bml matrix.
  !!
  !! \f$ A \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  subroutine scale_one_double_complex(alpha, a)

    complex(C_DOUBLE_COMPLEX), target, intent(in) :: alpha
    type(bml_matrix_t), intent(inout) :: a

    call bml_scale_inplace_C(c_loc(alpha), a%ptr)

  end subroutine scale_one_double_complex

  !> Scale a bml matrix.
  !!
  !! \f$ C \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  !! \param c The matrix
  subroutine scale_two_single_real(alpha, a, c)

    use bml_introspection_m

    real(C_FLOAT), target, intent(in) :: alpha
    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: c

    call bml_scale_C(c_loc(alpha), a%ptr, c%ptr)

  end subroutine scale_two_single_real

  !> Scale a bml matrix.
  !!
  !! \f$ C \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  !! \param c The matrix
  subroutine scale_two_double_real(alpha, a, c)

    use bml_introspection_m

    real(C_DOUBLE), target, intent(in) :: alpha
    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: c

    call bml_scale_C(c_loc(alpha), a%ptr, c%ptr)

  end subroutine scale_two_double_real

  !> Scale a bml matrix.
  !!
  !! \f$ C \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  !! \param c The matrix
  subroutine scale_two_single_complex(alpha, a, c)

    use bml_introspection_m

    complex(C_FLOAT_COMPLEX), target, intent(in) :: alpha
    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: c

    call bml_scale_C(c_loc(alpha), a%ptr, c%ptr)

  end subroutine scale_two_single_complex

  !> Scale a bml matrix.
  !!
  !! \f$ C \leftarrow \alpha A \f$
  !!
  !! \param alpha The factor
  !! \param a The matrix
  !! \param c The matrix
  subroutine scale_two_double_complex(alpha, a, c)

    use bml_introspection_m

    complex(C_DOUBLE_COMPLEX), target, intent(in) :: alpha
    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: c

    call bml_scale_C(c_loc(alpha), a%ptr, c%ptr)

  end subroutine scale_two_double_complex

end module bml_scale_m
