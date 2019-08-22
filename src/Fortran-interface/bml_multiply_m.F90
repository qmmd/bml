!> Matrix multiplication.
module bml_multiply_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  public :: bml_multiply
  public :: bml_multiply_x2

contains

  !> Multiply two matrices.
  !!
  !! \ingroup multiply_group
  !!
  !! \f$ C \leftarrow \alpha A \times B + \beta C \f$
  !!
  !! The optional scaling factors \f$ \alpha \f$ and \f$ \beta \f$
  !! default to \f$ \alpha = 1 \f$ and \f$ \beta = 0 \f$.
  !!
  !! \param a Matrix \f$ A \f$.
  !! \param b Matrix \f$ B \f$.
  !! \param c Matrix \f$ C \f$.
  !! \param alpha The factor \f$ \alpha \f$.
  !! \param beta The factor \f$ \beta \f$.
  !! \param threshold The threshold \f$ threshold \f$.
  subroutine bml_multiply(a, b, c, alpha, beta, threshold)

    type(bml_matrix_t), intent(in) :: a, b
    type(bml_matrix_t), intent(inout) :: c
    real(C_DOUBLE), optional, intent(in) :: alpha
    real(C_DOUBLE), optional, intent(in) :: beta
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: alpha_
    real(C_DOUBLE) :: beta_
    real(C_DOUBLE) :: threshold_

    if(present(alpha)) then
      alpha_ = alpha
    else
      alpha_ = 1
    end if
    if(present(beta)) then
      beta_ = beta
    else
      beta_ = 0
    end if
    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0
    end if
    call bml_multiply_c(a%ptr, b%ptr, c%ptr, alpha_, beta_, threshold_)

  end subroutine bml_multiply

  !> Square a matrix.
  !!
  !! \ingroup multiply_group
  !!
  !! \f$B \leftarrow A \times A \f$
  !!
  !! \param a Matrix \f$ A \f$.
  !! \param b Matrix \f$ B \f$.
  !! \param threshold The threshold \f$ threshold \f$.
  subroutine bml_multiply_x2(a, b, threshold, trace)

    use bml_types_m
    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: b
    double precision, optional, intent(in) :: threshold
    double precision, allocatable, intent(inout) :: trace(:)

    double precision :: threshold_

    type(C_PTR) :: trace_ptr
    double precision, pointer :: a_trace_ptr(:)

    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = 0
    end if

    trace_ptr = bml_multiply_x2_c(a%ptr, b%ptr, threshold_)
    call c_f_pointer(trace_ptr, a_trace_ptr, [2])
    trace = a_trace_ptr
    call bml_free(trace_ptr)

  end subroutine bml_multiply_x2

end module bml_multiply_m
