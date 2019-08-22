!> Matrix trace.
module bml_trace_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  public :: bml_trace
  public :: bml_trace_mult

contains

  !> Calculate the trace of a matrix.
  !!
  !! \f$ \leftarrow \mathrm{Tr} \left[ A \right] \f$
  !!
  !! \param a The matrix.
  function bml_trace(a) result(tr_a)

    class(bml_matrix_t), intent(in) :: a
    real(C_DOUBLE) :: tr_a

    tr_a = bml_trace_C(a%ptr)

  end function bml_trace

  !> Calculate the trace of a matrix multiplication.
  !!
  !! \param a The matrix a.
  !! \param b The matrix b.
  function bml_trace_mult(a, b) result(tr_mult)

    class(bml_matrix_t), intent(in) :: a, b
    real(C_DOUBLE) :: tr_mult

    tr_mult= bml_trace_mult_C(a%ptr, b%ptr)

  end function bml_trace_mult

end module bml_trace_m
