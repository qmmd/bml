module bml_add_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  !> \addtogroup add_group_Fortran
  !! @{

  !> Add two matrices.
  interface bml_add
     module procedure add_two
  end interface bml_add

  interface bml_add_deprecated
     module procedure add_two_deprecated
  end interface bml_add_deprecated

  !> Add identity matrix to a matrix.
  interface bml_add_identity
     module procedure add_identity_one
  end interface bml_add_identity
  !> @}

  public :: bml_add
  public :: bml_add_deprecated
  public :: bml_add_norm
  public :: bml_add_identity
  public :: bml_scale_add_identity

contains

  !> Add two matrices.
  !!
  !! \f$ A \leftarrow \alpha A + \beta B \f$
  !!
  !! The optional scalars \f$ \alpha \f$ and \f$ \beta \f$ default to
  !! 1.
  !!
  !! \param a Matrix \f$ A \f$
  !! \param b Matrix \f$ B \f$
  !! \param alpha Factor \f$ \alpha \f$
  !! \param beta Factor \f$ \beta \f$
  !! \param threshold \f$ threshold \f$
  subroutine add_two(a, b, alpha, beta, threshold)

    type(bml_matrix_t), intent(inout) :: a
    type(bml_matrix_t), intent(in) :: b
    real(C_DOUBLE), intent(in) :: alpha
    real(C_DOUBLE), intent(in) :: beta
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0.0_C_DOUBLE
    end if
    call bml_add_C(a%ptr, b%ptr, alpha, beta, threshold_)

  end subroutine add_two

  !> \deprecated Add two matrices.
  !!
  !! \f$ A \leftarrow \alpha A + \beta B \f$
  !!
  !! The optional scalars \f$ \alpha \f$ and \f$ \beta \f$ default to
  !! 1.
  !!
  !! \param alpha Factor \f$ \alpha \f$
  !! \param a Matrix \f$ A \f$
  !! \param beta Factor \f$ \beta \f$
  !! \param b Matrix \f$ B \f$
  !! \param threshold \f$ threshold \f$
  subroutine add_two_deprecated(alpha, a, beta, b, threshold)

    real(C_DOUBLE), intent(in) :: alpha
    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), intent(in) :: beta
    type(bml_matrix_t), intent(in) :: b
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0.0_C_DOUBLE
    end if
    call add_two(a, b, alpha, beta, threshold_)

  end subroutine add_two_deprecated

  !> Add two matrices and calculate trnorm.
  !!
  !! \f$ A \leftarrow \alpha A + \beta B \f$
  !!
  !! The optional scalars \f$ \alpha \f$ and \f$ \beta \f$ default to
  !! 1.
  !!
  !! \param a Matrix \f$ A \f$
  !! \param b Matrix \f$ B \f$
  !! \param alpha Factor \f$ \alpha \f$
  !! \param beta Factor \f$ \beta \f$
  !! \param threshold \f$ threshold \f$
  function bml_add_norm(a, b, alpha, beta, threshold) result(trnorm)

    real(C_DOUBLE), intent(in) :: alpha
    class(bml_matrix_t), intent(in) :: a
    real(C_DOUBLE), intent(in) :: beta
    class(bml_matrix_t), intent(in) :: b
    real(C_DOUBLE), optional, intent(in) :: threshold
    real(C_DOUBLE) :: trnorm

    real(C_DOUBLE) :: threshold_

    if(present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0.0_C_DOUBLE
    end if

    trnorm = bml_add_norm_C(a%ptr, b%ptr, alpha, beta, threshold_)

  end function bml_add_norm

  !> Add a scaled identity matrix to a bml matrix.
  !!
  !! \f$ A \leftarrow A + \beta \mathrm{Id} \f$
  !!
  !! The optional scalars \f$ \alpha \f$ and \f$ \beta \f$ default to
  !! 1.
  !!
  !! \param a Matrix A
  !! \param beta Factor \f$ \alpha \f$
  !! \param threshold \f$ threshold \f$
  subroutine add_identity_one(a, alpha, threshold)

    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), intent(in) :: alpha
    real(C_DOUBLE), intent(in) :: threshold

    call bml_add_identity_C(a%ptr, alpha, threshold)

  end subroutine add_identity_one

  !> Add a scaled identity matrix to a scaled a bml matrix.
  !!
  !! \f$ A \leftarrow \alpha A + \beta \mathrm{Id} \f$
  !!
  !! The optional scalars \f$ \alpha \f$ and \f$ \beta \f$ default to
  !! 1.
  !!
  !! \param a Matrix A
  !! \param alpha Factor \f$ \alpha \f$
  !! \param beta Factor \f$ \beta \f$
  !! \param threshold \f$ threshold \f$
  subroutine bml_scale_add_identity(a, alpha, beta, threshold)

    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), intent(in) :: alpha
    real(C_DOUBLE), intent(in) :: beta
    real(C_DOUBLE), intent(in) :: threshold

    call bml_scale_add_identity_C(a%ptr, alpha, beta, threshold)

  end subroutine bml_scale_add_identity

end module bml_add_m
