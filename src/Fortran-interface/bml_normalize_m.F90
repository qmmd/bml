!> Normalize and calculate gershgorin bounds for a matrix
module bml_normalize_m

  implicit none

  private

  interface bml_normalize
     module procedure bml_normalize
  end interface bml_normalize

  interface bml_gershgorin
     module procedure bml_gershgorin
  end interface bml_gershgorin

  public :: bml_normalize
  public :: bml_gershgorin
  public :: bml_gershgorin_partial

contains

  !> Normalize matrix given gershgorin bounds.
  !!
  !! \ingroup normalize_group_F
  !!
  !! \param a Matrix
  !! \param mineval Calculated min
  !! \param maxeval Calculated max
  subroutine bml_normalize(a, mineval, maxeval)

    use bml_types_m
    use bml_c_interface_m

    type(bml_matrix_t), intent(inout) :: a
    double precision, intent(in) :: mineval, maxeval

    call bml_normalize_C(a%ptr, mineval, maxeval)

  end subroutine bml_normalize

  !> Calculate gershgorin bounds for a matrix
  !!
  !! \ingroup normalize_group_F
  !!
  !! \param a Matrix
  !! \param a_gbnd Calculated min and max
  subroutine bml_gershgorin(a, a_gbnd)

    use bml_c_interface_m
    use bml_types_m
    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    double precision, allocatable, intent(inout) :: a_gbnd(:)
    type(C_PTR) :: ag_ptr
    real(C_DOUBLE), pointer :: a_gbnd_ptr(:)

    ag_ptr = bml_gershgorin_C(a%ptr)
    call c_f_pointer(ag_ptr, a_gbnd_ptr, [2])
    a_gbnd = a_gbnd_ptr
    call bml_free(ag_ptr)

  end subroutine bml_gershgorin

  !> Calculate gershgorin bounds for a partial matrix
  !!
  !! \ingroup normalize_group_F
  !!
  !! \param a Matrix
  !! \param nrows Number of rows used
  !! \param a_gbnd Calculated min and max
  subroutine bml_gershgorin_partial(a, a_gbnd, nrows)

    use bml_c_interface_m
    use bml_types_m
    use bml_allocate_m

    type(bml_matrix_t), intent(in) :: a
    integer, intent(in) :: nrows

    double precision, allocatable, intent(inout) :: a_gbnd(:)

    type(C_PTR) :: ag_ptr
    double precision, pointer :: a_gbnd_ptr(:)

    ag_ptr = bml_gershgorin_partial_C(a%ptr, nrows)
    call c_f_pointer(ag_ptr, a_gbnd_ptr, [2])
    a_gbnd = a_gbnd_ptr
    call bml_free(ag_ptr)

  end subroutine bml_gershgorin_partial

end module bml_normalize_m
