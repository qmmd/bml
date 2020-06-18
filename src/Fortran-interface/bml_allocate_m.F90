!> Matrix allocation functions.
module bml_allocate_m

  use bml_c_interface_m
  use bml_types_m
  use bml_interface_m

  implicit none
  private

  public :: bml_allocated
  public :: bml_banded_matrix
  public :: bml_clear
  public :: bml_free
  public :: bml_identity_matrix
  public :: bml_noinit_matrix
  public :: bml_random_matrix
  public :: bml_update_domain
  public :: bml_zero_matrix
  public :: bml_block_matrix

contains

  !> Clear a matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param a The matrix.
  subroutine bml_clear(a)

    type(bml_matrix_t), intent(inout) :: a

    call bml_clear_C(a%ptr)

  end subroutine bml_clear

  !> Free a c_ptr.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param cptr The C pointer.
  subroutine bml_free(cptr)

    type(C_PTR), intent(inout) :: cptr

    call bml_free_C(cptr)

  end subroutine bml_free

  !> Check if matrix is allocated.
  !!
  !!\param a The matrix.
  !!\return true/false
  function bml_allocated(a)

    type(bml_matrix_t), intent(in) :: a
    logical :: bml_allocated
    integer :: aflag

    aflag = bml_allocated_C(a%ptr)

    if (aflag .gt. 0) then
      bml_allocated = .true.
    else
      bml_allocated = .false.
    end if

  end function bml_allocated

  !> Create the zero matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param matrix_type The matrix type.
  !! \param element_type  Element type of the matrix.
  !! \param element_precision The precision of the matrix elements.
  !! \param n The matrix size.
  !! \param m The extra arg.
  !! \param a The matrix.
  !! \param distrib_mode The matrix distribution mode.
  subroutine bml_zero_matrix(matrix_type, element_type, element_precision, &
       & n, m, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: n, m
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_zero_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & n, m, get_dmode_id(distrib_mode_))

  end subroutine bml_zero_matrix

  subroutine bml_block_matrix(matrix_type, element_type, element_precision, &
       & nb, mb, m, bsizes, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: nb, mb, m
    integer, allocatable, intent(in) :: bsizes(:)
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_block_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & nb, mb, m, bsizes, get_dmode_id(distrib_mode_))

  end subroutine bml_block_matrix

  !> Create a matrix without initializing.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param matrix_type The matrix type.
  !! \param element_type  Element type of the matrix.
  !! \param element_precision The precision of the matrix elements.
  !! \param n The matrix size.
  !! \param m The extra arg.
  !! \param a The matrix.
  !! \param distrib_mode The matrix distribution mode.
  subroutine bml_noinit_matrix(matrix_type, element_type, element_precision, &
       & n, m, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: n, m
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_noinit_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & n, m, get_dmode_id(distrib_mode_))

  end subroutine bml_noinit_matrix

  !> Create a banded matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param matrix_type The matrix type.
  !! \param element_type  Element type of the matrix.
  !! \param element_precision The precision of the matrix.
  !! \param n The matrix size.
  !! \param m The extra arg.
  !! \param a The matrix.
  !! \param distrib_mode The matrix distribution mode.
  subroutine bml_banded_matrix(matrix_type, element_type, element_precision, &
       & n, m, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: n, m
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_banded_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & n, m, get_dmode_id(distrib_mode_))

  end subroutine bml_banded_matrix

  !> Create a random matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param matrix_type The matrix type.
  !! \param element_type  Element type of the matrix.
  !! \param element_precision The precision of the matrix.
  !! \param n The matrix size.
  !! \param m The extra arg.
  !! \param a The matrix.
  !! \param distrib_mode The matrix distribution mode.
  subroutine bml_random_matrix(matrix_type, element_type, element_precision, &
       & n, m, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: n, m
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_random_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & n, m, get_dmode_id(distrib_mode_))

  end subroutine bml_random_matrix

  !> Create the identity matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param matrix_type The matrix type.
  !! \param element_type  Element type of the matrix.
  !! \param element_precision The precision of the matrix.
  !! \param n The matrix size.
  !! \param m The extra arg.
  !! \param a The matrix.
  !! \param distrib_mode The matrix distribution mode.
  subroutine bml_identity_matrix(matrix_type, element_type, element_precision, &
       & n, m, a, distrib_mode)

    character(len=*), intent(in) :: matrix_type, element_type
    character(len=*), optional, intent(in) :: distrib_mode
    integer, intent(in) :: element_precision
    integer(C_INT), intent(in) :: n, m
    type(bml_matrix_t), intent(inout) :: a

    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    call bml_deallocate(a)
    a%ptr = bml_identity_matrix_C(get_matrix_id(matrix_type), &
         & get_element_id(element_type, element_precision), &
         & n, m, get_dmode_id(distrib_mode_))

  end subroutine bml_identity_matrix

  !> Update domain of a matrix.
  !!
  !! \ingroup allocate_group_Fortran
  !!
  !! \param a The matrix.
  !! \param globalPartMin First part on each rank
  !! \param globalPartMax Last part on each rank
  !! \param nnodesInPart Number of nodes in each part
  subroutine bml_update_domain(a, globalPartMin, globalPartMax, nnodesInPart)

    integer(C_INT), target, intent(in) :: globalPartMin(*)
    integer(C_INT), target, intent(in) :: globalPartMax(*)
    integer(C_INT), target, intent(in) :: nnodesInPart(*)
    type(bml_matrix_t), intent(inout) :: a

    call bml_update_domain_C(a%ptr, c_loc(globalPartMin), c_loc(globalPartMax), &
         c_loc(nnodesInPart))

  end subroutine bml_update_domain

end module bml_allocate_m
