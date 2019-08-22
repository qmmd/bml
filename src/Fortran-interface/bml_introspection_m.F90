!> Introspection procedures.
module bml_introspection_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  public :: bml_get_N
  public :: bml_get_M
  public :: bml_get_type
  public :: bml_get_element_precision
  public :: bml_get_precision
  public :: bml_get_element_type
  public :: bml_get_row_bandwidth
  public :: bml_get_bandwidth
  public :: bml_get_distribution_mode
  public :: bml_get_sparsity

contains

  !> Return the matrix size.
  !!
  !!\param a The matrix.
  !!\return The matrix size.
  function bml_get_n(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_get_n

    bml_get_n = bml_get_N_C(a%ptr)

  end function bml_get_n

  !> Return the max non-zero elements per row.
  !!
  !!\param a The matrix.
  !!\return The max non-zeroes per row.
  function bml_get_m(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_get_m

    bml_get_m = bml_get_M_C(a%ptr)

  end function bml_get_m

  !> Get the bandwidth of non-zero elements in a given row.
  !!
  !! @param a The matrix.
  !! @param i The row.
  !! @returns The bandwidth of non-zero elements (bandwidth) on that row.
  function bml_get_row_bandwidth(a, i)

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), intent(in) :: i
    integer :: bml_get_row_bandwidth

    bml_get_row_bandwidth = bml_get_row_bandwidth_C(a%ptr, i-1)

  end function bml_get_row_bandwidth

  !> Get the bandwidth of non-zero elements of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bandwidth of non-zero elements (bandwidth) of the matrix.
  function bml_get_bandwidth(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_get_bandwidth

    bml_get_bandwidth = bml_get_bandwidth_C(a%ptr)

  end function bml_get_bandwidth

  !> Get the type of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bml format type of the matrix.
  function bml_get_type(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_get_type_num
    character(20) :: bml_get_type

    bml_get_type_num = bml_get_type_C(a%ptr)

    select case(bml_get_type_num)
    case(0)
      bml_get_type = "Unformatted"
    case(1)
      bml_get_type = "dense"
    case(2)
      bml_get_type = "ellpack"
    case default
      stop 'Unknown matrix type in bml_get_type'
    end select

  end function bml_get_type

  !> Get the precision/type index of the elements of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bml type and precision index for the matrix elements.
  !! 0 = not initialized
  !! 1 = single real
  !! 2 = double real
  !! 3 = single complex
  !! 4 = double complex
  function bml_get_precision(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_get_precision

    bml_get_precision = bml_get_precision_C(a%ptr)

  end function bml_get_precision

  !> Get the precision of the elements of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bml format precision the matrix elements.
  function bml_get_element_precision(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_precision_id
    integer :: bml_get_element_precision

    bml_precision_id = bml_get_precision_C(a%ptr)

    select case(bml_precision_id)
    case(0)
      stop 'Type/precision elements not initialized'
    case(1)
      bml_get_element_precision = kind(1.0)
    case(2)
      bml_get_element_precision = kind(1.0d0)
    case default
      stop 'Unknown elements type/precision'
    end select

  end function bml_get_element_precision


  !> Get the elements type of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bml format type of all the elements of the matrix.
  function bml_get_element_type(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: bml_precision_id
    character(20) :: bml_get_element_type

    bml_precision_id = bml_get_precision_C(a%ptr)

    select case(bml_precision_id)
    case(0)
      stop 'Type/precision elements not initialized'
    case(1)
      bml_get_element_type = "real"
    case(2)
      bml_get_element_type = "real"
    case(3)
      bml_get_element_type = "complex"
    case(4)
      bml_get_element_type = "complex"
    case default
      stop 'Unknown elements type/precision'
    end select

  end function bml_get_element_type

  !> Get the distribution mode of a matrix.
  !!
  !! @param a The matrix.
  !! @returns The bml distribution mode of the matrix.
  function bml_get_distribution_mode(a)

    type(bml_matrix_t), intent(in) :: a
    integer :: dmode
    character(20):: bml_get_distribution_mode

    dmode = bml_get_distribution_mode_C(a%ptr)

    select case(dmode)
    case(0)
      bml_get_distribution_mode = "sequential"
    case(1)
      bml_get_distribution_mode = "distributed"
    case(3)
      bml_get_distribution_mode = "graph_distributed"
    case default
      stop 'Unknown distribution type in bml_get_distribution_mode'
    end select

  end function bml_get_distribution_mode


  !> Get the sparsity of a bml matrix.
  !!
  !! @param a The matrix.
  !! @params threshold The threshold parameter.
  !! @returns The sparsity of the matrix.
  function bml_get_sparsity(a, threshold) result(sparsity)

    type(bml_matrix_t), intent(in) :: a
    real(C_DOUBLE), intent(in) :: threshold
    real(C_DOUBLE) :: sparsity

    sparsity = bml_get_sparsity_C(a%ptr, threshold)

  end function bml_get_sparsity

end module bml_introspection_m
