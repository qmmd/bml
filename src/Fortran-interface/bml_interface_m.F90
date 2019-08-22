!> Interface module.
module bml_interface_m

  use bml_c_interface_m
  use bml_types_m

  implicit none
  private

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix type is unitialized.
  integer, parameter :: bml_matrix_type_uninitialized_enum_id = 0

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix type is dense.
  integer, parameter :: bml_matrix_type_dense_enum_id = 1

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix type is ellpack.
  integer, parameter :: bml_matrix_type_ellpack_enum_id = 2

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix type is ellblock.
  integer, parameter :: bml_matrix_type_ellblock_enum_id = 3

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix precision is unitialized.
  integer, parameter :: bml_matrix_precision_uninitialized_id = 0

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix precision is single real.
  integer, parameter :: bml_matrix_precision_single_real_enum_id = 1

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix precision is double real.
  integer, parameter :: bml_matrix_precision_double_real_enum_id = 2

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix precision is single complex.
  integer, parameter :: bml_matrix_precision_single_complex_enum_id = 3

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Matrix precision is double complex.
  integer, parameter :: bml_matrix_precision_double_complex_enum_id = 4

  !> The dense matrix element order.
  integer, parameter :: BML_DENSE_COLUMN_MAJOR = 1

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Distribution mode is sequential.
  integer, parameter :: bml_distribution_mode_sequential_enum_id = 0

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Distribution mode is distributed.
  integer, parameter :: bml_distribution_mode_distributed_enum_id = 1

  !> The enum values of the C API. Keep this synchronized with the
  !! enum in bml_types.h.
  !!
  !! Distribution mode is graph distributed.
  integer, parameter :: bml_distribution_mode_graph_distributed_enum_id = 2

  public :: get_matrix_id, get_element_id, get_dmode_id
  public :: BML_DENSE_COLUMN_MAJOR

contains

  !> Convert the matrix type and precisions strings into enum values.
  !!
  !! @param type_string The string used in the Fortran API to identify
  !! the matrix type and precision.
  !! @return The corresponding integer value matching the enum values
  !! in bml_matrix_types_t and bml_matrix_precision_t.
  function get_matrix_id(type_string) result(id)

    character(len=*), intent(in) :: type_string
    integer(C_INT) :: id
    integer dummy

    select case (trim(type_string))
    case(BML_MATRIX_DENSE)
      id = bml_matrix_type_dense_enum_id
    case(BML_MATRIX_ELLPACK)
      id = bml_matrix_type_ellpack_enum_id
    case(BML_MATRIX_ELLBLOCK)
      id = bml_matrix_type_ellblock_enum_id
    case default
      print *, "unknown matrix type"//trim(type_string)
      error stop
    end select

  end function get_matrix_id

  function get_element_id(element_type, element_kind) result(id)

    character(len=*), intent(in) :: element_type
    integer, intent(in) :: element_kind
    integer(C_INT) :: id

    select case (trim(element_type))

    case (BML_ELEMENT_REAL)
      select case (element_kind)
      case (C_FLOAT)
        id = bml_matrix_precision_single_real_enum_id
      case (C_DOUBLE)
        id = bml_matrix_precision_double_real_enum_id
      case default
        print "(A,1X,I0)", "Unknown element kind:", element_kind
      end select

    case (BML_ELEMENT_COMPLEX)
      select case (element_kind)
      case (C_FLOAT_COMPLEX)
        id = bml_matrix_precision_single_complex_enum_id
      case (C_DOUBLE_COMPLEX)
        id = bml_matrix_precision_double_complex_enum_id
      case default
        print "(A,1X,I0)", "Unknown element kind:", element_kind
      end select

    case default
      print "(A,1X,A)", "Unknown element type:", element_type

    end select

  end function get_element_id

  !> Convert the distribution mode strings into enum values.
  !!
  !! @param type_string The string used in the Fortran API to identify
  !! the distribution mode.
  !! @return The corresponding integer value matching the enum values
  !! in bml_distribution_mode_t.
  function get_dmode_id(dmode_string) result(id)

    character(len=*), intent(in) :: dmode_string
    integer(C_INT) :: id

    select case (trim(dmode_string))
    case(BML_DMODE_SEQUENTIAL)
      id = bml_distribution_mode_sequential_enum_id
    case(BML_DMODE_DISTRIBUTED)
      id = bml_distribution_mode_distributed_enum_id
    case(BML_DMODE_GRAPH_DISTRIBUTED)
      id = bml_distribution_mode_graph_distributed_enum_id
    case default
      print *, "unknown distribution mode"//trim(dmode_string)
      error stop
    end select

  end function get_dmode_id

end module bml_interface_m
