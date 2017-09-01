module bml_import_m

  use bml_c_interface_m
  use bml_types_m
  use bml_interface_m
  use bml_introspection_m

  implicit none
  private

  interface bml_import_from_dense
     module procedure bml_import_from_dense_single
     module procedure bml_import_from_dense_double
     module procedure bml_import_from_dense_single_complex
     module procedure bml_import_from_dense_double_complex
  end interface bml_import_from_dense

  public :: bml_import_from_dense

contains

  !> Convert a dense matrix into a bml matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param matrix_type The matrix type
  !! \param a_dense The dense matrix
  !! \param a The bml matrix
  !! \param threshold The matrix element magnited threshold
  !! \param m The extra arg
  subroutine bml_import_from_dense_single(matrix_type, a_dense, a, threshold, m, distrib_mode)

    character(len=*), intent(in) :: matrix_type
    real(C_FLOAT), target, intent(in) :: a_dense(:, :)
    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), optional, intent(in) :: threshold
    integer, optional, intent(in) :: m
    character(len=*), optional, intent(in) :: distrib_mode

    integer(C_INT) :: n_
    integer(C_INT) :: m_
    real(C_DOUBLE) :: threshold_
    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    if (present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0
    end if

    n_ = size(a_dense, 1, C_INT)

    if (matrix_type /= BML_MATRIX_DENSE) then
      if (.not. present(m)) then
        write(*, *) "missing parameter m; number of non-zeros per row"
        error stop
      end if
    end if

    if (present(m)) then
      m_ = m
    else
      m_ = n_
    end if

    call bml_deallocate(a)
    associate(a_ptr => a_dense(lbound(a_dense, 1), lbound(a_dense, 2)))
      a%ptr = bml_import_from_dense_C(get_matrix_id(matrix_type), &
           & get_element_id(BML_ELEMENT_REAL, C_FLOAT), &
           & BML_DENSE_COLUMN_MAJOR, n_, m_, &
           & c_loc(a_ptr), threshold_, get_dmode_id(distrib_mode_))
    end associate

  end subroutine bml_import_from_dense_single

  !> Convert a dense matrix into a bml matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param matrix_type The matrix type
  !! \param a_dense The dense matrix
  !! \param a The bml matrix
  !! \param threshold The matrix element magnited threshold
  !! \param m the extra arg
  subroutine bml_import_from_dense_double(matrix_type, a_dense, a, threshold, m, distrib_mode)

    character(len=*), intent(in) :: matrix_type
    real(C_DOUBLE), target, intent(in) :: a_dense(:, :)
    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), optional, intent(in) :: threshold
    integer, optional, intent(in) :: m
    character(len=*), optional, intent(in) :: distrib_mode

    integer(C_INT) :: n_
    integer(C_INT) :: m_
    real(C_DOUBLE) :: threshold_
    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    if (present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0
    end if

    n_ = size(a_dense, 1, C_INT)

    if (matrix_type /= BML_MATRIX_DENSE) then
      if (.not. present(m)) then
        write(*, *) "missing parameter m; number of non-zeros per row"
        error stop
      end if
    end if

    if (present(m)) then
      m_ = m
    else
      m_ = n_
    end if

    call bml_deallocate(a)
    associate(a_ptr => a_dense(lbound(a_dense, 1), lbound(a_dense, 2)))
      a%ptr = bml_import_from_dense_C(get_matrix_id(matrix_type), &
           & get_element_id(BML_ELEMENT_REAL, C_DOUBLE), &
           & BML_DENSE_COLUMN_MAJOR, n_, m_, &
           & c_loc(a_ptr), threshold_, get_dmode_id(distrib_mode_))
    end associate

  end subroutine bml_import_from_dense_double

  !> Convert a dense matrix into a bml matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param matrix_type The matrix type
  !! \param a_dense The dense matrix
  !! \param a The bml matrix
  !! \param threshold The matrix element magnited threshold
  !! \param m The extra arg
  subroutine bml_import_from_dense_single_complex(matrix_type, a_dense, a, threshold, m, distrib_mode)

    character(len=*), intent(in) :: matrix_type
    complex(C_FLOAT_COMPLEX), target, intent(in) :: a_dense(:, :)
    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), optional, intent(in) :: threshold
    integer, optional, intent(in) :: m
    character(len=*), optional, intent(in) :: distrib_mode

    integer(C_INT) :: n_
    integer(C_INT) :: m_
    real(C_DOUBLE) :: threshold_
    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    if (present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0
    end if

    n_ = size(a_dense, 1, C_INT)

    if (matrix_type /= BML_MATRIX_DENSE) then
      if (.not. present(m)) then
        write(*, *) "missing parameter m; number of non-zeros per row"
        error stop
      end if
    end if

    if (present(m)) then
      m_ = m
    else
      m_ = n_
    end if

    call bml_deallocate(a)
    associate(a_ptr => a_dense(lbound(a_dense, 1), lbound(a_dense, 2)))
      a%ptr = bml_import_from_dense_C(get_matrix_id(matrix_type), &
           & get_element_id(BML_ELEMENT_COMPLEX, C_FLOAT_COMPLEX), &
           & BML_DENSE_COLUMN_MAJOR, n_, m_, &
           & c_loc(a_ptr), threshold_, get_dmode_id(distrib_mode_))
    end associate

  end subroutine bml_import_from_dense_single_complex

  !> Convert a dense matrix into a bml matrix.
  !!
  !! \ingroup convert_group_Fortran
  !!
  !! \param matrix_type The matrix type
  !! \param a_dense The dense matrix
  !! \param a The bml matrix
  !! \param threshold The matrix element magnited threshold
  !! \param m the extra arg
  subroutine bml_import_from_dense_double_complex(matrix_type, a_dense, a, threshold, m, distrib_mode)

    character(len=*), intent(in) :: matrix_type
    complex(C_DOUBLE_COMPLEX), target, intent(in) :: a_dense(:, :)
    type(bml_matrix_t), intent(inout) :: a
    real(C_DOUBLE), optional, intent(in) :: threshold
    integer, optional, intent(in) :: m
    character(len=*), optional, intent(in) :: distrib_mode

    integer(C_INT) :: n_
    integer(C_INT) :: m_
    real(C_DOUBLE) :: threshold_
    character(len=20) :: distrib_mode_

    if (present(distrib_mode)) then
      distrib_mode_ = distrib_mode
    else
      distrib_mode_ = bml_dmode_sequential
    endif

    if (present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0
    end if

    n_ = size(a_dense, 1, C_INT)

    if (matrix_type /= BML_MATRIX_DENSE) then
      if (.not. present(m)) then
        write(*, *) "missing parameter m; number of non-zeros per row"
        error stop
      end if
    end if

    if (present(m)) then
      m_ = m
    else
      m_ = n_
    end if

    call bml_deallocate(a)
    associate(a_ptr => a_dense(lbound(a_dense, 1), lbound(a_dense, 2)))
      a%ptr = bml_import_from_dense_C(get_matrix_id(matrix_type), &
           & get_element_id(BML_ELEMENT_COMPLEX, C_DOUBLE_COMPLEX), &
           & BML_DENSE_COLUMN_MAJOR, n_, m_, &
           & c_loc(a_ptr), threshold_, get_dmode_id(distrib_mode_))
    end associate

  end subroutine bml_import_from_dense_double_complex

end module bml_import_m
