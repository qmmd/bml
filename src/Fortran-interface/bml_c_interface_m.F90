!> Contains all Fortran interface defintions to the C-API of the BML-library.
!!
!! \note This module also exports the content of the iso_c_binding intrinsic
!! module, so for all other modules in the Fortran interface, it should be
!! enough to import this module when interacting with the C-API.
!!
module bml_c_interface_m

  use, intrinsic :: iso_c_binding

  implicit none

  ! Note: According to Sec. 15.3.7.2.6: "any dummy argument without
  ! the value attribute corresponds to a formal parameter of the
  ! prototype that is of a pointer type, and the dummy argument is
  ! interoperable with an entity of the referenced type (ISO/IEC
  ! 9899:1999, 6.2.5, 7.17, and 7.18.1) of the formal parameter, ..."
  !
  ! In other words, a type(C_PTR) dummy argument is interoperable with
  ! the void** type.

  interface

     subroutine bml_add_C(a, b, alpha, beta, threshold) bind(C, name="bml_add")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      real(C_DOUBLE), value, intent(in) :: alpha
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
    end subroutine bml_add_C

    function bml_add_norm_C(a, b, alpha, beta, threshold) &
        & bind(C, name="bml_add_norm")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      real(C_DOUBLE), value, intent(in) :: alpha
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
      real(C_DOUBLE) :: bml_add_norm_C
    end function bml_add_norm_C

    function bml_allocated_C(a) &
        & bind(C, name="bml_allocated")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_allocated_C
    end function bml_allocated_C

    subroutine bml_add_identity_C(a, beta, threshold) &
        & bind(C, name="bml_add_identity")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value :: a
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
    end subroutine bml_add_identity_C

    subroutine bml_scale_add_identity_C(a, alpha, beta, threshold) &
        & bind(C, name="bml_scale_add_identity")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value :: a
      real(C_DOUBLE), value, intent(in) :: alpha
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
    end subroutine bml_scale_add_identity_C

    subroutine bml_allGatherVParallel_C(a) &
        & bind(C, name="bml_allGatherVParallel")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
    end subroutine bml_allGatherVParallel_C

    subroutine bml_adjungate_triangle_C(a, triangle) &
        & bind(C, name="bml_adjungate_triangle")
      import :: C_PTR, C_CHAR
      type(C_PTR), value :: a
      character(C_CHAR), intent(in) :: triangle(*)
    end subroutine bml_adjungate_triangle_C

    subroutine bml_transpose_triangle_C(a, triangle) &
        & bind(C, name="bml_transpose_triangle")
      import :: C_PTR, C_CHAR
      type(C_PTR), value :: a
      character(C_CHAR), value, intent(in) :: triangle
    end subroutine bml_transpose_triangle_C

    function bml_banded_matrix_C(matrix_type, matrix_precision, n, m, &
        & distrib_mode) bind(C, name="bml_banded_matrix")
      import :: C_INT, C_PTR
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: m
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_banded_matrix_C
    end function bml_banded_matrix_C

    function bml_import_from_dense_C(matrix_type, matrix_precision, order, &
        & n, m, a, threshold, distrib_mode) &
        bind(C, name="bml_import_from_dense")
      import :: C_INT, C_PTR, C_DOUBLE
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: order
      integer(C_INT), value, intent(in) :: n, m
      type(C_PTR), value, intent(in) :: a
      real(C_DOUBLE), value, intent(in) :: threshold
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_import_from_dense_C
    end function bml_import_from_dense_C

    function bml_export_to_dense_C(a, order) &
        & bind(C, name="bml_export_to_dense")
      import :: C_INT, C_PTR
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: order
      type(C_PTR) :: bml_export_to_dense_C
    end function bml_export_to_dense_C

    subroutine bml_copy_C(a, b) bind(C, name="bml_copy")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
    end subroutine bml_copy_C

    function bml_copy_new_C(a) bind(C, name="bml_copy_new")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR) :: bml_copy_new_C
    end function bml_copy_new_C

    subroutine bml_reorder_C(a, perm) bind(C, name="bml_reorder")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: perm
    end subroutine bml_reorder_C

    subroutine bml_save_domain_C(a) bind(C, name="bml_save_domain")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
    end subroutine bml_save_domain_C

    subroutine bml_restore_domain_C(a) bind(C, name="bml_restore_domain")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
    end subroutine bml_restore_domain_C

    subroutine bml_update_domain_C(a, globalPartMin, globalPartMax, &
        nnodesInPart) bind(C, name="bml_update_domain")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: globalPartMin
      type(C_PTR), value, intent(in) :: globalPartMax
      type(C_PTR), value, intent(in) :: nnodesInPart
    end subroutine bml_update_domain_C

    subroutine bml_diagonalize_C(a, eigenvalues, eigenvectors) &
        & bind(C, name="bml_diagonalize")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value :: eigenvalues
      type(C_PTR), value :: eigenvectors
    end subroutine bml_diagonalize_C

    subroutine bml_deallocate_C(a) bind(C, name="bml_deallocate")
      import :: C_PTR
      type(C_PTR) :: a
    end subroutine bml_deallocate_C

    subroutine bml_clear_C(a) bind(C, name="bml_clear")
      import :: C_PTR
      type(C_PTR) :: a
    end subroutine bml_clear_C

    function bml_get_C(a, i, j) bind(C, name="bml_get")
      import :: C_PTR, C_INT, C_FLOAT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      integer(C_INT), value, intent(in) :: j
      type(C_PTR) :: bml_get_C
    end function bml_get_C

    function bml_get_N_C(a) bind(C, name="bml_get_N")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_N_C
    end function bml_get_N_C

    function bml_get_M_C(a) bind(C, name="bml_get_M")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_M_C
    end function bml_get_M_C

    function bml_get_type_C(a) bind(C, name="bml_get_type")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_type_C
    end function bml_get_type_C

    function bml_get_precision_C(a) bind(C, name="bml_get_precision")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_precision_C
    end function bml_get_precision_C

    function bml_get_distribution_mode_C(a) &
        & bind(C, name="bml_get_distribution_mode")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_distribution_mode_C
    end function bml_get_distribution_mode_C

    function bml_get_row_C(a, i) bind(C, name="bml_get_row")
      import :: C_PTR, C_INT, C_double
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      type(C_PTR) :: bml_get_row_C
    end function bml_get_row_C

    subroutine bml_get_diagonal_C(a, diagonal) bind(C, name="bml_get_diagonal")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: diagonal
    end subroutine bml_get_diagonal_C

    subroutine bml_initF_C(fcomm) bind(C, name="bml_initF")
      import :: C_PTR, C_INT
      integer(C_INT), value, intent(in) :: fcomm
    end subroutine bml_initF_C

    subroutine bml_shutdownF_C() bind(C, name="bml_shutdownF")
    end subroutine bml_shutdownF_C

    function bml_sum_squares_C(a) bind(C, name="bml_sum_squares")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      real(C_DOUBLE) :: bml_sum_squares_C
    end function bml_sum_squares_C

    function bml_sum_squares_submatrix_C(a, core_size) &
      bind(C, name="bml_sum_squares_submatrix")
      import :: C_PTR, C_DOUBLE, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: core_size
      real(C_DOUBLE) :: bml_sum_squares_submatrix_C
    end function bml_sum_squares_submatrix_C

    function bml_sum_squares2_C(a, b, alpha, beta, threshold) &
      bind(C, name="bml_sum_squares2")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      real(C_DOUBLE), value, intent(in) :: alpha
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
      real(C_DOUBLE) :: bml_sum_squares2_C
    end function bml_sum_squares2_C

    function bml_fnorm_C(a) bind(C, name="bml_fnorm")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      real(C_DOUBLE) :: bml_fnorm_C
    end function bml_fnorm_C

    function bml_fnorm2_C(a, b) bind(C, name="bml_fnorm2")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      real(C_DOUBLE) :: bml_fnorm2_C
    end function bml_fnorm2_C

    subroutine bml_normalize_C(a, mineval, maxeval) &
      bind(C, name="bml_normalize")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value :: a
      real(C_DOUBLE), value, intent(in) :: mineval
      real(C_DOUBLE), value, intent(in) :: maxeval
    end subroutine bml_normalize_C

    function bml_gershgorin_C(a) bind(C, name="bml_gershgorin")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR) :: bml_gershgorin_C
    end function bml_gershgorin_C

    function bml_gershgorin_partial_C(a, nrows) bind(C, name="bml_gershgorin_partial")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: nrows
      type(C_PTR) :: bml_gershgorin_partial_C
    end function bml_gershgorin_partial_C

    function bml_get_row_bandwidth_C(a, i) &
        & bind(C, name="bml_get_row_bandwidth")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      integer(C_INT) :: bml_get_row_bandwidth_C
    end function bml_get_row_bandwidth_C

    function bml_get_bandwidth_C(a) bind(C, name="bml_get_bandwidth")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT) :: bml_get_bandwidth_C
    end function bml_get_bandwidth_C

    function bml_get_sparsity_C(a, threshold) bind(C, name="bml_get_sparsity")
      import :: C_PTR, C_DOUBLE, C_INT
      type(C_PTR), value, intent(in) :: a
      real(C_DOUBLE), value, intent(in) :: threshold
      real(C_DOUBLE) :: bml_get_sparsity_C
    end function bml_get_sparsity_C

    function bml_identity_matrix_C(matrix_type, matrix_precision, n, m, &
        & distrib_mode) bind(C, name="bml_identity_matrix")
      import :: C_INT, C_PTR
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: m
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_identity_matrix_C
    end function bml_identity_matrix_C

    function bml_inverse_C(a) bind(C, name="bml_inverse")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR) :: bml_inverse_C
    end function bml_inverse_C

    subroutine bml_multiply_C(a, b, c, alpha, beta, threshold) &
        & bind(C, name="bml_multiply")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      type(C_PTR), value, intent(in) :: c
      real(C_DOUBLE), value, intent(in) :: alpha
      real(C_DOUBLE), value, intent(in) :: beta
      real(C_DOUBLE), value, intent(in) :: threshold
    end subroutine bml_multiply_C

    function bml_multiply_x2_C(x, x2, threshold) &
        & bind(C, name="bml_multiply_x2")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: x
      type(C_PTR), value, intent(in) :: x2
      real(C_DOUBLE), value, intent(in) :: threshold
      type(C_PTR) :: bml_multiply_x2_C
    end function bml_multiply_x2_C

    function bml_random_matrix_C(matrix_type, matrix_precision, n, m, &
        & distrib_mode) bind(C, name="bml_random_matrix")
      import :: C_INT, C_PTR
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: m
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_random_matrix_C
    end function bml_random_matrix_C

    subroutine bml_scale_C(alpha, a, b) bind(C, name="bml_scale")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: alpha
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value :: b
    end subroutine bml_scale_C

    subroutine bml_scale_inplace_C(alpha, a) bind(C, name="bml_scale_inplace")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: alpha
      type(C_PTR), value :: a
    end subroutine bml_scale_inplace_C

    subroutine bml_set_row_C(a, i, row, threshold) bind(C, name="bml_set_row")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      type(C_PTR), value :: row
      real(C_DOUBLE), value :: threshold
    end subroutine bml_set_row_C

    subroutine bml_set_element_C(a, i, j, element) bind(C, name="bml_set_element")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      integer(C_INT), value, intent(in) :: j
      type(C_PTR), value :: element
    end subroutine bml_set_element_C

    subroutine bml_set_element_new_C(a, i, j, element) bind(C, name="bml_set_element_new")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i
      integer(C_INT), value, intent(in) :: j
      type(C_PTR), value :: element
    end subroutine bml_set_element_new_C


    subroutine bml_set_diagonal_C(a, diagonal, threshold) bind(C, name="bml_set_diagonal")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value :: diagonal
      real(C_DOUBLE), value :: threshold
    end subroutine bml_set_diagonal_C

    subroutine bml_matrix2submatrix_index_C(a, b, nodelist, nsize, &
      chlist, vsize, dj_flag) bind(C, name="bml_matrix2submatrix_index")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      type(C_PTR), value, intent(in) :: nodelist
      integer(C_INT), value, intent(in) :: nsize
      type(C_PTR), value, intent(in) :: chlist
      type(C_PTR), value, intent(in) :: vsize
      integer(C_INT), value, intent(in) :: dj_flag
    end subroutine bml_matrix2submatrix_index_C

    subroutine bml_matrix2submatrix_index_graph_C(b, nodelist, nsize, &
      chlist, vsize, dj_flag) bind(C, name="bml_matrix2submatrix_index_graph")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: b
      type(C_PTR), value, intent(in) :: nodelist
      integer(C_INT), value, intent(in) :: nsize
      type(C_PTR), value, intent(in) :: chlist
      type(C_PTR), value, intent(in) :: vsize
      integer(C_INT), value, intent(in) :: dj_flag
    end subroutine bml_matrix2submatrix_index_graph_C

    subroutine bml_matrix2submatrix_C(a, b, chlist, lsize) &
      bind(C, name="bml_matrix2submatrix")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      type(C_PTR), value, intent(in) :: chlist
      integer(C_INT), value, intent(in) :: lsize
    end subroutine bml_matrix2submatrix_C

    subroutine bml_submatrix2matrix_C(a, b, chlist, lsize, llsize, &
      threshold) bind(C, name="bml_submatrix2matrix")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      type(C_PTR), value, intent(in) :: chlist
      integer(C_INT), value, intent(in) :: lsize
      integer(C_INT), value, intent(in) :: llsize
      real(C_DOUBLE), value :: threshold
    end subroutine bml_submatrix2matrix_C

    subroutine bml_adjacency_C(a, xadj, adjncy, base_flag) &
      bind(C, name="bml_adjacency")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: xadj
      type(C_PTR), value, intent(in) :: adjncy
      integer(C_INT), value, intent(in) :: base_flag
    end subroutine bml_adjacency_C

    subroutine bml_adjacency_group_C(a, hindex, nnodes, xadj, adjncy, base_flag) &
      bind(C, name="bml_adjacency_group")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: hindex
      integer(C_INT), value, intent(in) :: nnodes
      type(C_PTR), value, intent(in) :: xadj
      type(C_PTR), value, intent(in) :: adjncy
      integer(C_INT), value, intent(in) :: base_flag
    end subroutine bml_adjacency_group_C

    function bml_group_matrix_C(a, hindex, ngroups, threshold) &
      bind(C, name="bml_group_matrix")
      import :: C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: hindex
      integer(C_INT), value, intent(in) :: ngroups
      real(C_DOUBLE), value, intent(in) :: threshold
      type(C_PTR) :: bml_group_matrix_C
    end function bml_group_matrix_C

    subroutine bml_threshold_C(a, threshold) bind(C, name="bml_threshold")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value :: a
      real(C_DOUBLE), value, intent(in) :: threshold
    end subroutine bml_threshold_C

    subroutine bml_print_bml_vector_C(v, i_l, i_u) &
        & bind(C, name="bml_print_bml_vector")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: v
      integer(C_INT), value, intent(in) :: i_l
      integer(C_INT), value, intent(in) :: i_u
    end subroutine bml_print_bml_vector_C

    subroutine bml_print_bml_matrix_C(a, i_l, i_u, j_l, j_u) &
        & bind(C, name="bml_print_bml_matrix")
      import :: C_PTR, C_INT
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i_l
      integer(C_INT), value, intent(in) :: i_u
      integer(C_INT), value, intent(in) :: j_l
      integer(C_INT), value, intent(in) :: j_u
    end subroutine bml_print_bml_matrix_C

    subroutine bml_print_dense_matrix_C(n, matrix_precision, order, a, i_l, &
        & i_u, j_l, j_u) bind(C, name="bml_print_dense_matrix")
      import :: C_PTR, C_INT
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: order
      type(C_PTR), value, intent(in) :: a
      integer(C_INT), value, intent(in) :: i_l
      integer(C_INT), value, intent(in) :: i_u
      integer(C_INT), value, intent(in) :: j_l
      integer(C_INT), value, intent(in) :: j_u
    end subroutine bml_print_dense_matrix_C

    subroutine bml_read_bml_matrix_C(a, filename) &
        & bind(C, name="bml_read_bml_matrix")
      import :: C_PTR, C_CHAR
      type(C_PTR), value, intent(in) :: a
      character(C_CHAR), intent(in) :: filename(*)
    end subroutine bml_read_bml_matrix_C

    function bml_trace_C(a) bind(C, name="bml_trace")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      real(C_DOUBLE) :: bml_trace_C
    end function bml_trace_C

    function bml_traceMult_C(a, b) bind(C, name="bml_traceMult")
      import :: C_PTR, C_DOUBLE
      type(C_PTR), value, intent(in) :: a
      type(C_PTR), value, intent(in) :: b
      real(C_DOUBLE) :: bml_traceMult_C
    end function bml_traceMult_C

    function bml_transpose_new_C(a) bind(C, name="bml_transpose_new")
      import :: C_PTR
      type(C_PTR), value, intent(in) :: a
      type(C_PTR) :: bml_transpose_new_C
    end function bml_transpose_new_C

    subroutine bml_write_bml_matrix_C(a, filename) &
        & bind(C, name="bml_write_bml_matrix")
      import :: C_PTR, C_CHAR
      type(C_PTR), value, intent(in) :: a
      character(C_CHAR), intent(in) :: filename(*)
    end subroutine bml_write_bml_matrix_C

    function bml_zero_matrix_C(matrix_type, matrix_precision, n, m, &
        & distrib_mode) bind(C, name="bml_zero_matrix")
      import :: C_INT, C_PTR
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: m
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_zero_matrix_C
    end function bml_zero_matrix_C

    function bml_noinit_matrix_C(matrix_type, matrix_precision, n, m, &
        & distrib_mode) bind(C, name="bml_noinit_matrix")
      import :: C_INT, C_PTR
      integer(C_INT), value, intent(in) :: matrix_type
      integer(C_INT), value, intent(in) :: matrix_precision
      integer(C_INT), value, intent(in) :: n
      integer(C_INT), value, intent(in) :: m
      integer(C_INT), value, intent(in) :: distrib_mode
      type(C_PTR) :: bml_noinit_matrix_C
    end function bml_noinit_matrix_C

  end interface

end module bml_c_interface_m
