!> Submatrix routines used for graph partitioning
module bml_submatrix_m

  use bml_types_m
  use bml_c_interface_m

  implicit none
  private

  public :: bml_matrix2submatrix_index
  public :: bml_matrix2submatrix
  public :: bml_submatrix2matrix
  public :: bml_adjacency
  public :: bml_adjacency_group
  public :: bml_group_matrix

contains

  !> Determine element indeces for submatrix, given a set of nodes.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Hamiltonian matrix
  !! \param b Graph as a matrix
  !! \param nodelist List of nodes to define submatrix
  !! \param nsize Number of nodes
  !! \param core_halo_index Indeces of core+halo
  !! \param vsize Sizes of core_halo_index and cores
  !! \param double_jump_flag Flag 0=no 1=yes
  subroutine bml_matrix2submatrix_index(b, nodelist, nsize, &
       core_halo_index, vsize, double_jump_flag, a)

    type(bml_matrix_t), optional, intent(in) :: a
    type(bml_matrix_t), intent(in) :: b
    integer(C_INT), target, intent(in) :: nodelist(*)
    integer(C_INT), target, intent(inout) :: core_halo_index(*)
    integer(C_INT), target, intent(inout) :: vsize(*)
    integer, intent(in) :: nsize
    logical, intent(in) :: double_jump_flag

    integer(C_INT) :: cflag

    if (double_jump_flag .eqv. .true.) then
       cflag = 1
    else
       cflag = 0;
    endif

    if (present(a)) then
       call bml_matrix2submatrix_index_C(a%ptr, b%ptr, &
            c_loc(nodelist), nsize, c_loc(core_halo_index), &
            c_loc(vsize), cflag)
    else
       call bml_matrix2submatrix_index_graph_C(b%ptr, &
            c_loc(nodelist), nsize, c_loc(core_halo_index), &
            c_loc(vsize), cflag)
    endif

  end subroutine bml_matrix2submatrix_index

  !> Create contracted submatrix from a set of element indeces.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Matrix
  !! \param b Submatrix
  !! \param core_halo_index Indeces of core+halo
  !! \param lsize Number of indeces in core_halo_index
  subroutine bml_matrix2submatrix(a, b, core_halo_index, lsize)

    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: b
    integer(C_INT), target, intent(in) :: core_halo_index(*)
    integer(C_INT), intent(in) :: lsize

    call bml_matrix2submatrix_C(a%ptr, b%ptr, c_loc(core_halo_index), lsize)

  end subroutine bml_matrix2submatrix

  !> Assemble a contracted submatrix into the final matrix.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Submatrix
  !! \param b Matrix
  !! \param core_halo_index Indeces of core+halo
  !! \param lsize Number of indeces in core_halo_index
  !! \param llsize Number of cores
  subroutine bml_submatrix2matrix(a, b, core_halo_index, lsize, &
       llsize, threshold)

    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: b
    integer(C_INT), target, intent(in) :: core_halo_index(*)
    integer(C_INT), intent(in) :: lsize, llsize
    real(C_DOUBLE), optional, intent(in) :: threshold

    real(C_DOUBLE) :: threshold_

    if (present(threshold)) then
       threshold_ = threshold
    else
       threshold_ = 0
    end if

    call bml_submatrix2matrix_C(a%ptr, b%ptr, c_loc(core_halo_index), lsize, &
         llsize, threshold_)

  end subroutine bml_submatrix2matrix

  !> Assemble adjacency vectors for a matrix.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Matrix
  !! \param xadj Start indeces for each row
  !! \param adjncy Indices of non-zero values
  !! \param base_flag Return 0- or 1-based
  subroutine bml_adjacency(a, xadj, adjncy, base_flag)

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), target, intent(in) :: xadj(*)
    integer(C_INT), target, intent(in) :: adjncy(*)
    integer(C_INT), intent(in) :: base_flag

    call bml_adjacency_C(a%ptr, c_loc(xadj), c_loc(adjncy), base_flag)

  end subroutine bml_adjacency

  !> Assemble adjacency for a matrix based on groups.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Matrix
  !! \param hindex Indeces of nodes in matrix
  !! \param nnodes Number of nodes
  !! \param xadj Start indeces for each row
  !! \param adjncy Indices of non-zero values
  !! \param base_flag Return 0- or 1-based
  subroutine bml_adjacency_group(a, hindex, nnodes, xadj, adjncy, base_flag)

    type(bml_matrix_t), intent(in) :: a
    integer(C_INT), target, intent(in) :: hindex(*)
    integer(C_INT), target, intent(in) :: xadj(*)
    integer(C_INT), target, intent(in) :: adjncy(*)
    integer(C_INT), intent(in) :: nnodes, base_flag

    call bml_adjacency_group_C(a%ptr, c_loc(hindex), nnodes, c_loc(xadj), c_loc(adjncy), base_flag)

  end subroutine bml_adjacency_group

  !> Assemble group-based matrix from given matrix.
  !!
  !! \ingroup submatrix_group_F
  !!
  !! \param a Original matrix
  !! \param b Group-based matrix
  !! \param hindex Indeces of group nodes in matrix
  !! \param ngroups Number of groups
  !! \param threshold Threshold to determine positions
  subroutine bml_group_matrix(a, b, hindex, ngroups, threshold)

    type(bml_matrix_t), intent(in) :: a
    type(bml_matrix_t), intent(inout) :: b
    integer(C_INT), target, intent(in) :: hindex(*)
    integer(C_INT), target, intent(in) :: ngroups
    real(C_DOUBLE), target, intent(in) :: threshold

    call bml_deallocate(b)
    b%ptr = bml_group_matrix_C(a%ptr, c_loc(hindex), ngroups, threshold)

  end subroutine bml_group_matrix

end module bml_submatrix_m
