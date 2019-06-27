!> Initialize bml.
module bml_init_m

  use bml_c_interface_m
  use bml_types_m
  use, intrinsic :: iso_C_binding

  implicit none
  private

  public :: bml_initF

contains

  !> Initialize from Fortran when using MPI.
  !!
  !! \param fcomm MPI communicator from Fortran
  subroutine bml_initF(fcomm)

    integer(C_INT), intent(in), optional :: fcomm
    integer :: arg_comm

    if(present(fcomm)) then
       arg_comm = fcomm
    else
       arg_comm = 0
    end if
    call bml_initF_C(arg_comm)

  end subroutine bml_initF

end module bml_init_m
