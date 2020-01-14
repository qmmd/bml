
module prec

  implicit none

  private

  integer, parameter, public :: sp = kind(1.0)
  integer, parameter, public :: dp = kind(1.0d0)
  character(20), parameter, public :: bml_real = bml_element_real
  character(20), parameter, public :: bml_complex = BML_ELEMENT_COMPLEX


end module prec
