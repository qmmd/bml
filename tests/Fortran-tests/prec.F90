
module prec

  use bml

  implicit none

  private

  integer, parameter, public :: sp = kind(1.0)
  integer, parameter, public :: dp = kind(1.0d0)
  character(20), parameter, public :: bml_real = bml_element_real
  character(20), parameter, public :: bml_complex = bml_element_complex


end module prec
