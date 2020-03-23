FORTRAN TESTS
=============

The tests are driven by a general executable created when the code is
compiled with `BML_TESTING=yes`. This driver is called bml-testf
compiled with the `testf.F90` source.

Every low level source code of the type name_typed.F90 is
pre-processed using the `/scripts/convert-template.in` to change to
the particular element kind and presicion. Two dummy varibles are
used:

  - `DUMMY_KIND`: That gets replaced with either `real` or `complex`
  - `DUMMY_PREC` or `_MP`: That gets replaced with `SP/_SP` of `DP/_DP` (defined in prec.F90)

There is are `exmple_template*` files that can be used as starting
point to add a particular test.

# Conventions and rules

The general driver take four variables (this can be extended as
needed). These variables are:

  - `test_name`: The name of the test
  - `matrix_type`: The matrix format (matrix format and matrix type
                   are the same thing)
  - `element_type`: The element "kind" and "precision". For example
                    double_real, which gets converted to real(8) at
                    the lowest level.

NOTE: Try to be as explicit as possible in naming the variables.
