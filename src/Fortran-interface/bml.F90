!> Main matrix library module.
!!
!! Use this modules in order to use the library.
module bml
  use bml_add_m
  use bml_adjungate_triangle_m
  use bml_allocate_m
  use bml_copy_m
  use bml_diagonalize_m
  use bml_error_m
  use bml_export_m
  use bml_getters_m
  use bml_import_m
  use bml_init_m
  use bml_introspection_m
  use bml_inverse_m
  use bml_multiply_m
  use bml_normalize_m
  use bml_norm_m
  use bml_parallel_m
  use bml_scale_m
  use bml_setters_m
  use bml_shutdown_m
  use bml_submatrix_m
  use bml_threshold_m
  use bml_trace_m
  use bml_transpose_m
  use bml_transpose_triangle_m
  use bml_types_m
  use bml_utilities_m
end module bml
