
# global reference to numpy, scipy, main; initiliazed on load
scipy <- np <- mbdmpy <- mbdmb_pycode <- NULL

# load python modules on load
.onLoad <- function(libname, pkgname) {
  # superassignment to update global python module references
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
  np <<- reticulate::import("numpy", delay_load = TRUE)

  # hack defining global string with py compute_A function
  mbdmb_pycode <<- "
  def compute_A(I, omega_list, eta):
    A = I
    for i in range(0, len(omega_list)):
      A = A - omega_list[i]*eta[i]
    return A
  "
}








