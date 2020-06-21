
# global reference to scypi; initiliazed on load
scipy <- NULL

# load python modules on load
.onLoad <- function(libname, pkgname) {
  # superassignment to update global python module references
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
}








