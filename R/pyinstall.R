##############################################################
# Installation function for python dependencies
##############################################################

install_pydep <- function(method = "auto", conda = "auto") {
  reticulate::py_install("scipy", method = method, conda = conda)
}
