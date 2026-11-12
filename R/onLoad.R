.onLoad <- function(libname, pkgname){

  #Set the options
  .setOptions()

  # Call to load c++ libraries
  olr_c(1:10)
}
.setOptions <- function(){
  if (is.null(getOption("coda.base.basis"))) {
    options("coda.base.basis" = FALSE)
  }
}
