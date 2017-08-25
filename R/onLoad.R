.onLoad <- function(libname, pkgname){

  #Set the options
  .setOptions()
}
.setOptions <- function(){
  options("coda.base.basis" = FALSE)
}
