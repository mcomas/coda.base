set.coda = function(x){
  class(x) = c('coda', class(x))
  x
}

#' Printing coordinates
#'
#' The function hides the basis attribute. An option is included to
#' show such basis.
#' @param x coordinates
#' @param ... parameters passed to print function
#' @param basis boolean to show or not the basis with the output
#' @export
print.coda = function(x, ..., basis = getOption('coda.base.basis')){
  x.print = x
  print.methods.list = utils::methods('print')
  orig_class = setdiff(class(x.print), 'coda')
  class(x.print) = orig_class
  print.method = stats::na.omit(match(paste0('print.',orig_class), print.methods.list))[1]
  if(!basis) attr(x.print, 'basis') = NULL
  if(is.na(print.method)){
    print.default(x.print, ...)
  }else{
    utils::getAnywhere(print.methods.list[print.method])$objs[[1]](x.print, ...)
    if(basis){
      B = attr(x.print, 'basis')
      dimnames(B) = list(paste0('P', 1:nrow(B)), colnames(x))
      cat(' Basis:\n')
      print(B)
    }
  }
}
