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
      cat(' Basis:\n')
      print(B)
    }
  }
}

#' Import data from a codapack workspace
#'
#'
#' @param fname cdp file name
#' @export
read_cdp = function(fname){
  file = jsonlite::read_json(fname)
  ldat = lapply(file$dataframes, function(df){
    vars = df$variables
    vnames = sapply(1:length(vars), function(i) vars[[i]]$n)
    vtype = sapply(1:length(vars), function(i) vars[[i]]$t)
    ldat_ = lapply(1:length(vars), function(i){
      if(vtype[i] == 2){
        values = as.numeric(unlist(vars[[i]]$a))
        values[is.nan(values)] = NA
        sel_zero = sapply(vars[[i]]$a, names) == 'l'
        if(sum(sel_zero) > 0){
          dl = rep(0, length(values))
          dl[sel_zero] = values[sel_zero]
          values[sel_zero] = 0
          return(cbind(values, dl))
        }else{
          return(cbind(values))
        }
      }else{
        values = unlist(vars[[i]]$a)
        return(cbind(values))
      }
    })
    dat = mapply(function(dat, nm){
      if(ncol(dat) == 1){
        colnames(dat) = nm
      }else{
        colnames(dat) = c(nm, paste0(nm, ".dl"))
      }
      dat
    }, ldat_, vnames, SIMPLIFY = FALSE)
    as.data.frame(dat)
  })
  if(length(ldat) == 1){
    return(ldat[[1]])
  }else{
    names(ldat) = sapply(file$dataframes, function(df) df$name)
    return(ldat)
  }
}
