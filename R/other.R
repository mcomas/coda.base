#' Import data from a codapack workspace
#'
#'
#' @param fname cdp file name
#' @export
read_cdp = function(fname){
  jsonlite_available = requireNamespace("jsonlite")
  if(!jsonlite_available){
    stop("To import CoDaPack's workspace, jsonlite package must be installed.")
  }
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
