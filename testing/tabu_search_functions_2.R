library(coda.base)
partial_pb_tabu_search = function(X, lI, iter, tabu_size){
  M = cov(log(X))

  steps = rep(0, iter)
  Xb = sapply(lI, function(I) Reduce(`*`, lapply(I, function(i)X[,i])))
  B0 = pb_basis(Xb, method = 'constrained')
  BAL = as.integer(sign(B0[,1]))

  evaluate_balance = function(bal){
    L = bal == -1L; iL = unlist(lI[L]); nL = length(iL)
    R = bal == +1L; iR = unlist(lI[R]); nR = length(iR)

    ( nR/nL * sum(M[iL,iL]) + nL/nR * sum(M[iR,iR]) - 2 * sum(M[iR,iL]) ) / (nL+nR)
  }

  neighbours = function(bal){
    bal_out = lapply(which(bal != 0L), function(i, bal){
      r_bal = bal
      r_bal[i] = 0L
      r_bal
    }, bal)
    bal_out = bal_out[sapply(bal_out, function(b) sum(b == 1L) > 0 & sum(b == -1L) > 0)]

    bal_left = lapply(which(bal == 0L), function(i, bal){
      r_bal = bal
      r_bal[i] = -1L
      r_bal
    }, bal)
    bal_right = lapply(which(bal == 0L), function(i, bal){
      r_bal = bal
      r_bal[i] = +1L
      r_bal
    }, bal)
    c(bal_out, bal_left, bal_right)
  }
  balance_key = function(bal){
    paste(bal, collapse='')
  }
  format_balance = function(bal){
    str_bal = ifelse(bal == 0, "·", ifelse(bal > 0, '+', '-'))
    paste(str_bal, collapse='|')
  }

  TABU_LIST = list()

  BEST = BAL
  BEST_TABU_SIZE =  0
  BEST_EV = evaluate_balance(BAL)
  BEST_ITER = 0
  # cat(sprintf("Current variance: %0.5f\n", BEST_EV))
  # cat(sprintf("Balance: %s\n\n", format_balance(BEST)))
  ### ITERACIONS TABU
  for(i in 1:iter){
    steps[i] = evaluate_balance(BAL)
    if(length(TABU_LIST) == tabu_size)  TABU_LIST = TABU_LIST[-1]
    TABU_LIST[[length(TABU_LIST)+1]] = balance_key(BAL)

    BAL_N = neighbours(BAL)
    BAL_KEYS = sapply(BAL_N, balance_key)
    BAL_N = BAL_N[!BAL_KEYS %in% TABU_LIST]
    if(length(BAL_N) == 0) break  # all possibilities explored
    BAL_N_EV = sapply(BAL_N, evaluate_balance)

    i_next = which.max(BAL_N_EV)
    BAL = BAL_N[[i_next]]
    BAL_EV = BAL_N_EV[i_next]
    # cat(sprintf("Current variance: %0.5f ➔ ", BAL_EV))
    # cat(sprintf("Balance: %s\n\n", format_balance(BAL)))
    if(BAL_EV > BEST_EV){
      BEST_ITER = i
      BEST = BAL
      BEST_TABU_SIZE = length(TABU_LIST)
      BEST_EV = BAL_EV
    }
  }

  list(iter_best = BEST_ITER,
       tabu_size = BEST_TABU_SIZE,
       steps = steps,
       dim = ncol(Xb)-1,
       lI = lI,
       variance = BEST_EV,
       balance = sbp_basis(cbind(BEST), silent = TRUE))

}
pb_tabu_search = function(X, iter = 100){
  max_steps = 0
  best = partial_pb_tabu_search(X, lapply(1:ncol(X), identity), iter = iter, tabu_size = ncol(X))
  if(best$iter_best > max_steps) max_steps = best$iter_best
  PB = list(best)
  SOLS = list()
  for(k in 1:(ncol(X)-2)){
    if(any(best$balance == 0)){
      lI_top = c(list(unlist(best$lI[best$balance != 0])),
                 lapply(which(best$balance == 0), function(i) best$lI[[i]]))

      top = partial_pb_tabu_search(X, lI_top, iter = iter, tabu_size = length(lI_top))
      if(top$iter_best > max_steps) max_steps = top$iter_best
      SOLS[[length(SOLS) + 1]] = top
    }
    if(sum(best$balance < 0) > 1){
      lI_left = best$lI[best$balance < 0]

      left = partial_pb_tabu_search(X, lI_left, iter = iter, tabu_size = length(lI_left))
      if(left$iter_best > max_steps) max_steps = left$iter_best
      SOLS[[length(SOLS) + 1]] = left
    }
    if(sum(best$balance > 0) > 1){
      lI_right = best$lI[best$balance > 0]

      right = partial_pb_tabu_search(X, lI_right, iter = iter, tabu_size = length(lI_right))
      if(right$iter_best > max_steps) max_steps = right$iter_best
      SOLS[[length(SOLS) + 1]] = right
    }
    i.best = which.max(sapply(SOLS, function(sol) sol$variance))
    best = SOLS[[i.best]]
    PB[[length(PB)+1]] = best
    SOLS = SOLS[-i.best]
  }

  SBP = sapply(PB, function(pb){
    b = rep(0, ncol(X))
    b[unlist(pb$lI[pb$balance>0])] = 1
    b[unlist(pb$lI[pb$balance<0])] = -1
    b
  })
  attr(SBP, 'max_steps') = max_steps
  SBP
}
