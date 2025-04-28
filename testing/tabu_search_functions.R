pb_tabu_search = function(X, iter = 20, tabu_size = 20){
  M = cov(log(X))

  steps = rep(0, iter)
  B0 = pb_basis(X, method = 'constrained')
  BAL = as.integer(sign(B0[,1]))

  evaluate_balance = function(bal){
    L = bal == -1L; nL = sum(L)
    R = bal == +1L; nR = sum(R)

    ( nR/nL * sum(M[L,L]) + nL/nR * sum(M[R,R]) - 2 * sum(M[R,L]) ) / (nL+nR)
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
       balance = sbp_basis(cbind(BEST), silent = TRUE))

}
