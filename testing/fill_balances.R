P = matrix(c(+1,+1,-1,-1,-1,-1, 0, 0, 0,
             +1,-1, 0, 0, 0, 0, 0, 0, 0,
              0, 0,+1,-1, 0, 0, 0, 0, 0), ncol=3)
RES = sbp_basis(P)
iRES = sign(RES)


fill_balances = function(iRES){
  # print(iRES)
  if(ncol(iRES)==0){
    # print(iRES)
    return(sign(ilr_basis(nrow(iRES))))
  }
  # Detecting free components
  free_comp = rowSums(iRES!=0) == 0
  if(sum(free_comp) > 0){
    iRES = cbind(iRES, 1-2*free_comp)
  }

  iROOT = which.max(colSums(iRES!=0))
  iBAL = which(iRES[,iROOT] > 0)

  BAL = matrix(0,nrow(iRES),nrow(iRES)-1)
  BAL[,1:ncol(iRES)] = iRES

  ## fill positive
  pBAL = iRES[+iBAL,-iROOT,drop=FALSE]
  pDIM = nrow(pBAL)-1
  pBAL = pBAL[,colSums(pBAL!=0)>0,drop=FALSE]
  pN = ncol(pBAL)
  # print(pN)
  # print(pDIM)
  # print(pBAL)
  # print(BAL)
  if(pN < pDIM){
    # print(pBAL)
    pBAL = Recall(pBAL)
    BAL[+iBAL,(ncol(iRES)+1):(ncol(iRES)+pDIM-pN)] = pBAL[,(pN+1):pDIM]
  }else{
    # print(pBAL)
    # print(BAL[+iBAL,(ncol(iRES)+1):(ncol(iRES)+pDIM-pN)])
    # print(BAL)
    # print(iRES)
    # print(pDIM)
    # print(pN)
    # BAL[+iBAL,(ncol(iRES)+1):(ncol(iRES)+pDIM-pN)] = pBAL
  }

  # dPos = length(iPOS)
  # Pos = matrix(0, nrow = dPos, ncol = dPos-1)
  # Pos[,1:sum(iPos)] = iRES[+iPOS,,drop=FALSE]
  # iPos = rowSums(Pos!=0) > 0

  ## fill negative
  nBAL = iRES[-iBAL,-iROOT,drop=FALSE]
  nDIM = nrow(nBAL)-1
  nBAL = nBAL[,colSums(nBAL!=0)>0,drop=FALSE]
  nN = ncol(nBAL)
  if(ncol(nBAL) < nDIM){
    # print(nBAL)
    nBAL = Recall(nBAL)
    BAL[-iBAL,(ncol(iRES)+1+pDIM-pN):(1+pDIM+nDIM)] = nBAL[,(nN+1):nDIM]
  }else{
    # BAL[-iBAL,(ncol(iRES)+1+pDIM-pN):(1+pDIM+nDIM)] = nBAL
  }

  BAL
  # Neg = iRES[-lPOS,-iROOT,drop=FALSE]
  # iNeg = rowSums(Neg!=0) > 0
  #
  # B = matrix(0, nrow(iRES), sum(free_comp)-1)
  # B[free_comp,] = sign(ilr_basis(sum(free_comp)))
  # iRES = cbind(iRES, 1-2*free_comp,B)
  #
  #
  # patterns = apply(iRES+1,1,paste,collapse='')
  # # patterns_tab = table(patterns)
  # not_included_pattern = paste(rep(1,ncol(iRES)), collapse='')
  # if(not_included_pattern %in% patterns){
  #   iRES = cbind(iRES, 2*(patterns != not_included_pattern)-1 )
  # }
  #
  # iROOT = which.max(colSums(iRES!=0))
  # bal = iRES[,iROOT]
  # bal_p = iRES[bal>0,-iROOT,drop=FALSE]
  # i_p = which(colSums(bal_p!=0)!=0)
  # bal_p = bal_p[i_p,,drop=FALSE]
  # if(nrow(bal_p) > 2){
  #   new_bals = matrix(0,nrow=nrow(iRES),ncol=nrow(iRES)-1)
  #   bal_p_filled = fill(bal_p[,-iROOT,drop=FALSE])
  # }
  # bal_n = iRES[bal<0,,drop=FALSE]
  # bal_n = bal_n[,colSums(bal_n!=0)!=0,drop=FALSE]
  # if(nrow(bal_n) > 2){
  #   bal_n_filled = fill(bal_n[,-iROOT,drop=FALSE])
  # }
  # iRES
}
fill_balances(iRES)

## EROR
X = data.frame(a=1:2, b=2:3, c=4:5, d=5:6, e=10:11, f=100:101, g=1:2)
iRES = sbp_basis(b1 = a~b,
                 b2 = b1~c,
                 b3 = b2~d, data = X) |> sign()
fill_balances(iRES)

