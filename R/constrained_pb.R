fBalChip<-function(Xcoda){
  # given  a coda set Xcoda
  # this function calls fBPMaxOrthNewChip for a searching using the algortihm for Constrained PCs
  # and it returns ALL balances with D parts
  # that maximizes the variance
  # The balances are sorted by the percentatge of variance
  #
  # Returns a list: balances and variance
  # Bres= balances
  # Vres= variance of balances

  numbal=ncol(Xcoda)-1
  # call the recursive function
  res<-fBPMaxOrthNewChip(Xcoda)
  Bres<-res$bal
  balname<-paste("bal",1:nrow(Bres),sep="")
  rownames(Bres)<-balname
  colnames(Bres)<-colnames(Xcoda)
  Vres<-res$varbal
  #
  # sort by expl var
  vopt<-res$varbal
  # sort variance
  vsopt<-sort(vopt,decreasing = TRUE,index.return=TRUE)
  #
  # assign variance explained already ordered
  Vres<-vsopt$x
  #
  # assign balances same order
  Bres<-Bres[vsopt$ix,]
  #
  # return results: balances and variances
  return(list(bal=Bres,varbal=Vres))
}


fBPMaxOrthNewChip<-function(Y,angle=TRUE)
{

  # recursion: given a coda set Y
  # return list of principal balances basis
  # that maximizises the variance
  # searching by the NO-FULL and COMPLETING the
  # SBP using loop (UP: fBUpChi.r) and recursive (DOWN) schemes
  # both based on Chipman procedure

  numpart=ncol(Y)
  numbal=ncol(Y)-1
  B=c()
  #B=matrix(0,numbal,numpart) to save balances
  V=c()
  #V=matrix(0,1,numbal) to save variances

  #first optimal in data set Y
  res<-fBalChipman(Y,angle=angle)
  B<-res$bal
  V<-res$varbal
  # if necessary GO UP to complete
  if (sum(B==0)>0){
    res<-fBPUpChi(Y,B)
    B=rbind(B,res$bal)
    V=cbind(V,res$varbal)
  }
  # control number of balances added
  if (is.vector(B)) B<-matrix(B,1,length(B))
  numbaladd<-nrow(B)-1

  ### GO DOWN THE CURRENt LIST AND THE FIRST
  ## first go down from the first optimal balance

  usenum<-(B[1,]>0)
  useden<-(B[1,]<0)
  # GO DOWN from numerator of the first optimal balance
  if(sum(usenum)>1){
    resP<-fBPMaxOrthNewChip(Y[,usenum],angle=angle)
    Bx<-matrix(0,length(resP$varbal),numpart)
    Bx[,usenum]<-resP$bal
    B<-rbind(B,Bx)
    V<-cbind(V,resP$varbal)
  }# end if
  # GO DOWN from denominator of the first optimal balance
  if(sum(useden)>1){
    resP<-fBPMaxOrthNewChip(Y[,useden],angle=angle)
    Bx<-matrix(0,length(resP$varbal),numpart)
    Bx[,useden]<-resP$bal
    B<-rbind(B,Bx)
    V<-cbind(V,resP$varbal)
  }# end if

  # REVISIT list of balances added GO UP so as to complete the SBP if necessary GO DOWN by the POSITIVE

  if (numbaladd > 0){
    for (k in 2:(1+numbaladd)){
      usepos=(B[k,]>0)
      if (sum(usepos)>1) {
        resP<-fBPMaxOrthNewChip(Y[,usepos],angle=angle)
        Bx<-matrix(0,length(resP$varbal),numpart)
        Bx[,usepos]<-resP$bal
        B<-rbind(B,Bx)
        V<-cbind(V,resP$varbal)
      }#end if2
    }# end for
  }# end if1

  # return results
  #
  V<-as.matrix(V,1,lenght(V))
  #
  return(list(bal=B,varbal=V))

}


fBPUpChi<-function(Yp,b)
{

  # given coda set Yp
  # and given a balance with some zero
  # return list of PARENT principal balances basis
  # that maximizises the variance
  # searching by the NO-FULL {0, -1,+1} and COMPLETING the
  # SBP using a loop (UP) scheme by the CHIPMAN procedure

  npart=ncol(Yp)
  nbal=ncol(Yp)-1
  # to save baances and variances
  Bal=c()
  VarB=c()

  usezero<-sum(b==0)

  #log-transfo data
  lYp=as.matrix(log(Yp))
  # while it is not the full balance go up
  k=0
  while (usezero>0){
    # new balance
    k<-k+1
    # non-zero in the will be in the denominator
    den<-sum(b!=0)
    # for only one zero we get the full
    if (usezero==1){

      b[b!=0]<--sqrt(1/((den+1)*den))
      b[b==0]<-sqrt(den/(den+1))

      VarB<-cbind(VarB,var(as.vector(lYp%*%b)))
      Bal<-rbind(Bal,b)
      usezero<-0
    }
    # for more than one zero we explore other {0,+1} combinations
    else{
      # create the combination by CHIPMAN procedure
      # search the maximum balance
      clrC<-log(Yp[,b==0]) - rowMeans(log(Yp[,b==0]))
      # PCs
      #
      pcClr<-prcomp(clrC)
      #first PC: PC1
      bx<-pcClr$rotation[,1]

      #
      # look for change of sign
      if (abs(min(bx))>max(bx)){bx<--bx}
      # force zeros to the other sign
      bx[bx<0]<-0
      # matrix of {0,+1} possibilities
      M<-matrix(0,sum(bx>0),length(bx))
      # sort
      bxsort<-sort(bx,decreasing = TRUE,index.return=TRUE)
      # index
      col<-bxsort$ix
      # create M
      for (i in 1:nrow(M)){
        M[i,col[1:i]]<-abs(bx[col[1:i]])
      }
      # sign
      M<-sign(M)

      # create a balance
      balax<-b
      # old non-zero to denominator
      balax[b!=0]<--1
      # search the max variance
      VarSBPx<-matrix(0,1,nrow(M))
      balsx<-c()
      for (i in 1:nrow(M))
      {
        # old non-zero to denominator
        balax[b!=0]<--1
        # take one possibility
        balax[b==0]<-M[i,]
        # create the coefficients
        num<-sum(balax==1)
        balax[balax==1]<-sqrt(den/((den+num)*num))
        balax[balax==-1]<--sqrt(num/((den+num)*den))
        balsx=rbind(balsx,balax)

        VarSBPx[i]<-var(as.vector(lYp%*%balax))
      }
      VarB=cbind(VarB,max(VarSBPx))
      Bal=rbind(Bal,balsx[VarSBPx==max(VarSBPx),])

      rm(M)
      usezero<-sum(Bal[k,]==0)
      b<-Bal[k,]
      # end else GO UP
    }

    # end GO UP while
  }
  # return results
  return(list(bal=Bal,varbal=VarB))
  # end function
}

fBalChipman<-function(C,angle=TRUE){
  # given a coda set C
  # This function returns the balance with D parts
  # that is close (angle) to the first PC
  # If "angle=FALSE" the Max the var of scores
  #
  # columns
  col<-dim(C)[2]
  nbal<-col-1
  # clr-transfo
  clrC<-log(C) - rowMeans(log(C))

  # PCs
  #
  pcClr<-prcomp(clrC)
  #first PC: PC1
  pcClr1<-pcClr$rotation[,1]

  balsig<-sign(pcClr1)
  bal<-matrix(0,nbal,col)
  colnames(bal)<-colnames(C)
  # balances associated to the PCs
  # first bal
  bal[1,pcClr1==max(pcClr1)]<-1
  bal[1,pcClr1==min(pcClr1)]<--1
  numbal=1

  # other bal
  if (col>2){
    numbal=numbal+1
    while (numbal<col){
      bal[numbal,]<-bal[numbal-1,]
      useonly<-(bal[numbal-1,]==0)
      bal[numbal,abs(pcClr1)==max(abs(pcClr1[useonly]))]<-balsig[abs(pcClr1)==max(abs(pcClr1[useonly]))]
      numbal=numbal+1
    }#end while
  }#end if

  # coefficients & angle
  VarSBP<-rep(0,nbal)
  for (f in 1:nbal) {
    den<-sum(bal[f,]==-1)
    num<-sum(bal[f,]==1)
    bal[f,bal[f,]==1]<-sqrt(den/((den+num)*num))
    bal[f,bal[f,]==-1]<--sqrt(num/((den+num)*den))
    # variance of the balance
    VarSBP[f]<-abs(sum(bal[f,]*pcClr1))
  }
  # log-trasnform
  lC<-as.matrix(log(C))
  mvar=var(as.vector(lC%*%bal[VarSBP==max(VarSBP),]))

  if (!angle) {

    # calculate variance in the balance direction
    VarSBP<-rep(0,nbal)

    for (i in 1:nbal)
    {
      Proj<-as.vector(lC%*%(bal[i,]))
      VarSBP[i]<-var(Proj)
    }# end for
    mvar=max(VarSBP)
  }# end if

  # return results
  return(list(bal=bal[VarSBP==max(VarSBP),],varbal=mvar))


}
