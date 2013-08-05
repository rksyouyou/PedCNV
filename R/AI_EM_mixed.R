
getInitial <- function(residual,phi,S){
     #.Call("getInitial", residual, phi, S, PACKAGE = "CNVLMM")
    temp <- .Call("getInitial", residual, phi, S)
#    cat('initial values: (sig2) ', temp[["sig2"]], ' (sig2g)', temp[["sig2g"]],'\n')
    unlist(temp)
}


doEM_REML <- function(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q, threshold){
    res = .Call("doEM_REML", curTheta, curK, y, X, yy, Xy, XX, phi, phiInv, trphiInv, S, q, threshold)
    return(res)
}

doAI_REML<-function(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q, itrmax,threshold){
    res = .Call("doAI_REML", curTheta, curK, y, X, yy, Xy, XX, phi, phiInv, trphiInv, S, q, itrmax, threshold)
    return(res)
}

doPolyGenic<-function(envirX,snp,pheno,phi,H0, thresEM,itrmax=20,thresAI,phiInv){

  S  <- length(pheno)
  if(H0) X  <- envirX else
  X  <- cbind(envirX,snp)
  y  <- matrix(pheno,ncol=1)
  yy <- t(y)%*%y
  Xy <- t(X)%*%y
  XX <- t(X)%*%X
  trphiInv<- sum(diag(phiInv))
  q  <- ncol(X)
  residual <- y - X%*%matrix(solve(XX)%*%Xy,ncol=1)

  inis     <- getInitial(residual,phi,S)

  if(inis[2]==0) {
      return(NA)
  } else {
  ## MME method to find the initial values
      updates  <- doEM_REML(inis[2],inis[1]/inis[2],y,X,yy,Xy,XX,phi,phiInv, trphiInv, S, q, thresEM)
  }
######

      curTheta <- updates[1]
      curK     <- updates[2]
      curTheta <- inis[2]
      curK     <- inis[1]/inis[2]

      print('AI algorithm for REML')
      results <- doAI_REML(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q,itrmax,thresAI)

      return(results)
}


## Test Part################

## Rao's score test statistic with the most probable copy numbers
RSTe<- function(envirX,snp,pheno,alpha,Invcov){

  S  <- length(pheno)
  e <- pheno-envirX%*%alpha
  nve <- t(snp)%*%Invcov%*%e
  nvn <- t(snp)%*%Invcov%*%snp
  nvx <- t(snp)%*%Invcov%*%envirX
  xvx <- t(envirX)%*%Invcov%*%envirX
  Trs <- nve%*%solve(nvn-nvx%*%solve(xvx)%*%t(nvx))%*%t(nve)
  pv <- 1-pchisq(Trs,1)
  return(list(STEs=Trs,STEp=pv))

}

##  Rao's score test statistic with the probe intensity measurements


RSTim<- function(envirX,signal,pheno,alpha,InvW,Invcov){

    temp <- scale(signal)
    R <- cor(t(temp))
    S  <- length(pheno)
    N1 <- matrix(1,S,1)
    vz <- Invcov%*%envirX

    s1 <- envirX%*%solve(t(envirX)%*%vz)%*%t(vz)
    InvWN1 <- InvW%*%N1
    s2 <- N1%*%solve(t(N1)%*%InvWN1)%*%t(InvWN1)
    t3 <- t(pheno)%*%t(diag(S)-s1)%*%Invcov%*%(diag(S)-s2)
    t1 <- (diag(S)-s2)%*%R
    t2 <- sum(diag(t1%*%t(diag(S)-s2)%*%Invcov%*%(diag(S)-s1)))
    u <- t3%*%signal
    psi <- t(signal)%*%t(diag(S)-s2)%*%(diag(S)-s2)%*%signal
    v <- t2/sum(diag(t1))*psi
    Trs <- u%*%solve(v)%*%t(u)
    df <- qr(v)$rank
    pv<- 1-pchisq(Trs,df)

    return(list(STIMs=Trs,STIMp=pv,df=df))
}



