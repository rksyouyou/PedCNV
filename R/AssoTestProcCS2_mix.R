

Inde <- function(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S){
    print('The individuals are independent, LM is used!')
    ## changed
    if(!H0){
        res <- lm(pheno~envirX+clusRes-1)
        len <- length(res$coefficients)
        alpha <- res$coefficients[1:(len-1)]
        bet <- res$coefficients[len]
    }else{
        res <- lm(pheno~envirX-1)
        alpha <- res$coefficients
        bet <- 0
    }
    zalpha <- envirX%*%alpha
    nbet <- clusRes*bet
    e <- pheno-zalpha-nbet
    sig2 <- as.numeric(var(e))
    return(list(alpha=alpha,bet=bet,sig2=sig2,e=e))
}


CNVtypeAnay <- function(pheno,pX,envirX,phi,S,FM,N,threshold, bet, alpha, sig2, sig2g,H0,thresAI,thresEM,itermax){

  seed<- sample(c(1:1000000),1)
  ## print(paste0("seed is ", seed))
  set.seed(seed)
  p <- rep(1/N,N)
  Ind.num <- dim(pX)[1]
  clusRes <- rep(NA,Ind.num)
  del <- list()
  Invdel <- list()
  mu <- list()
  LDdel <- rep(NA,N)
  Ddel <- rep(NA,N)
  mu1 <- rep(NA,N)
  q <- dim(pX)[2]


  ires <- ClusProc(signal=pX,N=N,varSelection='RAW')
  part1 <- data.frame(X=pX)
  part2 <- data.frame(Kn=ires$silWidth$silRes)
  temp <- part2[match(rownames(part1),rownames(part2)),]
  KX <- cbind(X=part1,Kn=temp[,1])


  clustemp <- KX$Kn+N
  t <- rep(NA,3)
  for(i in 1:N)
    mu1[i] <- mean(apply(as.matrix(KX[KX$Kn==i,c(1:q)]),2,mean))
  temp <- cbind(mu1,(N+1):(2*N),as.numeric(factor(mu1)))
  for(i in 1:N) t[i] <- temp[which(temp[,3]==i),2]
  for(i in 1:S)
    for(j in 1:N)
      if(clustemp[i]==t[j]) KX$Kn[i] <- j
  ## initial value for the means and covariance
  for(i in 1:N){
    mu[[i]] <- apply(as.matrix(KX[KX$Kn==i,c(1:q)]),2,mean)
    del[[i]] <- cov(as.matrix(KX[KX$Kn==i,c(1:q)]))
    Ddel[i] <- det(del[[i]])
    LDdel[i] <- log(Ddel[i])
    Invdel[[i]] <- solve(del[[i]])
  }
  
  ## variables preparing
  iter <- 0
  logLold <- Inf
  pin <- matrix(NA,S,N)
  P <- table(KX$Kn)/S
  clusFam <- rep(NA,FM)
  phiInv <- solve(phi)
  
  while(1){

    iter=iter+1
    for(i in 1:S)
      for(n in 1:N)
        pin[i,n]  <- Ddel[n]^(-1/2)*exp(-1/2*(pX[i,]-mu[[n]])%*%Invdel[[n]]%*%as.matrix(pX[i,]-mu[[n]]))*P[n]
    ## set the method to be used
    for(i in 1:S)
      clusRes[i] <- which(pin[i,]==max(pin[i,]))
    for(i in 1:N)
      P[i] <- mean(clusRes==i)
    Sdata <- data.frame(pX=pX,clusRes=clusRes)
    ## updata the signal model
    for(i in 1:N){
      mu[[i]] <- apply(as.matrix(Sdata[Sdata$clusRes==i,c(1:q)]),2,mean)
      del[[i]] <- cov(as.matrix(Sdata[Sdata$clusRes==i,c(1:q)]))
      Ddel[i] <- det(del[[i]])
      LDdel[i] <- log(Ddel[i])
      Invdel[[i]] <- solve(del[[i]])
    }
    ## updata the phenotype model
    pheno[which(is.na(pheno))] <- 0

    if(all(phi[lower.tri(phi)]==0,phi[upper.tri(phi)]==0)){
        temp <- Inde(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S)
        sig2g <- 0
        bet <- temp$bet
        sig2 <- temp$sig2
        alpha <- temp$alpha
        e <- temp$e
    }else{
        res <- doPolyGenic(envirX=envirX,snp=clusRes,pheno=pheno,phi=phi,phiInv=phiInv,thresEM=thresEM,thresAI=thresAI,H0=H0)

        if(length(res)==1) {
            temp <- Inde(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S)
            sig2g <- 0
            bet <- temp$bet
            sig2 <- temp$sig2
            alpha <- temp$alpha
            e <- temp$e
        }else {
            print('The individuals are correlated, LMM is used!')
            sig2 <- res$sig2
            sig2g <- res$sig2g
            len <- length(res$para)
            if(!H0){
                alpha <- res$para[1:(len-1)]
                bet <- res$para[len]}else{
                    alpha <- res$para[1:len]
                    bet <- 0
                }
            zalpha <- envirX%*%alpha
            nbet <- clusRes*bet
            e <- pheno-zalpha-nbet
        }        
    }

    V <- sig2*diag(S)+sig2g*phi
    ## loglikelihood calculation
    logL_sig <- 1
    for(i in 1:S)
        logL_sig <- logL_sig-1/2*log(Ddel[clusRes[i]])-1/2*(pX[i,]-mu[[clusRes[i]]])%*%solve(del[[clusRes[i]]])%*%as.matrix(pX[i,]-mu[[clusRes[i]]])
    logL <-logL_sig-1/2*determinant(V)$modulus-1/2*t(e)%*%solve(V)%*%e+sum(log(P[clusRes]))
    ##if(abs(logL-logLold)<threshold) ttt <- 1 else ttt <- 0
    #print(ttt)
    if(abs(logL-logLold)<threshold) break
    if(iter>itermax) break
    logLold <- logL
    print(paste("The logliklihood without of costants is now",logL,"."))
}

  return(list(clusRes=clusRes,bet=bet,alpha=alpha,sig2=sig2,sig2g=sig2g,logL1=round(logL,6)))
  
}




##' This function tests the association of CNV with continuous trait based on CNV-linear mixed model. Two statistics are provided for different strategies with the intensity measurement.
##'
##' @title CNV association testing
##' @param signal The matrix of intensity measurements. The rownames must be cosistent with the Individual ID in PED file.
##' @param ped The PED file which follows the format defined in PLINK. 
##' @param envirX The matrix of enviromental variables. The intercept should be included in it if it's necessary in further analysis.
##' @param phi The matrix of correlation between individuals. 
##' @param N Number of clusters one wants to fit to the data. N needs to be larger than 1 and if it is 1, error will be returned.
##' @param varSelection Factor. For specifying how to handle the intensity values. It must take value on 'RAW', 'PC1', 'PC.9' and 'MEAN'. If the value is 'RAW', then the raw intensity value will be used. If it is 'PC.9' (the default), then the first several PCA scores which account for 90% of all the variance will be used. If the value is 'PC1', then the first PCA scores will be used.If the value is 'MEAN', the mean of all the probes will be used. 
##' @param H0 Logicals. If the value is TRUE (the default), then the parameters  will be estimated under null hypotheis only. If the value is FALSE, then the parameters will be estimated under both null hypothesis and alternative hypothesis. 
##' @param threshold Optinal number of convergence threshold. The iteration stops if the absolute difference of loglikelihood between successive iterations is less than it. The default threshold 1e-05 will be used if it's missing.
##' @param itermax Optinal. The iteration stops if the time of iteration is large than this value. The default number 8 will be used if it's missing.
##' @param thresEM Optional number of convergence threshold in the EM (expectation-maximization method) procedure. The default threshold 0.005 will be used if it's missing.
##' @param thresAI Optional number of convergence threshold in the AI (average information method) procedure. The default threshold 1e-05 will be used if it's missing.
##' @return It returns object of class 'asso'. 'asso' is a list containing at least model estimation results under H0, the results under H1 will be included if the value of {H0} is FALSE.
##' \item{asso.test}{The association test statistics and p-value from score test. One is based on the estimate of the most probable copy numbers, denoted by T1. Another one is based on the probe intensity measurement, denoted by T2.}
##' \item{para.H0}{The parameter estimations for the best fit under H0.}
##' \item{clus.H0}{The clustering assignment for each individual under H0.}
##' \item{para.H1}{The parameter estimations for the best fit under H1.}
##' \item{clus.H1}{The clustering assignment for each individual under H1.}
##' @author Meiling Liu \email{meiling.sta@@gmail.com} and Sungho Won \email{sunghow@@gmail.com}
##' @examples
##' # Load data and correlation matrix
##' data(dat)
##' data(phi)
##' signal <- dat$Inx
##' ped <- dat$ped
##' envirX <- dat$x
##' # Fit the data under the assuption that there are 3 clusters
##' # fit.pc <- AssoTestProc(signal=signal,ped=ped,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' @export

AssoTestProc <- function(signal,ped,envirX,phi,N,varSelection=c('RAW','PC.9','PC1','MEAN'),H0=TRUE,threshold=1e-05,itermax=8,thresEM=0.005,thresAI=1e-05){  


    rn_signal <- row.names(signal)
    rn_envirX <- row.names(envirX)
    iid <- ped[,2]
    if(sum(!(rn_signal%in%iid))>0) {
        cat('The row name of signal data is not coincident with individual ID!\n')
        return(0)
    }
    if(sum(!(rn_envirX%in%iid))>0) {
        cat('The row name of covariant data is not coincident with individual ID!\n')
        return(0)
    }
    len <- sum(!duplicated(ped[,1]))
    pheno <- as.numeric(ped[,6])
    S <- length(pheno)
    p <- rep(1/N,N)
    INsnp <- rbinom(S,2,p)
    X  <- cbind(envirX,INsnp)
    y  <- matrix(pheno,ncol=1)
    yy <- t(y)%*%y
    Xy <- t(X)%*%y
    XX <- t(X)%*%X
    phiInv  <- solve(phi)
    trphiInv<- sum(diag(phiInv))
    q  <- ncol(X)
    residual <- y - X%*%matrix(solve(XX)%*%Xy,ncol=1)
    inis     <- getInitial(residual,phi,S)
    in_sig2g <- inis[2]
    in_sig2 <- inis[1]
    pca <- princomp(signal)
    pX1 <- signal%*%loadings(pca)[,1]
    cut <- 1
    if(varSelection=='PC.9'){
        prop <- 0.9
        pca <- princomp(signal) ## PCA
        sds <- pca$sdev
        vars <- sds^2
        varprop  <- vars/sum(vars)
        cumvars <- as.vector(cumsum(varprop))
        while(cumvars[cut]<prop)  cut=cut+1
        print(paste0("we use Comp.1 to Comp.",cut))
        Invcov <- matrix(0,nrow=cut,ncol=cut)
        diag(Invcov) <- 1/vars[1:cut]
        comptable <- data.frame(sdev=pca$sdev,vars=vars,cumu=cumvars)
        print("Variable Selection")
        print(t(comptable[1:(cut+1),]))
        coef <- pca$loadings[,1:cut]
        pX <- signal%*%coef
  }
  if(varSelection=='RAW') {
    pX <- signal
    cut <- ncol(pX)
  }
  if(varSelection=='MEAN')    pX <- as.matrix(apply(signal,1,mean))
  if(varSelection=='PC1'){
    pca <- princomp(signal)
    pX <- signal%*%loadings(pca)[,1]
  }
    ## under all
    if(!H0){
        fit <- lm(pheno~pX1+envirX-1)
        bet <- fit$coefficient[1]
        alpha <- fit$coefficient[-1]
        res_H1<- CNVtypeAnay(pX=pX,pheno=pheno,envirX=envirX,phi=phi,S=S,FM=len,N=N,threshold=threshold,bet=bet,alpha=alpha,sig2=in_sig2,sig2g=in_sig2g,H0=FALSE,thresEM=thresEM,thresAI=thresAI,itermax=itermax)
        logL <- res_H1$logL
    }
    
    fit <- lm(pheno~envirX-1)
    bet <- 0
    alpha <- fit$coefficient
    res_H0 <- CNVtypeAnay(pX=pX,pheno=pheno,envirX=envirX,phi=phi,S=S,FM=len,N=N,threshold=threshold,bet=bet,alpha=alpha,sig2=in_sig2,sig2g=in_sig2g,H0=TRUE,thresEM=thresEM,thresAI=thresAI,itermax=itermax)
    
    alpha_H0<- res_H0$alpha
    sig2_H0 <- res_H0$sig2
    sig2g_H0 <- res_H0$sig2g
    logL_H0 <- res_H0$logL
    clusRes <- res_H0$clusRes
  ##  Rao's score test statistic with the most probable copy numbers
    Invcov <- solve(sig2g_H0*phi+sig2_H0*diag(S))
    res <- RSTe(envirX=envirX,snp=clusRes,pheno=pheno,alpha=alpha_H0,Invcov=Invcov)
    STEs <- res$STEs
    STEp <- res$STEp
    res <- RSTim(envirX=envirX,signal=signal,pheno=pheno,alpha=alpha_H0,Invcov=Invcov,InvW=phiInv)
    STIMs <- res$STIMs
    STIMp <- res$STIMp
    if(H0){
        resfinal <- list(asso.test=list(T1s=STEs,T1p=STEp,T2s=STIMs,T2p=STIMp),para.H0=list(bet=res_H0$bet,alpha=res_H0$alpha,sig2=res_H0$sig2,sig2g=res_H0$sig2g),clus.H0=list(clusRes=res_H0$clusRes))
    } else
        {
      resfinal <- list(asso.test=list(T1s=STEs,T1p=STEp,T2s=STIMs,T2p=STIMp),para.H0=list(bet=res_H0$bet,alpha=res_H0$alpha,sig2=res_H0$sig2,sig2g=res_H0$sig2g),clus.H0=list(clusRes=res_H0$clusRes),para.H1=list(bet=res_H1$bet,alpha=res_H1$alpha,sig2=res_H1$sig2,sig2g=res_H1$sig2g),clus.H1=list(clusRes=res_H1$clusRes))
  }
 
    class(resfinal) <- 'asso'
    return(resfinal)
}







##' Prints formatted results from the association study returned by AssoTestProc.
##'
##' @title Prints association study results
##' @param x The association study results obtained from the AssoTestProc.
##' @param ... Usual arguments passed to the print function.
##' @author Meiling Liu \email{meiling.sta@@gmail.com} and Sungho Won \email{sunghow@@gmail.com}
##' @examples
##' # Load data and correlation matrix
##' data(dat)
##' data(phi)
##' signal <- dat$Inx
##' ped <- dat$ped
##' envirX <- dat$x
##' # Fit the data under the assuption that there are 3 clusters
##' # fit.pc <- AssoTestProc(signal=signal,ped=ped,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' # print(fit.pc)
##' @export
print.asso <- function(x, ...){

    cat('Pvalue of Association Test\n')
    cat('Using the most probable copy number:',x$asso.test$T1p,'.\n')
    cat('Using the probe intensity measures:',x$asso.test$T2p,'.\n')

    cat('\n\n Under H0:\n')
    cat('The coefficent:\n')
    temp <- c(bet=0,alpha=round(x$para.H0$alpha,4))
    print(temp,quote=FALSE,row.names=FALSE)
    cat('The estimated variance:\n')
    temp <- data.frame(sig2=x$para.H0$sig2,sig2g=x$para.H0$sig2g)
    print(temp,quote=FALSE,row.names=FALSE)
    ## cat('The loglikelihood function value:',round(x$H0$logL,4),'.\n')
    if(length(x$para.H1)>1){
        cat('\n\n Under H1:\n')
        cat('The coefficent:\n')
        temp <- c(bet=round(x$para.H1$bet,4),alpha=round(x$para.H1$alpha,4))
        print(temp,quote=FALSE,row.names=FALSE)
        cat('The estimated variance:\n')
        temp <- data.frame(sig2=x$para.H1$sig2,sig2g=x$para.H1$sig2g)
        print(temp,quote=FALSE,row.names=FALSE)
        ## cat('The loglikelihood function value:',round(x$H1$logL,4),'.\n')
    }
}

