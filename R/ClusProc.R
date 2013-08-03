ClusProc0 <- function(signal,Num,prop,threshold,thres_sil,thres_MAF,cut){

  seed<- sample(c(1:1000000),1)
  print(paste0("seed is ", seed))
  set.seed(seed)
  
  pX <- as.matrix(signal)
  S <- dim(pX)[1]
  clusRes <- rep(NA,S)
  ## the initial value 
  del <- list()
  mu <- list()

  kcl <- kmeans(pX,Num,nstart=10)
  KX <- data.frame(X=pX,Kn=kcl$cluster)

  
  for(i in 1:Num){
    mu[[i]] <- apply(as.matrix(KX[KX$Kn==i,c(1:cut)]),2,mean)
    del[[i]] <- cov(as.matrix(KX[KX$Kn==i,c(1:cut)]))
  }

  alpha=rep(1/Num,Num)
  iter <- 0
  logLold <- Inf
  pin <- matrix(NA,S,Num)
  M=rep(NA,Num)
  LDdel=rep(NA,Num)

  while(1){
    ## estimate n and alpha
    iter=iter+1
    for(i in 1:S)
      for(n in 1:Num)
        pin[i,n]  <- 1/(2*pi*det(del[[n]]))*exp(-1/2*(pX[i,]-mu[[n]])%*%solve(del[[n]])%*%as.matrix(pX[i,]-mu[[n]]))*alpha[n]
    ## set the method to be used
    for(i in 1:S)
      clusRes[i] <- which(pin[i,]==max(pin[i,]))
    for(i in 1:Num)
      alpha[i] <- mean(clusRes==i)
    Sdata <- data.frame(pX=pX,clusRes=clusRes)
    ## estimate the signal model
    for(i in 1:Num){
      mu[[i]] <- apply(as.matrix(Sdata[Sdata$clusRes==i,c(1:cut)]),2,mean)
      del[[i]] <- cov(as.matrix(Sdata[Sdata$clusRes==i,c(1:cut)]))
      LDdel[i] <- log(det(del[[i]]))
      M[i] <- sum(Sdata$clusRes==i)
    }          
    sum <- 0
    for(i in 1:S)
      sum <- sum+(pX[i,]-mu[[clusRes[i]]])%*%solve(del[[clusRes[i]]])%*%as.matrix(pX[i,]-mu[[clusRes[i]]])
    logL <- -1/2*(sum(M*LDdel)+sum)
    if(abs(logL-logLold)<threshold) break
    logLold <- logL
  }

  logL <- logL-1/2*S*cut*log(2*pi)
  print(paste("The logliklihood is",logL,"when clustering number is ",Num))

  sil <- silWidth(Sdata,thres_sil=thres_sil,thres_MAF=thres_MAF)
  return(list(logL=logL,sil=sil)) 
}

##' This function chooses the optimal number of clusters. The number of cluster which maximizes the silhouette width is selected.
##' no further details yet
##' @title This function returns the optimal number of clusters and some descriptive statistics for each object.
##' @param signal for specifying the probe intensity measures. It can be a numeric vector, matrix, data.frame or list.
##' @param N for applying the given number of clusters to the data. N needs to be larger than 1 and if it is 1, error will be returned.
##' @param prop a numeric value between 0 and 1. This is used to select the number of PC scores, and PC scores corresponding to the minimum number of eigenvalues whose cumulative proportion is larger than \textit{prop} will be selected. The default is 0.9 and it is meaningful if the value of {varSelection} is TRUE. Otherwise it is ignored. See {varSelection}.
##' @param threshold Convergence threshold. The iteration stops if the absolute difference of loglikelihood between successive iterations is less than {threshold}.
##' @param varSelection a factor. For specifying how to handle the intensity values. It must take value on 'RAW', 'PCA', 'PCA1' and 'MEAN'. If the value is 'RAW', then the raw intensity value will be used. If it is 'PCA', then the first several PCA scores which account for certain proportion (this value will be defined in parameter {prop}) of all the variance will be used. If the value is 'PCA1', then the first PCA scores will be used.If the value is 'MEAN', the mean of all the probes will be used. 
##' @return It returns object of class 'clust'. 'clust' is a list containing at least following componets:
##' @param clusNum The best clustering number given the parameter {N}.
##' @param logL The loglikelihood value of the model when the clustering number is {clusNum}.
##' @param silRes The mean of silhouette width for all the object when the clustering number is {clusNum}.
##' @param clusRes The clustering detail for each object when the cluster number is {clusNum}.
##' @param sil The sihouette width for each object when the cluster number is {clusNum}.
##' @author meiling
##' @export

ClusProc <- function(signal,N,prop=0.9,threshold=1e-05,varSelection=c('PCA','RAW','MEAN','PCA1'),thres_sil,adjust=TRUE,thres_MAF=0.01,scale=FALSE){

  sX0 <- as.matrix(signal)
  if(scale) sX <- scale(sX0) else sX <- sX0
  cut <- 1
  if(varSelection=='PCA'){
    ## sX0 <- scale(sX)
    pca <- princomp(sX) ## PCA
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
    pX <- sX%*%coef
  }
  if(varSelection=='RAW') {
    pX <- sX
    cut <- ncol(pX)
  }
  if(varSelection=='MEAN')    pX <- apply(sX,1,mean)
  if(varSelection=='PCA1'){
    pX <- prcomp(sX)$x[,1]
  }
  
    if(1%in%N) stop('The assigned clustering number must be larger than 1!')
    res <- list()
    for(i in 1:length(N)){
        Num <- N[i]
        res[[i]] <- ClusProc0(signal=pX,Num=Num,prop=prop,threshold=threshold,thres_sil=thres_sil,thres_MAF=thres_MAF,cut=cut)
    }

  Nlen <- length(N)
  silR <- rep(NA,Nlen)

  for(i in 1:Nlen)
    if(adjust) silR[i] <- res[[i]]$sil$silMean_adjust else silR[i] <- res[[i]]$sil$silMean
  if(adjust) {
    N <- rep(NA,Nlen)
    for(i in 1:Nlen)
      N[i] <- res[[i]]$sil$clusNum_adjust
  }

  n <- which(silR==max(silR))
  clusNum <- N[n]
  logL <- res[[n]]$logL
  sil <- res[[n]]$sil

  resfinal <- list(clusNum=clusNum,logL=logL,silWidth=sil,signal=signal,adjust=adjust)
  class(resfinal) <- 'clust'
  return(resfinal)
  
}


print.clust <- function(resfinal) {
    adjust <- resfinal$adjust
    if(adjust)   res <- data.frame(clusNum_adjust=resfinal$silWidth$clusNum_adjust,logL=round(resfinal$logL,4),silMean_adjust=round(resfinal$silWidth$silMean_adjust,4)) else
    res <- data.frame(clusNum=resfinal$silWidth$clusNum,logL=round(resfinal$logL,4),silMean=round(resfinal$silWidth$silMean,4)) 
  print(res,quote=FALSE,row.names=FALSE)
}

