ClusProc0 <- function(signal,Num,threshold,thresMAF,cut,itermax){

    thres_sil <- 0.01
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
        if(iter>itermax) break
        logLold <- logL
    }
    
    logL <- logL-1/2*S*cut*log(2*pi)
    print(paste("The logliklihood is",logL,"when clustering number is ",Num))
    
    sil <- silWidth(Sdata,thres_sil=thres_sil,thres_MAF=thresMAF)
    return(list(logL=logL,sil=sil)) 
}

##' This function chooses the optimal number of clusters and provides the assignments of each individuals under the optimum clustering number.
##' @title CNV clustering Procedure
##' @param signal The matrix of intensity measurements. The row names must be consistent with the Individual ID in fam file.
##' @param N Number of clusters one wants to fit to the data. N needs to be larger than 1 and if it is 1, error will be returned. The default value 2,3,...,6 will be used if it is missing.
##' @param  varSelection Factor. For specifying how to handle the intensity values. It must take value on 'RAW', 'PC.9', 'PC1'and 'MEAN'. If the value is 'RAW', then the raw intensity value will be used. If it is 'PC.9', then the first several PCA scores which account for 90\% of all the variance will be used. If the value is 'PC1', then the first PCA scores will be used. If the value is 'MEAN', the mean of all the probes will be used. 
##' @param threshold Optional number of convergence threshold. The iteration stops if the absolute difference of log likelihood between successive iterations is less than it. The default threshold 1e-05 will be used if it's missing.
##' @param itermax Optional. The iteration stops if the time of iteration is large than this value. The default number 8 will be used if it's missing.
##' @param adjust Logicals, If TRUE (default), the result will be adjusted by the silhouette score. See details. 
##' @param thresMAF The minor allele frequency threshold.
##' P-values are not calculated for CNVs of which minor allele frequencies are less than thresMAF.
##' @param scale Logicals. If TRUE, the signal will be scale by using sample mean and sample variance by columns before further data-processing.
##' @details 
##' \itemize{
##' \item{adjust}{If adjust is TRUE, the result will be adjusted by the silhouette score in the following criterion. For each individual, the silhouette scores are calculated for each group. The individual will assigned forcefully to the group which maximize the silhouette scores. }
##' }
##' @return It returns object of class 'clust'. 'clust' is a list containing following components:
##' \item{clusNum}{The optimal number of clusters among give parameter {N}.}
##' \item{silWidth}{Silhouette related results.}
##' @author Meiling Liu 
##' @examples
##' # Load data and correlation matrix
##' data(simudat)
##' # Fit the data under the given clustering numbers
##' clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
##' @export
ClusProc <- function(signal,N=2:6,varSelection=c('RAW','PC.9','PC1','MEAN'),threshold=1e-05,itermax=8,adjust=TRUE,thresMAF=0.01,scale=FALSE){

    sX0 <- as.matrix(signal)
    if(scale) sX <- scale(sX0) else sX <- sX0
    cut <- 1
    if(varSelection=='PC.9'){
        prop <- 0.9
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
    if(varSelection=='PC1'){
        pX <- prcomp(sX)$x[,1]
    }
  
    if(1%in%N) stop('The assigned clustering number must be larger than 1!')
    res <- list()
    for(i in 1:length(N)){
        Num <- N[i]
        res[[i]] <- ClusProc0(signal=pX,Num=Num,threshold=threshold,thresMAF=thresMAF,cut=cut,itermax=itermax)
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
    
    resfinal <- list(clusNum=clusNum,silWidth=sil,signal=signal,adjust=adjust)
    class(resfinal) <- 'clust'
    return(resfinal)
  
}

##' Prints formatted results returned by \code{\link{ClusProc}}.
##'
##' @title Prints clustering results
##' @param x The clustering results obtained from the \code{\link{ClusProc}}.
##' @param ... Usual arguments passed to the print function.
##' @author Meiling Liu and Sungho Won
##' @method print clust
##' @examples
##' # Load data and correlation matrix
##' data(simudat)
##' # Fit the data under the given clustering numbers
##' clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
##' print(clus.fit)
##' @export
print.clust <- function(x, ...) {
    adjust <- x$adjust
    if(adjust)   res <- data.frame(clusNum_adjust=x$silWidth$clusNum_adjust,silMean_adjust=round(x$silWidth$silMean_adjust,4)) else
    res <- data.frame(clusNum=x$silWidth$clusNum,silMean=round(x$silWidth$silMean,4)) 
    print(res,quote=FALSE,row.names=FALSE)
}

