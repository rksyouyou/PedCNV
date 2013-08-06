##' Makes formatted plots from the clustering result returned from ClusProc.
##'
##' @title Plots clustering result
##' @param x The clustering results obtained from ClusProc.
##' @param type Factor. For specifying the plot type. It must take value on 'histo', 'scat' and 'sil'. If the value is 'histo', then the histogram of the first PC score of the intensity measurement will be made. If the value is 'scat', then the scatter plot of the first PC score of the intensity measurement v.s the mean of the intensity measurement will be made. If the value is 'sil', then the silhouette plot will be made.
##' @param adjust Logicals. If TRUE (default), the silhouette-adjusted clustering result will be used. If FALSE, the initial clustering result will be used.
##' @param ... Usual arguments passed to the qplot function.
##' @author Meiling Liu
##' @method plot clust
##' @examples
##' # Load data and correlation matrix
##' # data(simudat)
##' # Fit the data under the given clustering numbers
##' # clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
##' # plot(clus.fit,type='histo')
##' @export
plot.clust <- function(x,type=c('histo','scat','sil'), adjust=TRUE, ...){

    if(type=='histo'){
        
        sX <- as.matrix(x$signal)
        pca <- princomp(sX)
        PCA1 <- sX%*%loadings(pca)[,1]

        if(adjust) {
            clusters <- matrix(x$silWidth$silRes_adjust$clus)
            rownames(clusters) <- rownames(x$sil$silRes_adjust)} else {
                clusters <- matrix(x$silWidth$silRes$clus)
                rownames(clusters) <- rownames(x$sil$silRes)}
        temp <- merge(PCA1,clusters,by='row.names')[,-1]
        colnames(temp) <- c('PCA1','clusters')
        temp[,2] <- factor(temp[,2])
        print(qplot(PCA1,fill=clusters,data=temp,geom="density")+geom_histogram(aes(y=..count..),binwidth=0.2))
  }
    
  if(type=='scat'){

      signal <- x$signal
      sX <- as.matrix((signal))
      pca <- princomp(sX)
      PCA1 <- sX%*%loadings(pca)[,1]
      segmean <- as.matrix(apply(sX,1,mean))

      if(adjust) {
          clusters <- matrix(x$silWidth$silRes_adjust$clus)
          rownames(clusters) <- rownames(x$sil$silRes_adjust)} else {
              clusters <- matrix(x$silWidth$silRes$clus)
              rownames(clusters) <- rownames(x$sil$silRes)}
      

      temp <- merge(cbind(segmean,PCA1),clusters,by='row.names')[,-1]
      colnames(temp) <- c('Mean','PCA1','clusters')
      temp[,3] <- factor(temp[,3])
      print(qplot(Mean,PCA1,color=clusters,data=temp))
  }

  if(type=='sil'){

      if(adjust){
          silRes <- x$silWidth$silRes_adjust
          silMean <- x$silWidth$silMean_adjust
          clusNo <- x$silWidth$clusNum_adjust
          clusAvg <- x$silWidth$clusAvg_adjust
          abandon_num <- length(x$silWidth$abandon_adjust)
      }else{
          silRes <- x$silWidth$silRes
          silMean <- x$silWidth$silMean
          clusNo <- x$silWidth$clusNum
          clusAvg <- x$silWidth$clusAvg
          abandon_num <- length(x$silWidth$abandon)
      }
      obsNo <- dim(silRes)[1]
      clusRes <- silRes$clus

      s <- rev(silRes[,"sil"])
      space <- c(0,rev(diff(cli <- silRes$clusRes)))
      space[space!=0] <- 5
      xlab <- expression("Silhouette width"* s[i])
      main <- paste("Silhouette plot")
      sub <- paste("Average silhouette width:",round(silMean,4))
      y <- barplot(s,width=1,space=space,xlim=c(min(0,min(s)),1),horiz=TRUE,col="grey",mgp=c(2.5,1,0),las=1,border=0,xlab=xlab)
      title(main=main,sub=sub,adj=0)
      mtext(paste("n=",obsNo,'; abandon=', abandon_num),adj=0)
      mtext(substitute(k ~ ~"clusters" ~ ~C[j], list(k=clusNo)),adj=1)
      mtext(expression(paste(j, " :  ", n[j], " | ", ave[i %in% Cj] ~ ~s[i])), adj = 1.04, line = -1.2)
      y <- rbind(rev(y),(clusRes))
      for (j in 1:clusNo) {
          yj <- mean(y[1,y[2,]==j])
          text(1,yj , paste(j, ":  ",table(clusRes)[j], " | ", format(clusAvg[j], digits = 1, nsmall = 2)), xpd = NA, adj = 0.8)
      }
  }
}
