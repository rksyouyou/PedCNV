plot.clust <- function(resfinal,type=c('histo','scat','sil'),save=c(TRUE,FALSE),adjust=TRUE){

    if(type=='histo'){
        
        sX <- as.matrix(resfinal$signal)
        pca <- princomp(sX)
        PCA1 <- sX%*%loadings(pca)[,1]

        if(adjust) {
            clusters <- matrix(resfinal$silWidth$silRes_adjust$clus)
            rownames(clusters) <- rownames(resfinal$sil$silRes_adjust)} else {
                clusters <- matrix(resfinal$silWidth$silRes$clus)
                rownames(clusters) <- rownames(resfinal$sil$silRes)}
        temp <- merge(PCA1,clusters,by='row.names')[,-1]
        colnames(temp) <- c('PCA1','clusters')
        temp[,2] <- factor(temp[,2])
        p <- qplot(PCA1,fill=clusters,data=temp,geom="density")+geom_histogram(aes(y=..count..),binwidth=0.2)
        if(save) ggsave(paste0("histo_Comp1_N",resfinal$clusNum,".png"))
  }
    
  if(type=='scat'){

      signal <- resfinal$signal
      sX <- as.matrix((signal))
      pca <- princomp(sX)
      PCA1 <- sX%*%loadings(pca)[,1]
      segmean <- as.matrix(apply(sX,1,mean))

      if(adjust) {
          clusters <- matrix(resfinal$silWidth$silRes_adjust$clus)
          rownames(clusters) <- rownames(resfinal$sil$silRes_adjust)} else {
              clusters <- matrix(resfinal$silWidth$silRes$clus)
              rownames(clusters) <- rownames(resfinal$sil$silRes)}
      

      temp <- merge(cbind(segmean,PCA1),clusters,by='row.names')[,-1]
      colnames(temp) <- c('Mean','PCA1','clusters')
      temp[,3] <- factor(temp[,3])
      p <- qplot(Mean,PCA1,color=clusters,data=temp)
      if(save) ggsave(paste0("scat_Comp1_N",resfinal$clusNum,".png"))
  }

  if(type=='sil'){

      if(adjust){
          silRes <- resfinal$silWidth$silRes_adjust
          silMean <- resfinal$silWidth$silMean_adjust
          clusNo <- resfinal$silWidth$clusNum_adjust
          clusAvg <- resfinal$silWidth$clusAvg_adjust
          abandon_num <- length(resfinal$silWidth$abandon_adjust)
      }else{
          silRes <- resfinal$silWidth$silRes
          silMean <- resfinal$silWidth$silMean
          clusNo <- resfinal$silWidth$clusNum
          clusAvg <- resfinal$silWidth$clusAvg
          abandon_num <- length(resfinal$silWidth$abandon)
      }
      obsNo <- dim(silRes)[1]
      clusRes <- silRes$clus

      if(save) pdf(paste0("silWidthN",clusNo,".pdf"))
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
      if(save) dev.off()

  }
}
