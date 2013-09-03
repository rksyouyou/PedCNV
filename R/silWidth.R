silWidth <- function(dat,thres_sil=0.05,thres_MAF=0.01){

    old.order <- row.names(dat)
    obsNo <- dim(dat)[1]
    n <- dim(dat)[2]
    sdat <- dat[with(dat,order(dat[,n])),]
    clusRes <- sdat[,n]
    clusNo <- length(levels(factor(clusRes)))
    x <- as.matrix(sdat[,-n])
    t <- as.matrix(table(clusRes))
    breakp <- append(0,cumsum(t))
    thres_n <- ceiling(obsNo*thres_MAF)
    res <- .Call("sil_inter", x_ = x, clusRes_ = clusRes, clusNo_ = clusNo, obsNo_ = obsNo, breakp_= breakp)
    sil <- matrix(res$sil)
    rownames(sil) <- rownames(sdat)
    d0 <- res$d0
    rownames(d0) <- rownames(sdat)
    clus_adjust <- rep(NA,obsNo)
    for(i in 1:obsNo){
        if(abs(sil[i])<thres_sil) clus_adjust[i] <- -1
        else  clus_adjust[i] <- which(d0[i,]==min(d0[i,]))
    }
    clus_adjust_selec <- as.matrix(table(clus_adjust))
    abandon_clus <- -1
    for(i in 2:length(clus_adjust_selec))
        if(clus_adjust_selec[i,1]<thres_n) abandon_clus <- c(abandon_clus,rownames(clus_adjust_selec)[i])
    abandon_clus <- as.numeric(abandon_clus)
    sil_adjust <- matrix(NA,obsNo,1)
    rownames(sil_adjust) <- rownames(sil)
    for(i in 1:obsNo){
        if(clus_adjust[i]%in%abandon_clus)  sil_adjust[i,] <- NA
        else sil_adjust[i,] <- abs(sil[i,])
    }
    abandon <- row.names(sdat[clus_adjust==-1,])
    abandon_adjust <- row.names(sdat[clus_adjust%in%abandon_clus,])
    silRes <- data.frame(clus=clus_adjust,sil=sil_adjust)
    silRes_adjust<- silRes[!silRes$clus%in%abandon_clus,]

    ## change the order to the oringinal order
    old.order.abandon <- old.order[!old.order%in%abandon]
    silRes <- silRes[old.order.abandon,]
    old.order.abandon_adjust <- old.order[!old.order%in%abandon_adjust]
    silRes_adjust <- silRes[old.order.abandon_adjust,]

        
    clusAvg_adjust<- tapply(silRes_adjust$sil,silRes_adjust$clus,mean)
    silMean_adjust<- mean(silRes_adjust$sil)
    clusNum_adjust <- length(clusAvg_adjust)
    silRes0 <-data.frame(clus=clusRes,sil=sil) 
    clusAvg<- tapply(silRes0$sil,silRes0$clus,mean)
    silMean<- mean(silRes0$sil)
    clusNum <- length(clusAvg)
    

    return(list(clusNum_adjust=clusNum_adjust,silMean_adjust=silMean_adjust,clusAvg_adjust=clusAvg_adjust,silRes_adjust=silRes_adjust,abandon_adjust=abandon_adjust,clusNum=clusNum,silMean=silMean,clusAvg=clusAvg,silRes=silRes0,abandon=abandon))

}
