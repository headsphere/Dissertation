# ---- jump ----

#Thanks to Gareth M James for making this code available
#http://www-bcf.usc.edu/~gareth/research/jump

jump <-
  function(data=NULL,K=10,y=NULL,plotjumps=T,rand=10,fits=NULL,B=0,dist=NULL,trace=F){
    if (!is.null(data)){
      # Compute the kmeans fit to the data
      if (is.null(fits))
        fits <- kmeans.rndstart(data,K,rand)
      if (is.null(y))
        y <- dim(data)[2]/2
      n <- nrow(data)
      p <- ncol(data)
      # Compute the distortion associated with the kmeans fit
      dist<- fits/(n*p)
    }
    # Call the compute.jump function to produce plots and calculate
    # maximum jump for each value of Y
    jump.results <- compute.jump(dist,y,plotjumps)
    jump.results$fits <- fits
    # Implement bootstrap routine
    if (B>0 & !is.null(data)){
      n <- nrow(data)
      boot.results <- matrix(0,length(y),K)
      bootdist <- rep(0,K)
      for (b in 1:B){
        if (trace)
          print(paste("Bootstrap Iteration ",b))
        # Make bootstrap data
        bootdata <- data[sample(1:n,replace=T),]
        # Get kmeans fit to the bootstrap data
        bootfits <- kmeans.rndstart(bootdata,K,rand)
        # Compute bootstrap distortion and maximum jumps
        for (k in 1:K)
          bootdist[k] <- sum(bootfits[[k]]$within)/(n*p)
        bootmaxjumps <- compute.jump(bootdist,y,plotjumps=F,printresults=F)$maxjump
        for (j in 1:length(y))
          boot.results[j,bootmaxjumps[j]] <-  boot.results[j,bootmaxjumps[j]]+1
      }
      # Calculate proportions of each number of clusters chosen
      jump.results$boot.result <- round(boot.results/B,3)
      for (j in 1:length(y))
        print(paste(jump.results$boot.result[j,jump.results$maxjump[j]]*100,"% of bootstrap iterations corresponding to ",jump.results$maxjump[j], "clusters with Y=",y[j]))
    }
    jump.results}

compute.jump <-
  function(dist,y,plotjumps=T,printresults=T){
    K <- length(dist)
    numb.y <- length(y)
    numbclust <- rep(0,numb.y)
    transdist <- matrix(0,numb.y,K+1)
    jumps <- matrix(0,numb.y,K)
    if (plotjumps)
      par(mfrow=c(numb.y,3))
    for (i in 1:numb.y){
      # Compute the transformed distortion
      transdist[i,] <- c(0,dist^(-y[i]))
      # Compute the jumps in transformed distortion
      jumps[i,] <- diff(transdist[i,])
      # Compute the maximum jump
      numbclust[i] <- order(-jumps[i,])[1]
      # Plot distortion, transformed distortion and jumps
      if (plotjumps){
        plot(1:K,dist,type='l',
             xlab="Number of Clusters",ylab="Distortion",main=paste("Y = ",y[i]))
        plot(0:K,transdist[i,],type='l',
             xlab="Number of Clusters",ylab="Transformed Distortion",main=paste("Y = ",y[i]))
        plot(1:K,jumps[i,],type='l',
             xlab="Number of Clusters",ylab="Jumps",main=paste("Y = ",y[i]))
        # Plot line and point to indicate maximum jump
        lines(rep(numbclust[i],2),c(0,jumps[i,numbclust[i]]),lty=3,lwd=3)
        points(numbclust[i],jumps[i,numbclust[i]],col=2,pch=19,cex=1.5)}
      # Report maximum jump
      if (printresults)
        print(paste("The maximum jump occurred at ",numbclust[i], "clusters with Y=",y[i]))
    }
    list(maxjump=numbclust,dist=dist,transdist=transdist[,-1],jumps=jumps)}

kmeans.rndstart <-
  function(x, K, rand = 10)
  {
    fits <- sum((t(x) - apply(x, 2, mean))^2)
    iter.max <- 10
    # Run kmeans for 2 to K clusters
    for (k in 2:K){
      Z=kmeans(x,k,nstart=rand)
      fits=c(fits,sum(Z$withinss))
    }
    fits
  }

# testdata <- matrix(c(rnorm(2000),3+rnorm(2000)),byrow=T,ncol=4)
# temp <- jump(testdata,y=c(1.5,2,2.5),rand=10,trace=F)