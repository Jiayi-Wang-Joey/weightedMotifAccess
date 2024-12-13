#' getBinData
#' 
#' break points into equal-width bins based on average value, merge bins that 
#' are below a min size, and compute mean and SD per sample per bin
#'
#' @param e A normalized matrix of log-transformed peak accessibility (e.g. 
#'   logCPM), with samples as columns.
#' @param minBinSize The minimum number of peaks in a bin
#' @param startingNBins The starting number of bins.
#'
#' @return A list
getBinData <- function(e, minBinSize=20, startingNBins=100){
    m <- rowMeans(e)
    b <- seq(min(m), max(m), length.out=1L+startingNBins)
    bins <- cut(m, breaks=b, include.lowest = TRUE)
    while(any((tt <- table(bins)) < minBinSize)){
        smallest_index <- which.min(tt)
        if(smallest_index==1){ # merge to the right
            merge_with <- 2L
        }else if(smallest_index==length(tt)){ # merge to the left
            merge_with <- length(tt)-1L
        }else{ # merge on the smallest side
            merge_with <- smallest_index + 
                c(-1,1)[which.min(tt[c(smallest_index-1,smallest_index+1)])]
        }
        levels(bins)[smallest_index] <- levels(bins)[merge_with]
    }
    binM <- mean(splitAsList(m,bins))
    binV <- aggregate(e,by=list(bin=bins),FUN=mean)[,-1]-binM
    binSd <- aggregate(e,by=list(bin=bins),FUN=sd)[,-1]
    list(binM=binM, be=binV, bsd=binSd, peak2bin=bins)
}

#' localWeightedSmooth
#' 
#' Apply local smoothing by doing a weighted average of +/- range
#'
#' @param x A numeric vector of values
#' @param weights A numeric vector of weights
#' @param range The smoothing range (+/- number of neighboring values)
#'
#' @return A smoothed `x`
localWeightedSmooth <- function(x, weights, range=1L){
    stopifnot(is.integer(range) && range>0)
    if(range==0) return(x)
    n <- seq_along(x)
    w2 <- c(rep(0,range),weights,rep(0,range))
    x2 <- c(rep(0,range),x,rep(0,range))
    totW <- rowSums(sapply(seq(from=-range,to=range), FUN=function(i){
        w2[range+i+n]
    }))
    wa <- x2*w2
    rowSums(sapply(seq(from=-range,to=range), FUN=function(i){
        wa[range+i+n]
    }))/totW
}

#' makeMonotonous
#' 
#' Makes a vector of values monotonous
#'
#' @param a A vector of numeric values
#'
#' @return A vector of numeric values
#' @examples
#' a <- c(1:5,3)
#' makeMonotonous(a)
makeMonotonous <- function(a){
    diffs <- sign(a[-1]-a[-length(a)])
    # if it's going down, reverse
    down <- sum(diffs) < 0
    if(down){
        diffs <- -diffs
        a <- -a
    }
    # find the longest monotonous stretch and split there
    rle <- Rle(diffs)
    wmr <- which.max(rle@lengths)
    offset <- 0L
    if(wmr>1) offset <- sum(rle@lengths[seq_len(wmr-1L)])
    splitAt <- offset+floor(rle@lengths[wmr]/2)
    a1 <- c(a[seq_len(splitAt)],Inf)
    a2 <- c(0,a[seq(from=splitAt+1L,length(a))])
    for(i in seq_along(a2)[-1]){ if(a2[i]<a2[i-1L]) a2[i] <- a2[i-1L] }
    for(i in rev(seq_along(a1))[-1]){ if(a1[i]>a1[i+1L]) a1[i] <- a1[i+1L] }
    a <- c(a1[-length(a1)], a2[-1])
    if(down) a <- -a
    a
}


#' binnedLoess
#'
#' @param e A normalized matrix of log-transformed peak accessibility (e.g. 
#'   logCPM), with samples as columns.
#' @param ... Passed to getBinData()
#'
#' @return A matrix of offsets (on the log-scale) per sample per peak
binnedLoess <- function(e, ..., 
    BPPARAM=BiocParallel::SerialParam(progress=TRUE)){
    ll <- getBinData(e, ...)
    m <- rowMeans(e)
    out <- BiocParallel::bplapply(colnames(e), BPPARAM=BPPARAM, \(i){
        # fit loess on bins
        mod <- loess(y~x, data=data.frame(x=as.matrix(ll$binM), y=ll$be[,i]),
            weights=1/(ll$bsd[,i]+1L), control=loess.control(surface="direct"))
        # predict on peaks
        predict(mod, newdata=data.frame(x=m))
    })
    matrix(unlist(out), ncol=ncol(ll$be), dimnames = dimnames(e))
}



#' monotonicSmoothed
#'
#' @param e A normalized matrix of log-transformed peak accessibility (e.g. 
#'   logCPM), with samples as columns.
#' @param ... Passed to getBinData()
#'
#' @return A matrix of offsets (on the log-scale) per sample per peak
monotonicSmoothed <- function(e, ..., range=2L){
    ll <- getBinData(e, ...)
    f <- function(m, a, w){
        a <- localWeightedSmooth(a, w, range=range)
        a <- makeMonotonous(a)
    }
    binOffset <- sapply(seq_len(ncol(ll$be)), \(i){
        f(ll$binM, ll$be[,i], w=1/(ll$bsd[,i]+1L))
    })
    m <- rowMeans(e)
    peakOffset <- apply(binOffset,2,\(a){
        y <- splinefun(ll$binM, a, method="monoH.FC")(m)
        y[which(m<min(ll$binM))] <- a[1]
        y[which(m>max(ll$binM))] <- a[length(a)]
        y
    })
    dimnames(peakOffset) <- dimnames(peakOffset)
    peakOffset
}

#' smoothedBinLm
#'
#' @param e A normalized matrix of log-transformed peak accessibility (e.g. 
#'   logCPM), with samples as columns.
#' @param ... Passed to getBinData()
#'
#' @return A matrix of offsets (on the log-scale) per sample per peak
smoothedBinLm <- function(e, ..., range=2L,
    BPPARAM=BiocParallel::SerialParam(progress=TRUE)){
    ll <- getBinData(e, ...)
    m <- rowMeans(e)
    out <- BiocParallel::bplapply(seq_len(ncol(ll$be)), BPPARAM=BPPARAM, \(i){
        # fit lm on bins
        co <- coef(lm.wfit(x=as.matrix(ll$binM), y=ll$be[,i],  w=1/(ll$bsd[,i]+1L)))
        # predict on peaks
        m*co
    })
    matrix(unlist(out), ncol=ncol(ll$be), dimnames = dimnames(e))
    
}