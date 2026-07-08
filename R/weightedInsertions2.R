#' weightedInsertions
#'
#' @param lf A named vector of paths to files, either (indexed) bam files or
#'   tabix-indexed fragment files.
#' @param peaks A GRanges of peaks.
#' @param motifMatches A GRanges of motif matches, with a motif_id metadata 
#'   column.
#' @param extension How much to extend the profile around the peak center.
#' @param shift Logical; whether to shift fragment ends. If not specified,
#'   shifting is performed for bam inputs.
#' @param smooth Smoothing span for lowess smoothing. If `smooth=TRUE`, Tukey's
#'   running median smoothing (3RS3R) is used.
#' @param ncores The number of threads to use, or a BiocParallelParam object.
#' @param verbose Logical; whether to print progress messages.
#' @param weightMode Either 'sub' (subtract a fraction of the minimum value in
#'   the profile before computing weights) or 'asIs' (use counts as is for 
#'   computing weights).
#' @param symm Logical; whether weights should be made symmetrical.
#' @param retFull Logical; whether to return the full matrices of counts around
#'   all matches.
#' @param res Eventual list of pre-computed counts around matches, bypassing 
#'   reading data. This is meant for debugging/development.
#'
#' @returns A list with the following slots:
#' \describe{
#' \item{\code{peakCounts}:}{A SummarizedExperiment object of insertion counts per peak.}
#' \item{\code{wInsCounts}:}{A matrix of weighted insertion counts per motif (rows) and sample (column).}
#' \item{\code{profiles}:}{A matrix of summed count at each position (row) around matches for each motif (column).}
#' \item{\code{weights}:}{A matrix of weights at each position (row) around matches for each motif (column).}
#' \item{\code{peakAnnotation}:}{A sparse logical matrix of motif (column) annotations for each peak (row).}
#' }
#' @export
#' @importFrom epiwraps bamChrChunkApply tabixChrApply views2Matrix
#' @importFrom BiocParallel bplapply bpstop MulticoreParam
#' @importFrom GenomicRanges sort reduce resize start end seqnames granges
#' @importFrom GenomicRanges countOverlaps overlapsAny coverage Views
#' @importFrom IRanges ranges
#' @importFrom stats smooth lowess
#' @examples
weightedInsertions <- function(lf, peaks, motifMatches, extension=200L,
                                 shift=NULL, ncores=1, verbose=TRUE,
                                 weightMode=c("sub","asIs"), smooth=1/30,
                                 symm=FALSE, retFull=FALSE, res=NULL){
  weightMode <- match.arg(weightMode)
  if(is.null(shift)){
    if(shift <- all(grepl("bam$", lf, ignore.case=TRUE)))
      message("`shift` argument not specified; as these appear ",
              "to be bam files, the fragments will be shifted.")
  }
  p2 <- sort(motifMatches)
  p2$motif_id <- factor(p2$motif_id)
  mow <- setNames(width(p2)[!duplicated(p2$motif_id)], unique(p2$motif_id))
  
  peaks <- sort(peaks)
  reduced <- reduce(c(peaks,p2))

  if(is.null(res)){
    
    if(verbose) message("Reading data and computing insertion coverage...")
    
    BP <- .wiNcores2BP(ncores, length(lf), progress=verbose)
    
    res <- bplapply(lf, BPPARAM=BP, FUN=\(f){
      if(grepl("bam$", f, ignore.case=TRUE)){
        rl <- epiwraps::bamChrChunkApply(f, paired=TRUE, progress=FALSE,
                                         FUN=.wiProcessChunk, shift=shift,
                                         extension=extension, p2=p2, peaks=peaks,
                                         reduced=reduced, BPPARAM=SerialParam())
      }else{
        rl <- epiwraps::tabixChrApply(f, FUN=.wiProcessChunk, shift=shift,
                                      extension=extension, p2=p2, peaks=peaks,
                                      reduced=reduced, BPPARAM=SerialParam())
      }
      rl <- rl[!sapply(rl, is.null)]
      bigmat <- do.call(rbind, lapply(rl, \(x) x$prof))
      peakCounts <- Reduce("+", lapply(rl, \(x) x$peakCounts))
      motif_id <- p2$motif_id[as.integer(row.names(bigmat))]
      depth <- Reduce("+", lapply(rl, \(x) x$nfrags))
      list( mat=bigmat, motif=motif_id, peakCounts=peakCounts, depth=depth )
    })
    bpstop(BP)
    
  }

  peakCounts <- vapply(res, FUN.VALUE=integer(length(peaks)),
                       FUN=\(x) x$peakCounts)
  peakCounts <- SummarizedExperiment(list(counts=peakCounts), rowRanges=peaks)
  peakCounts$depth <- vapply(res, FUN.VALUE=integer(1L), FUN=\(x) x$depth)
  
  if(verbose) message("Computing motif weight profiles")
  
  profiles <- .wiGetProfiles(res, extension, mow, norm=weightMode=="sub",
                             smooth.span=smooth, symm=symm)
  names(motifs) <- motifs <- colnames(profiles$weights)
  
  if(verbose) message("Computing weighted insertion counts")
  
  po <- to(distanceToNearest(p2, peaks))
  
  BP <- .wiNcores2BP(ncores, progress=verbose)
  wic <- matrix(unlist(bplapply(motifs, BPPARAM=BP, FUN=.wiGetWIC, res=res,
                                weights=profiles$weights, po=po)),
                ncol=length(res), byrow=TRUE)
  row.names(wic) <- motifs
  colnames(wic) <- names(res)
  
  if(retFull) profiles$fullRes <- res
  rm(res)
  
  pa <- sapply(split(p2, p2$motif_id), FUN=\(x) overlapsAny(peaks, x))
  pa <- as(pa, "sparseMatrix")
  rm(p2)
  gc()
  c(list(peakCounts=peakCounts, wInsCounts=wic, peakAnnotation=pa), profiles)
}

.wiNcores2BP <- function(ncores, max=Inf, ...){
  if(is.numeric(ncores) && ncores>=1 && round(ncores)==ncores){
    BP <- MulticoreParam(min(max,ncores), ...)
  }else if (inherits(ncores,"BiocParallelParam")){
    BP <- ncores
  }else{
    stop("Unrecognized value for the 'ncores' argument.")
  }
  BP
}

# get insertion coverages across motif matches
.wiProcessChunk <- function(x, peaks, p2, shift, reduced, extension=200L){
  chrs <- as.character(unique(seqnames(x)))
  w <- which(seqnames(p2) %in% chrs)
  if(length(w)==0 && !any(seqnames(peaks) %in% chrs)) return(NULL)
  frags <- length(x)
  if(shift) x <- resize(x, width(x)-8L, fix="center")
  x <- sort(GRanges(rep(seqnames(x),2), IRanges(c(start(x), end(x)), width=1L)))
  
  windows <- resize(granges(p2[w,]), width=2*extension+1L, fix="center")
  chrs <- intersect(as.character(unique(seqnames(windows))), chrs)

  x <- x[overlapsAny(x, reduced[seqnames(reduced) %in% chrs])]
  peakCounts <- countOverlaps(peaks, x, minoverlap=1L)
  w2 <- unlist(split(w,seqnames(windows),drop=FALSE)[chrs])
  
  v <- Views(coverage(x)[chrs],
             split(ranges(windows),seqnames(windows),drop=FALSE)[chrs])
  mat <- epiwraps:::views2Matrix(v, 0L)

  neg <- which(strand(p2[w2,])=="-")
  mat[neg,] <- mat[neg,rev(seq_len(ncol(mat)))]
  row.names(mat) <- w2
  rm(v, windows, x)
  gc(verbose=FALSE)
  list(prof=as(mat, "sparseMatrix"), peakCounts=peakCounts, nfrags=frags)
}

# based on insertion coverages across motif matches, get the summed profile and
# weight profile for each motif
.wiGetProfiles <- function(res, extension, mow, norm=TRUE, smooth.span=1/30, symm=FALSE){
  profiles <- Reduce("+", lapply(res, \(x){
    m <- t(rowsum(x$mat, x$motif))
    1000*m/sum(m)
  }))/length(res)
  meanProf <- rowMeans(profiles)
  weights <- vapply(setNames(colnames(profiles),colnames(profiles)),
              FUN.VALUE=vector("double", 2L*extension+1L),
              FUN=\(n){
    mhw <- floor(mow[n]/2)
    mhw <- seq(from=extension-mhw, to=extension+1L+mhw)
    x <- profiles[,n]
    x[mhw] <- wmed <- median(x[mhw])
    if(isTRUE(smooth.span)){
      x <- as.numeric(smooth(x))
    }else{
      x <- lowess(x, f=smooth.span)$y
    }
    x[mhw] <- pmin(mean(x),(wmed+mean(x[mhw]))/2)
    if(norm) x <- x-(1-abs(cor(meanProf, x)))*min(x)
    (1+2*extension)*x/sum(x)
  })
  if(symm) weights <- (weights+weights[rev(seq_len(nrow(weights))),])/2
  list(profiles=profiles, weights=weights)
}

# based on insertion coverages across motif matches and weight profiles, get 
# the sum weighted insertion counts for each sample/motif
.wiGetWIC <- function(mn, res, weights, po){
  wmo <- which(res[[1]]$motif==mn)
  unlist(lapply(res, \(x){
    fromPeak <- po[wmo]
    x <- x$mat[wmo,,drop=FALSE]
    x <- x %*% Diagonal(x=weights[,mn])
    o <- order(Matrix::rowSums(x), decreasing=TRUE)
    wmo2 <- o[!duplicated(fromPeak[o])]
    sum(x[wmo2,])
  }))
}
