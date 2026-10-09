#' bamChrChunkApply
#' 
#' Runs a function on reads/fragments from chunks of (chromosomes) of an 
#' indexed bam file. This is especially used by other functions to avoid 
#' loading all alignments into memory, or to parallelize reads processing.
#'
#' @param x A bam file.
#' @param FUN The function to be run, the first argument of which should be a
#'   `GRanges`
#' @param paired Logical; whether to consider the reads as paired (fragments, 
#'   rather than reads, will be returned)
#' @param strandMode Strandmode for paired data (see 
#'   \code{\link[GenomicAlignments]{readGAlignmentPairs}}).
#' @param keepSeqLvls An optional vector of seqLevels to keep
#' @param flgs `scanBamFlag` for filtering the reads
#' @param mapqFilter Integer of the minimum mapping quality for reads to be 
#'   included.
#' @param nChunks The number of chunks to use (higher will use less memory but
#'   increase overhead)
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that if used, memory usage will be high; in this context we recommend a 
#'   high `nChunks`.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param progress Logical; whether to show a progress bar.
#' @param ... Passed to `FUN`
#'
#' @return A list of whatever `FUN` returns
#' @export
#' @examples
#' # as an example we'll use the function to obtain fragment sizes:
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' fragLen <- bamChrChunkApply(bam, paired=TRUE, FUN=function(x) width(x))
#' quantile(unlist(fragLen))
bamChrChunkApply <- function(x, FUN, paired=FALSE, keepSeqLvls=NULL, nChunks=4,
                             strandMode=2, flgs=scanBamFlag(), exclude=NULL,
                             mapqFilter=NA_integer_, progress=TRUE,
                             BPPARAM=NULL, ...){
  if(!is.null(exclude)) stopifnot(is(exclude, "GRanges"))
  param <- .getBamChunkParams(x, flgs=flgs, keepSeqLvls=keepSeqLvls, 
                              nChunks=nChunks)
  f2 <- function(p, ...){
    if(paired){
      bam <- .initPairedBam(x)
      x <- readGAlignmentPairs(bam, param=p, strandMode=strandMode)
      x <- as(x[isProperPair(x)], "GRanges")
    }else{
      x <- GRanges(readGAlignments(x, param=p))
    }
    if(!is.null(exclude)) x <- x[!overlapsAny(x,exclude)]
    if(paired && length(x)==0)
      warning("Nothing found (in one of the chunks). If this is unexpected, it",
              " could be because your read mates don't have matching names, or",
              " have suffixes to their names. If this is the case, specify it ",
              'by using an input like:
  BamFile("aligned/test.bam", asMates=TRUE, qnameSuffixStart="/")',
              immediate.=TRUE)
    x <- FUN(x, ...)
    gc(full=TRUE, verbose=FALSE)
    x
  }
  if(is.null(BPPARAM) || BiocParallel::bpnworkers(BPPARAM)==1){
    if(progress) return(pblapply(param, FUN=f2, ...))
    return(lapply(param, FUN=f2, ...))
  }
  
  bplapply(param, FUN=f2, ..., BPPARAM=BPPARAM)
}

.getBamChunkParams <- function(x, flgs=scanBamFlag(), keepSeqLvls=NULL, 
                               nChunks=4, ...){
  seqs <- Rsamtools::scanBamHeader(x)
  if(is.null(seqs$targets)) seqs <- seqs[[1]]
  seqs <- seqs$targets
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls, names(seqs)))>0){
      stop(paste0(
        "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
head(paste(missingLvls, collapse=", "), 3)))
    }
    seqs <- seqs[keepSeqLvls]
  }
  # generate list of reading params
  if(isTRUE(nChunks) || nChunks>=2){
    if(is.numeric(nChunks) && nChunks>=2){
      stopifnot(round(nChunks)==nChunks)
      nChunks <- as.integer(nChunks)
      seqs <- sort(seqs, decreasing=TRUE)
      chrGroup <- head(rep(seq_len(nChunks), 
                           ceiling(length(seqs)/nChunks)), length(seqs))
      chrGroup <- split(seqs, chrGroup)
    }else{
      chrGroup <- split(seqs, names(seqs))
    }
    param <- lapply(chrGroup, FUN=function(x)
      ScanBamParam(flag=flgs, ..., which=GRanges(names(x), IRanges(1L,x))))
  }else{
    chrGroup <- list(wg=seqs)
    param <- list(wg=ScanBamParam(flag=flgs, ..., 
                                  which=GRanges(names(seqs), IRanges(1L,seqs))))
  }
  param
}

.initPairedBam <- function(path){
  bam <- BamFile(path, asMates=TRUE, yieldSize=10L)
  qn <- Rsamtools::scanBam(bam, param=ScanBamParam(what="qname"))[[1]]$qname
  suffix <- NA
  if(sum(grepl("/[0-3]$", qn))>=9) suffix <- "/"
  if(sum(grepl("\\.[0-3]$", qn))>=9) suffix <- "."
  BamFile(path, asMates=TRUE, qnameSuffixStart=suffix)
}


.align2cuts <- function(x, size=1L){
  sort(c(resize(x,size,fix="start"), resize(x,size,fix="end")))
}

#' tabixChrApply
#' 
#' Runs a function on reads/fragments from each chromosomes of a Tabix-indexed 
#' fragment file. This is especially used by other functions to 
#' avoid loading all alignments into memory, or to parallelize reads processing.
#'
#' @param x The path to a tabix-indexed bam file, or a TabixFile object.
#' @param fn The function to be run, the first argument of which should be a
#'   `GRanges`
#' @param keepSeqLvls An optional vector of seqLevels to keep
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that if used, memory usage will be high; in this context we recommend a 
#'   high `nChunks`.
#' @param only An optional GRanges of regions for which overlapping reads should
#'   be included. If set, all other reads are discarded.
#' @param progress Logical; whether to show a progress bar.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param ... Passed to `fn`
#'
#' @return A list of whatever `fn` returns
#' @export
#' @importFrom Rsamtools TabixFile seqnamesTabix bgzip indexTabix
#' @importFrom rtracklayer path import
#' @importFrom BiocParallel bpnworkers bplapply
#' @importFrom pbapply pblapply
#' @examples
#' # generate dummy regions and save them to a temp file:
#' frags <- tempfile(fileext = ".tsv")
#' d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
#' d$end <- d$start + 15L
#' write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t")
#' # tabix-index it
#' frags <- Rsamtools::bgzip(frags)
#' Rsamtools::indexTabix(frags, format = "bed")
#' # now we can do something chunk-wise, e.g. extract coverage:
#' res <- tabixChrApply(frags, fn=coverage)
#' # aggregate the chunk results into an RleList object:
#' reduceRleLists(res)
tabixChrApply <- function(x, fn, keepSeqLvls=NULL, exclude=NULL, only=NULL,
                          BPPARAM=NULL, progress=TRUE, ...){
  x <- TabixFile(x)
  if(!is.null(exclude)) stopifnot(is(exclude, "GRanges"))
  if(!is.null(only)) stopifnot(is(only, "GRanges"))
  
  f2 <- function(sn, postfn, ...){
    x <- rtracklayer::import(path(x), format="bed",
                             which=GRanges(sn, IRanges(1,5*10^8)))
    if(!is.null(only)){
      .comparableStyles(x, only)
      x <- x[overlapsAny(x, only)]
    }
    if(!is.null(exclude)){
      .comparableStyles(x, exclude)
      x <- x[!overlapsAny(x, exclude)]
    }
    postfn(x, ...)
  }
  
  slvls <- seqnamesTabix(x)
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls, slvls)>0))
      stop(paste0(
        "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
        head(paste(missingLvls, collapse=", "), 3)))
    slvls <- keepSeqLvls
  }
  
  if(is.null(BPPARAM) || BiocParallel::bpnworkers(BPPARAM)==1){
    if(progress) return(pblapply(slvls, FUN=f2, postfn=fn, ...))
    return(lapply(slvls, FUN=f2, postfn=fn, ...))
  }
  bplapply(slvls, FUN=f2, postfn=fn, ..., BPPARAM=BPPARAM)
}

.comparableStyles <- function(a,b,stopIfNot=TRUE){
  if(all(.getSeqLevelsStyle(a)==.getSeqLevelsStyle(b))) return(TRUE)
  msg <- paste("It seems your are providing objects for which the seqlevel ",
               "styles do not match.")
  if(stopIfNot) stop(msg)
  warning(msg)
  FALSE
}
.getSeqLevelsStyle <- function(x){
  if(is.character(x)){
    if(grepl("\\.bw$|\\.bigwig", x, ignore.case=TRUE)){
      return(seqlevelsStyle(rtracklayer::BigWigFile(x)))
    }else if(grepl("\\.bam$", x, ignore.case=TRUE)){
      return(seqlevelsStyle(Rsamtools::BamFile(x)))
    }
    return(tryCatch({
      seqlevelsStyle(rtracklayer::BEDFile(x))
    }, error=function(e)
      stop("Unknown filetype.")))
  }
  seqlevelsStyle(x)
}