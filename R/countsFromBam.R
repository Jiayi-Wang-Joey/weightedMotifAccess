#' peakCountsFromBAM
#' 
#' Creates a SummarizedExperiment of fragment (or insertion) counts from bam 
#' files that overlap given regions.
#'
#' @param bam_files A vector of paths to the bam files.
#' @param regions A `GRanges` of regions in which to counts.
#' @param paired Logical; whether the data is paired (assumed unpaired by 
#'   default). Use `paired="auto"` for automatic detection using the first bam 
#'   file.
#' @param ignore.strand Logical; whether to ignore strand for the purpose of
#'   counting overlaps (default TRUE).
#' @param randomAcc Logical, whether to use random access. This is disabled by
#'   default because the overhead of random access to a lot of regions is 
#'   typically worse than reading the entire file. However, if you need to get
#'   counts in few regions, enabling this will be faster. Note however that 
#'   when using random access, the output object will not contain depth 
#'   information.
#' @param ov.type Overlap type. See the `type` argument of 
#'   \code{link[GenomicRanges]{countOverlaps}}.
#' @param maxgap Maximum gap allowed for overlaps (see the corresponding 
#'   argument of \code{link[GenomicRanges]{countOverlaps}}).
#' @param minoverlap Minimum overlap (see the corresponding argument of 
#'   \code{link[GenomicRanges]{countOverlaps}}).
#' @param getMedianFragLength Logical; whether to compile the median fragment
#'   length per region. This is slightly slower. The log10-transformed, 
#'   (weighted mean across samples of the) median fragment length per region is
#'   stored in `rowData(results)$flbias`.
#' @param BPPARAM BiocParallel params for multithreading. Note that
#'   multithreading can lead to high memory usage.
#' @param extend The amount *by* which to extend single-end reads (e.g. 
#'   fragment length minus read length). If `paired=TRUE` and `type` is either
#'   'ends' or 'center', then the extension will be applied after taking the 
#'   (shifted) fragment ends or centers, resulting in ranges of width equal to 
#'   `extend`.
#' @param type Type of the coverage to compile. Either full (full read/fragment),
#'   start (count read/fragment start locations), end, center, or 'ends' (both
#'   ends of the read/fragment).
#' @param strand Strand(s) to capture (any by default).
#' @param strandMode The strandMode of the data (whether the strand is given by
#'   the first or second mate, which depends on the library prep protocol). See
#'   \link[GenomicAlignments]{strandMode} for more information. This parameter 
#'   has no effect unless one of the `strand`, `extend` parameters or a 
#'   strand-specific `shift` are used.
#' @param shift Shift (from 3' to 5') by which reads/fragments will be shifted.
#'   If `shift` is an integer vector of length 2, the first value will represent
#'   the shift for the positive strand, and the second for the negative strand.
#' @param includeDuplicates Logical, whether to include reads flagged as 
#'   duplicates.
#' @param includeSecondary Logical; whether to include secondary alignments
#' @param minMapq Minimum mapping quality (1 to 255)
#' @param minFragLength Minimum fragment length (ignored if `paired=FALSE`)
#' @param maxFragLength Maximum fragment length (ignored if `paired=FALSE`)
#' @param splitByChr Whether to process chromosomes separately, and if so by how
#'   many chunks. The should not affect the output, and is simply slightly 
#'   slower and consumes less memory. Can be a logical value (in which case each
#'   chromosome is processed separately), but we instead recommend giving a 
#'   positive integer indicating the number of chunks.
#' @param verbose Logical; whether to print progress messages
#' @param ... Passed to `ScanBamParam`
#' @importFrom S4Vectors from to metadata countSubjectHits splitAsList
#' @importFrom BiocParallel bpnworkers SerialParam
#' 
#' @return A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} 
#'   with a 'counts' assay.
#' @export
#'
#' @examples
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create regions of interest
#' peaks <- GRanges(c("seq1","seq1","seq2"), IRanges(c(400,900,500), width=100))
#' peakCountsFromBAM(bam, peaks, paired=FALSE)
peakCountsFromBAM <- function(
                      bam_files, regions, paired, extend=0L, shift=0L, 
                      type=c("full","center","start","end","ends"),
                      ov.type="any", maxgap=-1L, minoverlap=1L,
                      ignore.strand=TRUE, strandMode=1, includeDuplicates=TRUE, 
                      includeSecondary=FALSE, minMapq=1L, minFragLength=1L,
                      maxFragLength=5000L, splitByChr=NULL, randomAcc=FALSE,
                      getMedianFragLength=FALSE, BPPARAM=SerialParam(),
                      verbose=TRUE, ...){
  # check inputs
  stopifnot(is(regions, "GRanges"))
  stopifnot(all(vapply(bam_files, FUN.VALUE=logical(1L), FUN=file.exists)))
  stopifnot(all(grepl("\\.bam",bam_files,ignore.case=TRUE)))
  type <- match.arg(type)
  stopifnot(extend>=0L)
  if(is.null(paired)){
    if(verbose) message("`paired` not specified, assuming single-end reads. ",
                        "Set to paired='auto' to automatically detect.")
    paired <- FALSE
  }else if(paired=="auto"){
    paired <- testPairedEndBam(head(bam_files[1]))
    if(verbose) message("Detected ", ifelse(paired,"paired","unpaired")," data")
  }
  if(!paired && !(type %in% c("center","ends")) && verbose && extend!=0)
    message("`extend` argument ignored.")
  if(type=="ends" && !paired){
    if(verbose)
      warning("type='ends' typically only makes sense with paired-end data...")
    type <- "start"
  }
  # prepare flags for bam reading
  flgs <- scanBamFlag(isDuplicate=ifelse(includeDuplicates,NA,FALSE), 
                      isSecondaryAlignment=ifelse(includeSecondary,NA,FALSE),
                      isNotPassingQualityControls=FALSE)
  seqs <- Rsamtools::scanBamHeader(bam_files[[1]])[[1]]$targets
  seqs <- seqs[.checkMissingSeqLevels(names(seqs), seqlevelsInUse(regions),
                                      argName="regions")]
  
  total_depth <- vapply(bam_files, FUN.VALUE=integer(1), FUN=function(x){
    as.integer(sum(idxstatsBam(BamFile(x, asMates=paired))$mapped))
  })
  
  if(is.null(randomAcc))
    randomAcc <- length(regions)<1000 & maxFragLength<=10000

  if(!randomAcc){
    if(is.null(splitByChr)){
      if(bpnworkers(SerialParam())==1){
        splitByChr <- 3L
      }else{
        splitByChr <- 8L
      }
    }
    param <- .getBamChunkParams(bam_files[[1]], flgs=flgs,
                                keepSeqLvls=names(seqs),  nChunks=splitByChr)
  }else{
    # resize and merge regions for random access to avoid double-counting
    sizeExt <- sum(abs(c(extend+shift+maxFragLength,1L)))
    regions2 <- reduce(resize(regions, width(regions)+sizeExt, fix = "center"))
    param <- list(x=ScanBamParam(flag=flgs, which=regions2, ...))
  }
  
  cnts <- bplapply(bam_files, BPPARAM=BPPARAM, FUN=function(bamfile){
    if(verbose && bpnworkers(BPPARAM)==1)
      message("Reading file ",bamfile)
    res <- lapply(names(param), FUN=function(x){
      r <- .bam2bwGetReads(bamfile, paired=paired, param=param[[x]], type=type,
                          extend=extend, shift=shift, minFragL=minFragLength,
                          maxFragL=maxFragLength, strandMode=strandMode,
                          si=seqs)
      if(getMedianFragLength){
        o <- findOverlaps(r, regions, type=ov.type, maxgap=maxgap,
                          minoverlap=minoverlap, ignore.strand=ignore.strand)
        mfl <- rep(0, length(regions))
        if(length(o) > 0){
          y <- IRanges::median(splitAsList(width(r)[queryHits(o)],
                                           subjectHits(o)))
          mfl[unique(sort(subjectHits(o)))] <- as.numeric(y)
        }
        return(list(ov=countSubjectHits(o), reads=metadata(r)$reads, mfl=mfl))
      }
      list(ov=countOverlaps(regions, r, type=ov.type, maxgap=maxgap,
                            ignore.strand=ignore.strand, minoverlap=minoverlap),
           reads=metadata(r)$reads)
    })
    depth <- sum(vapply(res, FUN.VALUE=integer(1), FUN=function(x) x$reads))
    mfl <- NULL
    if(getMedianFragLength) mfl <- Reduce("+", lapply(res, \(x) x$mfl))
    res <- rowSums(vapply(res, \(x) x[[1]], integer(length(regions))))
    gc(full=TRUE, verbose=FALSE)
    list(ov=res, depth=depth, mfl=mfl)
  })
  
  depths <- vapply(cnts, FUN.VALUE=integer(1), FUN=function(x) x$depth)
  
  if(getMedianFragLength)
    mfl <- matrix(unlist(lapply(cnts, \(x) x$mfl)), ncol=length(bam_files))
  
  cnts <- matrix(unlist(lapply(cnts, \(x) x[[1]])), ncol=length(bam_files))
  
  if(getMedianFragLength){
    mfl <- rowSums(mfl)/rowSums(cnts)
    regions$flbias <- log10(1+mfl)
  }
  
  if(is.null(names(bam_files))){
    if(!any(duplicated(bn<-basename(bam_files)))){
      colnames(cnts) <- gsub("\\.bam$","",bn,ignore.case=TRUE)
    }else if(!any(duplicated(bn<-dirname(bam_files)))){
      colnames(cnts) <- bn
    }else{
      colnames(cnts) <- bam_files
    }
  }else{
    colnames(cnts) <- names(bam_files)
  }
  if(is.null(names(regions))){
    names(regions) <- row.names(cnts) <- as.character(granges(regions))
  }else{
    row.names(cnts) <- names(regions)
  }
  se <- SummarizedExperiment(list(counts=cnts), rowRanges=regions)
  se$total_depth <- total_depth
  se$depth <- depths
  se
}
