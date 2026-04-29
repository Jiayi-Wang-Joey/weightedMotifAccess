#' @title .resizeRanges
#'
#' @description
#' Resize the peaks
#'
#' @param peakRanges: a GRanges object of peak ranges
#' @param width: the re-defined size of each peak
#' @param fix: the fixed point for resizing
#' @return a GRange object with resized ranges
#' @author Jiayi Wang

.resizeRanges <- function(peakRanges,
    width = 300,
    fix = c("center", "start", "end", "summit"),
    ...) {

    fix <- match.arg(fix, choices = c("center", "start", "end", "summit"))
    # Sanity check
    if (!class(peakRanges) == "GRanges") {
        stop("peakRanges must be a GRanges object")
    }

    if (fix == "summit") {
        start(peakRanges) <- round(peakRanges$summit-width/2)
        end(peakRanges) <- start(peakRanges)+width-1
    } else {
        peakRanges <- resize(peakRanges, width = width, fix = fix)
    }

    return(peakRanges)
}

#' @title .sanityCheck
#'
#' @description
#' Sanity check to ensure the input arguments have the correct classes
#'
#' @param atacFrag: a list of data.tables that contain the ranges of fragments
#' @param peakRanges: a GRange object that contains the ranges of peaks
#' @param motifRanges:  a GRange object that contains the ranges of motifs,
#' check metadata columns
#' check seqnames to factor in datatable
.sanityCheck <- function(atacFrag,
                         ranges,
                         type = c("peaks", "motifs")) {
    lapply(atacFrag, function(frag) {
        if (!is.data.table(frag)) {
            stop("Each element in atacFrag should be a data.table")
        }
    })

    if (!class(ranges)=="GRanges") {
        stop("ranges should be a GRanges object")
    }

    type <- match.arg(type, choices = c("peaks", "motifs"))
    if (type=="motifs") {
        if (!("motif" %in% names(mcols(ranges)))) {
            stop("There is no motif names in metadata columns")
        }
    }

}
#' @title getGCcontent
#'
#' @description To calculate the GC content of each fragment or peak
#' @param gr: a GRanges object
#' @param genome: a BSgenome object, the corresponding genome
#' @return a GRanges objects with an additional metadata column gc that contains
#' GC content
#' @author  Jiayi Wang

.getGCContent <- function(gr, genome, contentOnly=FALSE) {
    # Sanity check
    if(is.data.table(gr) | is.data.frame(gr)){
      gr <- makeGRangesFromDataFrame(as.data.frame(gr))
    }

    if (!class(gr) == "GRanges") {
        stop("peakRanges must be a GRanges object")
    }
    seqs <- Biostrings::getSeq(x = genome, gr)
    mcols(gr)$gc <- letterFrequency(seqs, "GC",as.prob=TRUE)[,1]
    if(contentOnly){
      return(gr$gc)}
    else{
      return(gr)
    }
}

#' @title filterFrags
#'
#' @description remove fragments that are too short or too long
#' @param atacFrag: a list of data tables containing the fragments information
#' @param min and @param max the minimum and maximum limit of fragment length
#' @return a list of data tables of filtered fragments
#' @author Jiayi Wang
.filterFrags <- function(atacFrag, min = 30, max = 2000) {
    res <- lapply(atacFrag, function(frag) {
        frag[,width:=end-start+1]
        frag[width>=min & width<=max,]
        frag
    })
    res
}

#' @title .matchSeqlevels
#'
#' @description match the chromosomes between fragments and peaks/motifs
#' @param atacFrag: a list of data tables containing the fragments information
#' @param ranges: a genomic object of peaks/motifs
#' @author Jiayi Wang

.matchSeqlevels <- function(atacFrag, ranges) {
    frags <- rbindlist(atacFrag)
    fragSeq <- unique(frags$seqnames)
    rangeSeq <- GenomicRanges::seqnames(ranges)
    common <- intersect(fragSeq,rangeSeq)
    atacFrag <- lapply(atacFrag, function(frag) {
        frag <- frag[seqnames %in% common,]
        frag$seqnames <- factor(frag$seqnames)
        frag})
    ranges <- ranges[seqnames(ranges) %in% common,]
    # turn seqnames to factor
    list(atacFrag=atacFrag, ranges=ranges)
}

#' @title dtToGr
#'
#' @description
#' Convert a data table to GenomicRange object
#' @author Emanuel Sonder

.dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end",
                    strandCol="strand", stranded=FALSE, addMetaCols=TRUE){
  dt <- copy(dt)
  setnames(dt, seqCol, "seqnames", skip_absent = TRUE)

  if(stranded) strand <- dt[[strandCol]] else strand <- NULL

  gr <- GRanges(seqnames=dt[["seqnames"]],
                strand=strand,
                ranges=IRanges(start=dt[[startCol]], end=dt[[endCol]]))

  if(startCol==endCol)
  {
    gr <- GPos(seqnames=dt[["seqnames"]],
               strand=strand,
               pos=dt[[startCol]])
  }

  if(addMetaCols){
    metaCols <- dt[,setdiff(colnames(dt),
                            c(seqCol, startCol, endCol, strandCol,
                              "seqnames", "chr")),with=FALSE]
    mcols(gr) <- metaCols
  }

  return(gr)
}

.standardChromosomes <- function(gr, species) {
    gr <- keepStandardChromosomes(gr,
        species=species,
        pruning.mode="coarse")
    seqlevelsStyle(gr) <- "UCSC"
    gr
}

.processData <- function(data, readAll=FALSE, shift=FALSE,
                         subSample=NULL, seqLevelStyle="UCSC",
                         minFrag=30, maxFrag=3000){
  if(is.character(data)){
    if(grepl(".bam", basename(data), fixed=TRUE))
    {
      param <- Rsamtools::ScanBamParam(what=c('pos', 'qwidth', 'isize'))
      readPairs <- GenomicAlignments::readGAlignmentPairs(data, param=param)

      # get fragment coordinates from read pairs
      seqDat <- GRanges(seqnames(readPairs@first),
                        IRanges(start=pmin(GenomicAlignments::start(readPairs@first),
                                           GenomicAlignments::start(readPairs@last)),
                                end=pmax(GenomicAlignments::end(readPairs@first),
                                         GenomicAlignments::end(readPairs@last))),
                        strand=GenomicAlignments::strand(readPairs))
      seqDat <- granges(seqDat, use.mcols=TRUE)
      seqDat <- as.data.table(seqDat)
      setnames(seqDat, c("seqnames"), c("chr"))
    }
    else if(grepl(".bed", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data, stringsAsFactors=TRUE)
      else{

        readBed <- function(data){
          tryCatch(
            {
              seqDat <- fread(data, select=c(1:3,6),
                              col.names=c("chr", "start", "end", "strand"),
                              stringsAsFactors=TRUE)
              return(seqDat)},
            error = function(cond){
              seqDat <- fread(data, select=c(1:3),
                              col.names=c("chr", "start", "end"),
                              stringsAsFactors=TRUE)
              return(seqDat)
            })}
        seqDat <- readBed(data)
      }
    }
    else if(grepl(".tsv", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data, stringsAsFactors=TRUE)
      else{
        seqDat <- fread(data, select=c(1:3),
                        col.names=c("chr", "start", "end"),
                        stringsAsFactors=TRUE)}
      if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    }
    else if(grepl(".rds", basename(data), fixed=TRUE)){
      seqDat <- as.data.table(readRDS(data))
      if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    }
  }
  else{
    seqDat <- as.data.table(data)
    if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    seqDat$chr <- factor(seqDat$chr)
  }

  if(!is.null(subSample) & is.numeric(subSample)){
    message("Subsampling file")
    subSample <- as.integer(subSample)
    seqDat <- seqDat[sample(1:nrow(seqDat), min(nrow(seqDat), subSample)),]
  }

  # Match seqlevelstyle to reference
  if((sum(grepl("chr", levels(seqDat$chr)))==0 & seqLevelStyle=="UCSC") |
     (sum(grepl("chr", levels(seqDat$chr)))>0 & seqLevelStyle=="NCBI")){
    tmpgr <- GRanges(levels(seqDat[["chr"]]),
                      IRanges(seq_along(levels(seqDat[["chr"]])), width=2L))
    seqlevelsStyle(tmpgr) <- seqLevelStyle
    levels(seqDat$chr) <- seqlevels(tmpgr)
  }

  # Insert ATAC shift
  if(shift){
    seqDat[, start:=start+4L]
    seqDat[, end:=end-4L]
  }
  else if(shift){
    warning("Did not shift as no column named strand was not found")
  }

  seqDat[, start:=as.integer(start)]
  seqDat[, end:=as.integer(end)]
  seqDat <- subset(seqDat, end-start+1>=minFrag & end-start+1<=maxFrag)
  if("width" %in% colnames(seqDat)) seqDat$width <- NULL

  return(seqDat)
}

#' Mapping & aggregation of modalities with genomic coordinates to a set of reference
#' coordinates.
#'
#' Convencience function for mapping different modality scores with
#' cell type and further labels such as tfs to reference coordinates. The resulting
#' table will have dimension ref coords x byCols (or ref coord x byCols`[1]` x byCols`[2]`).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#'
#' @name genomicRangesMapping
#' @param refRanges GRanges object with reference coordinates
#' @param assayTable  List of [GenomicRanges::GRanges-class], [data.table::data.table], data.frames or paths to .bed /. bam files
#' containing assay data (e.g. fragments, motif scores, peaks) to be aggregated across the reference ranges provided.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param byCols Variables across which the assays will be aggregated.
#' Will be the columns of the resulting [Matrix::Matrix-class]. If byCols is a vector with two elements,
#' the first one will constitute the list elements, the second the columns of the matrices being the list elements.
#' If they are factors (preferred) aggregation will happen across all levels, otherwise across all unique entries.
#' @param scoreCol name of the score column (e.g. motif matching scores, atac fragment counts) to be aggregated.
#' If it is NULL, the number of assay ranges overlapping the reference ranges will be counted.
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate.
#' If it is NULL, the number of assay ranges overlapping the reference ranges will be counted.
#' @param minoverlap Minimal overlap between refRanges and the assay ranges
#' Passed to [GenomicRanges::findOverlaps()]
#' @param type Type of overlap to be used, can be one of "any", "start", "end", "within", "equal".
#' Passed to [GenomicRanges::findOverlaps()]. Default is "any".
#' @param shift If Tn5 insertion bias should be considered (only if strand column is provided).
#' @param chunk If data should be processed in chunks (determined by chromosomes in `refRanges`). Recommended for large data.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return [Matrix::Matrix-class] or list of Matrices with rows corresponding to the reference ranges and columns (and list elements) to byCols.
#' @import data.table
#' @import Matrix
#' @importFrom GenomicRanges findOverlaps GRanges GRangesList
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevelsStyle keepStandardChromosomes seqlevels
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam SnowParam
#' @importFrom GenomicAlignments readGAlignmentPairs start end strand
#' @importFrom Rsamtools ScanBamParam
#' @importFrom S4Vectors split
#' @author Emanuel Sonder
#' @export
genomicRangesMapping <- function(refRanges,
                                 assayTable,
                                 byCols=c("tf_name",
                                          "cellular_context"),
                                 scoreCol=NULL,
                                 aggregationFun=NULL,
                                 minoverlap=1,
                                 type=c("any", "start", "end",
                                        "within", "equal"),
                                 shift=FALSE,
                                 chunk=NULL,
                                 BPPARAM=SerialParam()){

  type <- match.arg(type, choices=c("any","start", "end",
                                    "within", "equal"))

  # TODO: - add warning for integer overflows - data.table size
  assayTable <- .processData(assayTable, readAll=TRUE, shift=shift,
                             seqLevelStyle=seqlevelsStyle(refRanges))
  colsDepth <- unique(assayTable[[byCols[1]]])

  if(is.null(chunk)){
    chunk <- fifelse(nrow(assayTable)>1e7, TRUE, FALSE)
  }

  if(chunk){
    chrLevelsRef <- levels(seqnames(refRanges))
    assayTable <- subset(assayTable, chr %in% chrLevelsRef)
    assayTable[,chr:=factor(chr, levels=chrLevelsRef)]
    assayTable[[byCols[1]]] <- factor(assayTable[[byCols[1]]])
    if(length(byCols)>1){
      assayTable[[byCols[2]]] <- factor(assayTable[[byCols[2]]])}

    assayTable <- split(assayTable, by="chr")
    refRangesList <- S4Vectors::split(refRanges, seqnames(refRanges))

    assayTable <- assayTable[names(refRangesList)]
    overlapTable <- mapply(genomicRangesMapping,
                           refRangesList, assayTable,
                           MoreArgs=list(byCols=byCols, scoreCol=scoreCol,
                                         aggregationFun=aggregationFun,
                                         chunk=FALSE, shift=shift,
                                         type=type,
                                         BPPARAM=BPPARAM),
                           SIMPLIFY=FALSE)

    # retrieve original order
    refRangesList <- Reduce("c", refRangesList[-1], refRangesList[[1]])
    ind <- GenomicRanges::findOverlaps(refRangesList, refRanges,
                                       select="first", type="equal")

    rbindFill  <- function(mat1, mat2){

      if(is.null(mat1)) mat1 <- Matrix::Matrix(nrow=0, ncol=0)
      if(is.null(mat2)) mat2 <- Matrix::Matrix(nrow=0, ncol=0)

      allCols <- union(colnames(mat1), colnames(mat2))
      diffCols1 <- setdiff(allCols, colnames(mat1))
      diffCols2 <- setdiff(allCols, colnames(mat2))

      # get missing columns
      mat1Missing <- Matrix::Matrix(0, nrow=nrow(mat1), ncol=length(diffCols1),
                                       dimnames=list(NULL, diffCols1))
      mat1 <- cbind(mat1, mat1Missing)

      mat2Missing <- Matrix::Matrix(0, nrow=nrow(mat2), ncol=length(diffCols2),
                                       dimnames=list(NULL, diffCols2))
      mat2 <- cbind(mat2, mat2Missing)

      mat1 <- mat1[,allCols, drop=FALSE]
      mat2 <- mat2[,allCols, drop=FALSE]
      rbind(mat1, mat2)
    }

    if(length(byCols)>1){
      overlapTable <- lapply(colsDepth, function(col){
                             tablesChr <- lapply(overlapTable,
                                                 function(tables) tables[[col]])
                             tablesChr <- Reduce("rbindFill", tablesChr[-1],
                                                              tablesChr[[1]])
                             tablesChr <- tablesChr[ind,,drop=FALSE]
                             tablesChr})
      names(overlapTable) <- colsDepth
    }
    else{
      overlapTable <- Reduce("rbindFill", overlapTable[-1], overlapTable[[1]])
      overlapTable <- overlapTable[ind,,drop=FALSE]
    }

    return(overlapTable)
  }

  seqNamesCol <- "chr"
  startCol <- "start"
  endCol <- "end"

  if(sum(!(byCols %in% colnames(assayTable)))>0 | is.null(byCols)){
    stop("byCols needed to be provided and column names of assayTable.")
  }

  # attribute generic names to dimensionalities
  if(length(byCols)==2)
  {
    setnames(assayTable, byCols, c("col_depth", "col_width"))
    byCols <- c("col_depth", "col_width")
    multiTf <- TRUE
  }
  else
  {
    setnames(assayTable, byCols, c("col_width"))
    byCols <- c("col_width")
    multiTf <- FALSE
  }

  # get dimensions of tables
  nRefs <- length(refRanges)

  if(is.factor(assayTable$col_width)){
    colsWidth <- levels(assayTable$col_width)
  }
  else{
    colsWidth <- unique(assayTable$col_width)
  }
  nColsWidth <- length(colsWidth)
  # convert to integer for speed-up
  assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=colsWidth))]

  if(is.factor(assayTable$col_depth)){
    colsDepth <- levels(assayTable$col_depth)
  }
  else{
    colsDepth <- unique(assayTable$col_depth)
  }

  # convert to GRanges for faster overlap finding
  if("width" %in% colnames(assayTable)) assayTable$width <- NULL
  if("strand" %in% colnames(assayTable)) assayTable$strand <- NULL

  assayRanges <- .dtToGr(assayTable, seqCol="chr", addMetaCols=TRUE)

  # find overlaps with ref. coordinates
  overlapTable <- as.data.table(GenomicRanges::findOverlaps(refRanges,
                                                            assayRanges,
                                                            type=type,
                                                            minoverlap=minoverlap,
                                                            ignore.strand=TRUE))
  rm(refRanges, assayRanges)

  # retrieve tf and cell type ids
  overlapTable <- cbind(overlapTable$queryHits,
                        assayTable[overlapTable$subjectHits,
                                   c(byCols, scoreCol),
                                   with=FALSE])
  rm(assayTable)

  threads <- floor(getDTthreads())/BPPARAM$workers

  if(multiTf)
  {
    setkey(overlapTable, V1, col_width)
    if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- split(overlapTable, by=c("col_depth"))

    overlapTable <- BiocParallel::bplapply(overlapTable, function(table,
                                                                  scoreCol,
                                                                  aggregationFun,
                                                                  nRefs,
                                                                  colsWidth,
                                                                  threads){

      data.table::setDTthreads(threads)

      if(is.null(scoreCol) | is.null(aggregationFun)){
        table <- table[,.(value=.N),
                       by=c("V1", "col_width")]}
      else{
        table <- table[,.(value=aggregationFun(scoreCol)),
                       by=c("V1", "col_width")]}

      nColsWidth <- length(colsWidth)

      # convert to sparse matrix
      table <- sparseMatrix(i=table$V1,
                            j=as.integer(table$col_width),
                            dims=c(nRefs, nColsWidth),
                            x=table$value)
      colnames(table) <- colsWidth
      return(table)}, scoreCol=scoreCol, aggregationFun=aggregationFun,
                      nRefs=nRefs, colsWidth=colsWidth, threads=threads,
      BPPARAM=BPPARAM)
  }
  else
  {
    # setkeys for speed-up
    overlapTable[,V1:=as.integer(V1)]
    setkey(overlapTable, col_width, V1)

    # overlap with ref. coordinates
    if(is.null(scoreCol) | is.null(aggregationFun)){
      overlapTable <- overlapTable[,.(scoreCol=.N),
                                   by=c("col_width", "V1")]}
    else{
      setnames(overlapTable, scoreCol, "scoreCol")
      overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
                                   by=c("col_width", "V1")]}

    # convert to sparse matrix
    overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1,
                                         j=overlapTable$col_width,
                                         dims=c(nRefs, nColsWidth),
                                         x=overlapTable$scoreCol)
    overlapTable <- Matrix::Matrix(overlapTable)

    colnames(overlapTable) <- colsWidth
  }

  # add combinations with zero overlaps
  missingDepthCols <- setdiff(colsDepth, names(overlapTable))
  if(length(missingDepthCols)>0){
    missingMat <- Matrix(0,nrow=nRefs, ncol=nColsWidth, doDiag=FALSE)
    colnames(missingMat) <- colsWidth
    missingTables <- replicate(length(missingDepthCols),
                               missingMat)
    names(missingTables) <- missingDepthCols
    overlapTable <- c(overlapTable, missingTables)
  }

  gc()
  return(overlapTable)
}