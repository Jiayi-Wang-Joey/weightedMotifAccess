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
#'
.matchSeqlevels <- function(atacFrag, ranges) {
    frags <- data.table::rbindlist(atacFrag)
    fragSeq <- unique(as.character(frags$seqnames))
    rangeSeq <- unique(as.character(GenomicRanges::seqnames(ranges)))
    common <- intersect(fragSeq, rangeSeq)

    atacFrag <- lapply(atacFrag, function(frag) {
        frag <- frag[as.character(seqnames) %in% common, ]
        frag$seqnames <- factor(frag$seqnames)
        frag
    })

    ranges <- ranges[as.character(GenomicRanges::seqnames(ranges)) %in% common]

    list(atacFrag = atacFrag, ranges = ranges)
}

# .matchSeqlevels <- function(atacFrag, ranges) {
#     frags <- rbindlist(atacFrag)
#     fragSeq <- unique(frags$seqnames)
#     rangeSeq <- GenomicRanges::seqnames(ranges)
#     common <- intersect(fragSeq,rangeSeq)
#     atacFrag <- lapply(atacFrag, function(frag) {
#         frag <- frag[seqnames %in% common,]
#         frag$seqnames <- factor(frag$seqnames)
#         frag})
#     ranges <- ranges[seqnames(ranges) %in% common,]
#     # turn seqnames to factor
#     list(atacFrag=atacFrag, ranges=ranges)
# }

#' @title dtToGr
#'
#' @description
#' Convert a data table to GenomicRange object
#' @author Emanuel Sonder

dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end",
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


# Functionality is from another project, not yet sure how it will be packaged
#' Mapping & aggregating modalities with genomic coordinates to reference
#' coordinates.
#'
#' Internal convencience function for mapping different modality scores with
#' cell type and further labels such as tfs to reference coordinates. The resulting
#' table will have dimension ref coords x byCols (or ref coord x byCols[1] x byCols[2]).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#'
#'@name .genomicRangesMapping
#'@param refRanges GRanges object with reference coordinates
#'@param assayTable modality table to be mapped to the reference coordinates.
#'Needs to containing genomic coordinates (see args: seqNamesCol, startCol, endCol), and byCols.
#'@param byCols will be the columns /depths of the resulting matrix with dimension
#' ref coords x byCols (or ref coord x byCols[1] x byCols[2]).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#' @param seqNamesCol name of the column in motif, atac and chIP-seq data.tables containing
#' the sequence information.
#' @param startCol name of the column in motif, atac and chIP-seq data.tables containing
#' the start coordinate.
#' @param endCol name of the column in motif, atac and chIP-seq data.tables containing
#' the end coordinate.
#' @param scoreCol name of the score column (e.g. motif matching scores, atac fragment counts)
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate
#' if multiple rows of the assayTable overlap a reference coordinate.
#' @param BPPARAM BiocParallel argument either SerialParam() or MulticoreParam(workers=n)
#' @author Emanuel Sonder
#' @export
.genomicRangesMapping <- function(refRanges,
                                  assayTable,
                                  byCols=c("tf_uniprot_id",
                                           "cell_type_id"),
                                  seqNamesCol="chr",
                                  startCol="start",
                                  endCol="end",
                                  scoreCol=NULL,
                                  calledInCol=NULL,
                                  aggregationFun=NULL,
                                  minoverlap=1,
                                  BPPARAM=SerialParam()){
  # TODO: - add warning for integer overflows - data.table size
  assayTable <- copy(assayTable)

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
  nColsWidth <- length(unique(assayTable$col_width))

  # convert to integer for speed-up
  levels <- unique(assayTable$col_width)
  assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=levels))]

  # convert to GRanges for faster overlap finding
  suppressWarnings(assayTable$width <- NULL)
  suppressWarnings(assayTable$strand <- NULL)
  assayRanges <- makeGRangesFromDataFrame(as.data.frame(assayTable),
                                          keep.extra.columns=TRUE,
                                          seqnames.field=seqNamesCol,
                                          start.field=startCol,
                                          end.field=endCol,
                                          ignore.strand=TRUE)

  # find overlaps with ref. coordinates
  overlapTable <- as.data.table(findOverlaps(refRanges, assayRanges,
                                             type="any",
                                             minoverlap=minoverlap,
                                             ignore.strand=TRUE))
  rm(refRanges, assayRanges)

  # retrieve tf and cell type ids
  overlapTable <- cbind(overlapTable$queryHits,
                        assayTable[overlapTable$subjectHits,
                                   c(byCols, scoreCol, calledInCol),
                                   with=FALSE])

  if(multiTf)
  {
    setkey(overlapTable, V1, col_width)
    if(!is.null(calledInCol)) setnames(overlapTable, calledInCol, "calledInCol")
    if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- split(overlapTable, by=c("col_depth"))

    overlapTable <- BiocParallel::bplapply(overlapTable, function(table){

      if(!is.null(scoreCol) | !is.null(aggregationFun)){
        table <- table[,.(value=aggregationFun(scoreCol)),
                       by=c("V1", "col_width")]}
      else{
        table <- table[,.(value=.N),
                       by=c("V1", "col_width")]}

      # one would need to construct a second table here for the neg labels

      # convert to sparse matrix
      table <- sparseMatrix(i=table$V1,
                            j=as.integer(table$col_width), # 11.07.2024 as.integer is not needed
                            dims=c(nRefs, nColsWidth),
                            x=table$value)
      colnames(table) <- levels
      return(table)},
      BPPARAM=BPPARAM)

  }
  else
  {
    # setkeys for speed-up
    overlapTable[,V1:=as.integer(V1)]
    setkey(overlapTable, col_width, V1)

    if(!is.null(scoreCol) | !is.null(aggregationFun)){
      setnames(overlapTable, scoreCol, "scoreCol")
      overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
                                   by=c("col_width", "V1")]}
    else{
      overlapTable <- overlapTable[,.(scoreCol=.N),
                                   by=c("col_width", "V1")]}

    # convert to sparse matrix
    overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1,
                                         j=overlapTable$col_width,
                                         dims=c(nRefs, nColsWidth),
                                         x=overlapTable$scoreCol)

    colnames(overlapTable) <- levels
  }

  gc()
  return(overlapTable)
}

#' @author Pierre-Luc
getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
                                  species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
                     TF=gsub(pat,"",names(motifs)),
                     grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
         universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
         PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
           universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
           modf$TF)),
         PWMatrix=do.call(TFBSTools::PWMatrixList,
                          setNames(universalmotif::convert_motifs(motifs,
                                                                  class="TFBSTools-PWMatrix"), modf$TF))
  )
}

#' Get motif matches
#'
#' @param genome A BSgenome respective to the genome used for the alignment.
#' @param peakpath The path to the merged ATAC-seq peaks from control and treatment conditions.
#' @param spec Species. Either "Hsapiens" or "Mmusculus".
#' @param seqStyle Either "ensembl" or "UCSC" depending on the format of the peak file.
#' @param srcFolder Folder where to find scripts and important objects
#'
#' @author Pierre-Luc
#'

getpmoi <- function(genome,
                    peaks,
                    spec=c("Homo sapiens", "Mus musculus"), minHits=50L,
                    seqStyle=c("ensembl", "NCBI","UCSC"), keepTop=NULL,
                    srcFolder, motifs=NULL, thresh=NULL){
  seqStyle <- match.arg(seqStyle)

  if(is.character(peaks)){
    peals <- sort(rtracklayer::import(peaks))
  }
  peaks <- keepStandardChromosomes(peaks,
                                   pruning.mode = "coarse")

  seqlevelsStyle(genome) <- seqStyle
  peak_seqs <- memes::get_sequence(peaks, genome)

  # Get the motifs in universal format required by memes
  if(is.null(motifs)){
    if(spec=="Homo sapiens"){
      spec <- "Hsapiens"}
    else if(spec=="Mus musculus"){
      spec <- "Mmusculus"}

    motifs <- getNonRedundantMotifs("universal", species = spec)
  }

  # Obtain the positions of motif instances which are later required as input for runATAC
  if(is.null(thresh)) thresh <- 1e-4
  pmoi <- memes::runFimo(peak_seqs,
                  motifs, thresh=thresh,
                  meme_path="/common/meme/bin/",
                  skip_matched_sequence=TRUE)
  pmoi <- pmoi[order(mcols(pmoi)$motif_id, -mcols(pmoi)$score)]
  pmoi <- split(pmoi, pmoi$motif_id)
  if(!is.null(minHits)) pmoi <- pmoi[which(lengths(pmoi)>=minHits)]
  if(!is.null(keepTop)) pmoi <- lapply(split(pmoi, pmoi$motif_id), n=keepTop, FUN=head)
  pmoi <- sort(unlist(GRangesList(pmoi)))
  return(pmoi)
}

#' @author Emanuel Sonder
.processData <- function(data, readAll=FALSE, shift=FALSE,
                         subSample=NULL, seqLevelStyle="UCSC"){
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
  if("width" %in% colnames(seqDat)) seqDat$width <- NULL

  return(seqDat)
}
