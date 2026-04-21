#' computeDeviationsWeighted
#'
#' @param counts A matrix of aggregated (weighted) counts, with motifs as rows
#'   and samples as columns.
#' @param annotations A matrix of motif annotations, with peaks as rows and
#'   motifs as columns.
#' @param bgcounts A matrix of background-weighted counts, with peaks as rows
#'   and samples as columns.
#' @param bg An optional `bcvBackground` object. Either `bg` or `bgcounts` must
#'   be given.
#'
#' @returns A list with two matrices, for deviations and z-scores.
#' @author Pierre-Luc Germain
#' @importFrom Matrix crossprod t sparseMatrix
#' @importFrom betterChromVAR computeBackgrounds
#' @export
computeDeviationsWeighted <- function(counts, annotations, bgcounts=NULL, bg=NULL){
  stopifnot(!is.null(bg) || !is.null(bgcounts))
  if(is.null(bg)){
    bg <- computeBackgrounds(bgcounts, getBackgroundBins(bgcounts))
  }else if(length(bg@depth)==0 || length(bg@expectation)==0){
    stop("Incomplete background; please run computeBackgrounds() first.")
  }else if(length(bg@depth) != ncol(counts)){
    stop("The background provided does not seem to fit `counts`")
  }else if(length(bg@expectation) != nrow(annotations)){
    stop("The background provided does not seem to fit `annotations`")
  }
  stopifnot(nrow(counts)==ncol(annotations))

  binMap <- bg@peak2bin
  bin2peakMat <- sparseMatrix(i=binMap, j=seq_along(binMap),
                              dims=c(nrow(bg@binBinProbs), length(binMap)))
  motifBinCounts <- Matrix::t(annotations) %*% Matrix::t(bin2peakMat)

  motif_bg_exp <- as.matrix(motifBinCounts %*% bg@E)

  numerator <- counts - motif_bg_exp

  # z-scores:
  z <- numerator / sqrt(pmax(0, as.matrix(motifBinCounts %*% bg@V)))

  # deviations:
  globalMotifAvg <- as.vector(Matrix::crossprod(annotations, bg@expectation))
  sf <- bg@depth / sum(bg@expectation)
  deviations <- numerator/outer(globalMotifAvg, sf)

  list(deviations=deviations, z=z)
}

.getInsertsPos <- function(atacFrag, motifData, stranded, shiftLeft){

  if(stranded){
    specColsMotif <- c("motif_center", "start",
                       "end", "motif_id", "motif_match_id", "strand")}
  else{
    specColsMotif <- c("motif_center", "start",
                       "end", "motif_id", "motif_match_id")}

  if(stranded) specColsFrag <- c("sample", "strand") else specColsFrag <- c("sample")

  # convert to granges for faster overlapping
  motifMarginRanges <- .dtToGr(motifData, startCol="start_margin",
                              endCol="end_margin", seqCol="chr",
                              stranded=stranded)
  atacStartRanges <- .dtToGr(atacFrag, startCol="start", endCol="start",
                            seqCol="chr", stranded=stranded)
  atacEndRanges <- .dtToGr(atacFrag, startCol="end", endCol="end",
                          seqCol="chr", stranded=stranded)

  startHits <- GenomicRanges::findOverlaps(atacStartRanges, motifMarginRanges,
                                           type="within", ignore.strand=TRUE)
  endHits <- GenomicRanges::findOverlaps(atacEndRanges, motifMarginRanges,
                                         type="within", ignore.strand=TRUE)

  # get overlapping insertion sites
  atacStartInserts <- atacFrag[queryHits(startHits),
                               c(specColsFrag, "start"), with=FALSE]
  atacEndInserts <- atacFrag[queryHits(endHits),
                             c(specColsFrag, "end"), with=FALSE]
  setnames(atacStartInserts, "start", "insert")
  if(stranded) setnames(atacStartInserts, "strand", "strand_insert")
  setnames(atacEndInserts, "end", "insert")
  if(stranded) setnames(atacEndInserts, "strand", "strand_insert")

  ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
              motifData[c(subjectHits(startHits),
                          subjectHits(endHits)), specColsMotif, with=FALSE])

  # count insertions around motif
  ai[,rel_pos:=insert-motif_center]
  ai[,type:=fifelse(insert>=start & insert<=end, 1,0)]
  ai[,ml:=end-start+1L, by=motif_id] # motif length

  if(stranded){
    # take strandedness of fragment into account
    if(nrow(ai)>0){
      if(data.table::first(ai$ml) %% 2!=0){
        ai[,rel_pos:=fifelse(strand_insert=="-", -1L*rel_pos, rel_pos)]
      }
      else{
        aiMotif <- subset(ai, type==1)
        if(shiftLeft){
          ai[,rel_pos:=fifelse(strand_insert=="-", -1L*rel_pos-1L, rel_pos)]
          ai[,rel_pos:=fifelse(strand=="-", -1L*rel_pos-1L, rel_pos)]
        }
        else{
          ai[,rel_pos:=fifelse(strand_insert=="-", -1L*rel_pos+1L, rel_pos)]
          ai[,rel_pos:=fifelse(strand=="-", -1L*rel_pos+1L, rel_pos)]
        }
      }}
  }

  return(ai)
}

#' Tn5 insertion counting
#'
#' Counts Tn5 insertions around provided motif-matches.
#' If requested also computes insertion footprint profiles and weighted insertion counts.
#'
#' @name getInsertionProfiles
#' @param atacData Named list of [GenomicRanges::GRanges-class], [data.table::data.table], data.frames or paths to .bed /. bam files
#' containing ATAC-seq fragment coordinates (i.e. chr/seqnames, start, end and optionally a strand column). List names will be used as sample names.
#' If a single object is provided and it contains a column named "sample", insertion counts will be computed for each sample separately.
#' @param motifRanges [GenomicRanges::GRanges-class] object containing coordinates of motif-matches.
#' Needs to contain a metadata column `motif_id` if insertion profiles should be computed for several motifs separately.
#' @param margin Margin around motif-matches to consider for computing Tn5 insertion events
#' @param shift If Tn5 insertion bias should be considered (only if strand column is provided).
#' @param calcProfile If insertion footprint profiles should be computed.
#' @param profiles Pre-computed insertion footprint profile as obtained by this function when running with `calcProfile=TRUE`.
#' Needs to contain a column with relative positions wrt to the center of the motif match (termed "rel_pos") and weight column (termed "w").
#' @param symmetric If transcription factor footprint profiles should be symmetric around the motif matches. Only used if `calcProfile=TRUE`.
#' @param stranded If insertion footprint profiles should be computed taking strandedness of fragments into account.
#' @param subSample If fragments should be sub-sampled for speed-up.
#' Default is no sub-sampling, if a number is provided the fragments of each file/object provided will be subsampled to that number.
#' @param simplified If return should be simplified and be provided as a [SummarizedExperiment::RangedSummarizedExperiment-class] object.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()].
#' @return If `simplified=TRUE` a [data.table::data.table] containing insertion counts within and in margins around motif matches and weighted insertion counts in case
#' an insertion profile is provided or if `calcProfile=TRUE`.
#' If `calcProfile=TRUE` also a footprint profile around the motif matches is returned, containing a weight ("w") column corresponding to relative insertion
#' frequency at the respective position relative to the motif center ("rel_pos").
#' If `simplified=FALSE` a [SummarizedExperiment::RangedSummarizedExperiment-class] object is returned with ranges corresponding to the `motifRanges` provided as an input.
#' Assays correspond to the different insertion counts computed, columns will correspond to the samples. If `calcProfile=TRUE` the profiles will be saved in the metadata.
#' @import data.table
#' @importFrom GenomicRanges findOverlaps GPos resize GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @author Emanuel Sonder
#' @export
getInsertionProfiles <- function(atacData,
                                 motifRanges,
                                 margin=200,
                                 shift=FALSE,
                                 calcProfile=TRUE,
                                 profiles=NULL,
                                 symmetric=FALSE,
                                 stranded=FALSE,
                                 subSample=NULL,
                                 simplified=FALSE,
                                 BPPARAM=SerialParam()){

  # prep motif data
  motifData <- .processData(motifRanges, shift=FALSE, readAll=FALSE)
  if(!("motif_id" %in% colnames(motifData))){
    message("Assuming all ranges are of the same type")
    motifData[,motif_id:=1L]
  }

  if(!calcProfile & !is.null(profiles) &
     length(setdiff(motifData$motif_id, names(profiles)))>0){
      warning("Not all motif-ranges have an insertion-profile provided.
      If wished to use a pre-computed profile provide one for all the motifs specified in the motifRanges arg.
      Switching to computing profiles for all (calcProfile=TRUE).")
    calcProfile <- TRUE
  }

  margin <- as.integer(margin)
  if(margin>0){
    motifMarginRanges <- as.data.table(GenomicRanges::resize(motifRanges,
                                                             width=2L*margin,
                                                             fix="center"))
  }
  else{
    motifMarginRanges <- as.data.table(motifRanges)
  }
  setnames(motifMarginRanges, c("start", "end"), c("start_margin", "end_margin"))
  motifData <- cbind(motifData, motifMarginRanges[,c("start_margin", "end_margin"), with=FALSE])

  # prep ATAC fragment data
  if(is.data.table(atacData)) atacData <- list(atacData)
  atacFrag <- lapply(atacData, .processData, shift=shift, subSample=subSample)
  if(!("sample" %in% colnames(atacFrag[[1]]))){
    names(atacFrag) <- names(atacData)
    atacFrag <- rbindlist(atacFrag, idcol="sample")
  }
  else{
    atacFrag <- rbindlist(atacFrag)
  }

  commonChr <- intersect(unique(motifData$chr),unique(atacFrag$chr))
  chrLevels <- commonChr
  if(length(chrLevels)==0){
    stop("No common chromosomes found in atacData and motifRanges")
  }

  atacFrag <- subset(atacFrag, chr %in% chrLevels)
  motifData <- subset(motifData, chr %in% chrLevels)

  motifLevels <- unique(motifData$motif_id)

  # convert to factors (memory usage)
  motifData[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
  motifData[,motif_id:=factor(motif_id, levels=motifLevels, ordered=TRUE)]

  # determine motif center
  motifData[,motif_center:=as.integer(floor((end-start)/2))+start]

  motifData[,end_margin:=fifelse(end_margin-motif_center<margin,end_margin+1L,end_margin)]
  motifData[,start_margin:=fifelse(start_margin-motif_center>-margin,start_margin-1L,start_margin)]

  distEnd <- motifData$end[1] - motifData$motif_center[1]
  distStart <- motifData$motif_center[1] - motifData$start[1]
  if(distStart>distEnd){
    shiftLeft <- TRUE}
  else{
    shiftLeft <- FALSE
  }

  atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
  medZero <- function(x, len){ median(c(rep(0L,max(0,len-length(x))),x)) }
  nSamples <- length(unique(atacFrag$sample))

  setorder(motifData, chr)
  setorder(atacFrag, chr)
  motifData[,motif_match_id:=1:nrow(motifData)]

  motifData <- split(motifData, by="chr")
  atacFrag <- split(atacFrag, by="chr")

  if(calcProfile){
   message("Computing insertion-profiles")
   atacProfiles <- BiocParallel::bpmapply(function(md, af, stranded, shiftLeft){
      atacInserts <- .getInsertsPos(af, md, stranded, shiftLeft)
      atacProfile <- atacInserts[,.(pos_count_global=.N),
                                 by=.(ml, rel_pos, motif_id, type)]
      return(atacProfile)},
    motifData, atacFrag, MoreArgs=list(stranded=stranded, shiftLeft=shiftLeft),
    SIMPLIFY=FALSE, BPPARAM=BPPARAM)
    atacProfiles <- rbindlist(atacProfiles, idcol="seqnames")

    atacProfiles <- atacProfiles[,.(pos_count_global=sum(pos_count_global)),
                                 by=.(rel_pos, motif_id, type, ml)]
    # get inserts within motif
    atacProfilesMotif <- subset(atacProfiles, type==1)
    atacProfilesMotif[,med_pos_count_global:=medZero(pos_count_global,
                                                     data.table::first(ml)),
                      by=.(motif_id)]
    atacProfilesMotif <- atacProfilesMotif[,
   .(pos_count_global=(pos_count_global+data.table::first(med_pos_count_global))/2),
                                           by=.(rel_pos, motif_id, type)]

    atacProfilesMargin <- subset(atacProfiles, type==0)

    # fill non covered positions
    atacProfiles <- rbind(atacProfilesMargin, atacProfilesMotif, fill=TRUE)
    allPos <- data.table(expand.grid(motifLevels, seq(-margin,margin)))

    colnames(allPos) <- c("motif_id", "rel_pos")
    allPos$motif_id <- factor(allPos$motif_id, levels=motifLevels, ordered=TRUE)
    allPos$pos_count_global <- 0

    colsProfile <- c("rel_pos", "motif_id", "pos_count_global")
    atacProfiles[,rel_pos:=as.integer(rel_pos)]

    atacProfiles <- rbind(atacProfiles[, colsProfile, with=FALSE],
                          allPos[!atacProfiles, on=c("rel_pos", "motif_id")])

    # calculate weights
    setorder(atacProfiles, motif_id, rel_pos)
    if(max(atacProfiles$pos_count_global)>0){
      atacProfiles[,w:=pos_count_global/max(pos_count_global), by=motif_id]
      atacProfiles[,w_smooth:=smooth(w), by=motif_id]}
    else{
      atacProfiles[,w:=1]
    }

    if(symmetric){
      atacProfiles[,w:=rev(w)+w, by=motif_id]
      if("w_smooth" %in% colnames(atacProfiles)){
        atacProfiles[,w_smooth:=rev(w_smooth)+w_smooth, by=motif_id]}
    }

    atacProfiles[,w:=w/sum(w), by=motif_id]
    if("w_smooth" %in% colnames(atacProfiles)){
      atacProfiles[,w_smooth:=w_smooth/sum(w_smooth), by=motif_id]}
  }
  else{
    atacProfiles <- profiles
    if(is.list(atacProfiles)){
      message("Skipped insertion-profiles computation. Using provided pre-computed ones")
      atacProfiles <- rbindlist(atacProfiles, idcol="motif_id")}
  }

  # get match scores
  motifScores <- BiocParallel::bpmapply(function(md,af,
                                                 stranded,
                                                 profiles,
                                                 shiftLeft){

    atacInserts <- .getInsertsPos(af, md, stranded, shiftLeft)

    if(!is.null(profiles)){
      atacInserts <- atacInserts[,.(pos_count=.N),
                                 by=.(motif_match_id, motif_id, sample,
                                      rel_pos, type)]
      atacInserts <- merge(atacInserts,
                           profiles[,c("rel_pos", "motif_id", "w"),with=FALSE],
                           by.x=c("motif_id","rel_pos"),
                           by.y=c("motif_id","rel_pos"), all.x=TRUE, all.y=FALSE)
      atacInserts[,score:=w*pos_count]
      atacInserts[,dev:=(pos_count/sum(pos_count)-w)^2/(w),
                   by=.(motif_match_id, motif_id, sample)]
      atacInsertSum <- atacInserts[,.(score=sum(score),
                                      chi2=sum(dev)+(1-sum(w)), # add deviation for positions with 0 inserts
                                      tot_count=sum(pos_count)),
                                   by=.(motif_match_id, motif_id, sample, type)]
    }
    else{
      atacInsertSum <- atacInserts[,.(tot_count=.N),
                                   by=.(motif_match_id, motif_id, sample, type)]
    }
    return(atacInsertSum)
  }, motifData, atacFrag,
  MoreArgs=list(stranded=stranded,
                profiles=atacProfiles,
                shiftLeft=shiftLeft),
  SIMPLIFY=FALSE,
  BPPARAM=BPPARAM)

  motifScores <- rbindlist(motifScores)

  motifData <- rbindlist(motifData)
  motifData[,chr:=chrLevels[chr]]
  motifScores <- cbind(motifScores,
                       motifData[motifScores$motif_match_id,
                                 c("start", "end", "chr"), with=FALSE])
  motifScores[,type:=factor(fifelse(type=="0", "margin", "within"),
                            levels=c("margin", "within"))]

  if("tot_count" %in% colnames(motifScores)){
  setnames(motifScores, c("tot_count"), INSERTFEATNAME)}
  if("score" %in% colnames(motifScores)){
    setnames(motifScores, c("score", "chi2"), c(WINSERTSFEATNAME, DEVFEATNAME))}

  if(simplified){
    scoreCols <- intersect(c(INSERTFEATNAME, WINSERTSFEATNAME, DEVFEATNAME),
                           colnames(motifScores))
    assayMats <- lapply(scoreCols, function(SCORECOL){
      ms <- motifScores[,.(score=sum(get(SCORECOL))), by=.(motif_id, sample,
                                                           motif_match_id,
                                                           chr, start, end)]
      am <- genomicRangesMapping(motifRanges, assayTable=ms,
                                 byCols="sample",
                                 scoreCol="score",
                                 aggregationFun=max,
                                 type="equal", #otw this does not work in case motif matches of different TFs are overlapping
                                 BPPARAM=BPPARAM)})

    names(assayMats) <- scoreCols
    res <- SummarizedExperiment(rowRanges=motifRanges,
                                assays=assayMats,
                                metadata=list(profiles=atacProfiles))
  }
  else{
    res <- list(motifScores, atacProfiles)
    names(res) <- c(RETSCORESNAME, REPROFILENAME)
  }

  return(res)
}


.getInsertionBgProfiles <- function(atacFrag,
    coords,
    #minWidth=30,
    #maxWidth=2000,
    symmetric=TRUE,
    libNorm=FALSE,
    chunk=TRUE){

    # prep motif data
    coordData <- as.data.table(coords)
    setnames(coordData, "seqnames", "chr")
    chrLevels <- unique(coordData$chr)
    motifLevels <- unique(coordData$id)

    # convert to factors (memory usage)

    coordData[,chr:=as.integer(factor(chr,
        levels=chrLevels, ordered=TRUE))]
    coordData[,id:=as.integer(factor(id, levels=motifLevels, ordered=TRUE))]

    # determine motif center
    coordData[,center:=floor((end-start)/2)+start]

    # convert to factors (memory usage)
    if("seqnames" %in% colnames(atacFrag)){
        setnames(atacFrag, "seqnames", "chr")
    }
    #atacFrag <- copy(atacFrag) #TODO: take out that copy
    atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]

    # if("seqnames" %in% colnames(atacFrag)){
    #   setnames(atacFrag, "seqnames", "chr")
    # }

    nSamples <- length(unique(atacFrag$sample))

    setorder(coordData, chr)
    setorder(atacFrag, chr)
    coordData[,coord_id:=1:nrow(coordData)]
    coordData <- split(coordData, by="chr")
    atacFrag <- split(atacFrag, by="chr")

    atacInserts <- mapply(function(md,af){

        #md[,motif_match_id:=1:nrow(md)]

        # convert to granges for faster overlapping
        coordRanges <- dtToGr(md, startCol="start", endCol="end", seqCol="chr")
        atacStartRanges <- dtToGr(af, startCol="start", endCol="start", seqCol="chr")
        atacEndRanges <- dtToGr(af, startCol="end", endCol="end", seqCol="chr")

        startHits <- findOverlaps(atacStartRanges,
            coordRanges, type="within") # check if type within faster or slower
        endHits <- findOverlaps(atacEndRanges, coordRanges, type="within")

        # get overlapping insertion sites
        atacStartInserts <- af[queryHits(startHits), c("sample", "start"), with=FALSE]
        atacEndInserts <-af[queryHits(endHits), c("sample", "end"), with=FALSE]
        setnames(atacStartInserts, "start", "insert")
        setnames(atacEndInserts, "end", "insert")

        ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
            rbindlist(list(
                md[subjectHits(startHits), c("center", "start",
                    "end", "id", "coord_id")],
                md[subjectHits(endHits), c("center", "start", "end", "id",
                    "coord_id")])))

        # count insertions around motif
        ai[,rel_pos:=abs(insert-center)]

        ai <- ai[,.(pos_count_sample=.N),  by=.(coord_id, rel_pos, sample, id)]
        ai
    },
        coordData,
        atacFrag,
        SIMPLIFY=FALSE)

    # combine insertion counts across chromosomes
    atacInserts <- rbindlist(atacInserts)
    atacProfiles <- atacInserts[,.(pos_count_sample=sum(pos_count_sample)), by=.(id, rel_pos, sample)]
    atacProfiles <- atacProfiles[,pos_count_global:=sum(pos_count_sample), by=.(id, rel_pos)]
    setorder(atacProfiles, id, rel_pos)
   # atacProfiles[,w:=smooth(pos_count_global/sum(pos_count_global),
  #      twiceit=TRUE), by=id]
    atacProfiles[,w:=smooth(pos_count_global, twiceit=TRUE), by=id]
    atacProfiles[,w:=w*length(w)/sum(w), by=id]

    atacProfiles <- atacProfiles[,.(w=first(w)), by=.(rel_pos, id)]
    if(symmetric) atacProfiles[,w:=rev(w)+w, by=id]
    atacProfiles[,w:=w/sum(w)]
    atacInserts <- merge(atacInserts,
        atacProfiles[,c("rel_pos", "id", "w"), with=FALSE],
        by.x=c("id", "rel_pos"),
        by.y=c("id", "rel_pos"))
    if (libNorm) {
        atacInserts[,pos_count_sample_norm:=
                (pos_count_sample/sum(pos_count_sample)), by=sample]
        atacInserts[,score:=w*pos_count_sample_norm]
    } else {
        atacInserts[,score:=w*pos_count_sample]
    }


    # calculate per sample motif match scores
    scores <- atacInserts[,.(score=sum(score)), by=.(coord_id, sample)]
    #scores[,score:=sum(score),by=id]

    # get back original coordinates
    coordData <- rbindlist(coordData)
    scores <- cbind(coordData[scores$coord_id,setdiff(colnames(coordData), "score"), with=FALSE],
        scores)
    scores[,chr:=chrLevels[seqnames]]
    #scores[,motif_name:=motifLevels[motif_id]]

    return(scores)
}



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
    isReplicated=FALSE,
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
    assayRanges <- makeGRangesFromDataFrame(as.data.frame(assayTable),
        keep.extra.columns=TRUE,
        seqnames.field=seqNamesCol,
        start.field=startCol,
        end.field=endCol,
        ignore.strand=TRUE)

    # find overlaps with ref. coordinates
    overlapTable <- as.data.table(findOverlaps(assayRanges, refRanges,
        type="within",
        minoverlap=minoverlap,
        ignore.strand=TRUE))
    rm(refRanges, assayRanges)


    # retrieve tf and cell type ids
    overlapTable <- cbind(overlapTable$subjectHits,
        assayTable[overlapTable$queryHits,
            c(byCols, scoreCol, calledInCol), with=FALSE])

    if(multiTf)
    {
        setkey(overlapTable, V1, col_width)
        if(!is.null(calledInCol)) setnames(overlapTable, calledInCol, "calledInCol")
        if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
        overlapTable <- split(overlapTable, by=c("col_depth"))

        #TODO:  BPPARAM= serialParam or multiCoreParam() as function arguments
        #TODO:  bplapply(fls[1:3], FUN, BPPARAM = MulticoreParam(), param = param)
        # overlap with ref. coordinates
        overlapTable <- BiocParallel::bplapply(overlapTable, function(table){

            if(is.null(aggregationFun))
            {
                # get number of max callers
                #table[, n_max:=length(unique(calledInCol)),
                #        by=c("col_width")]

                if(isReplicated)
                {
                    table[,rep:=tstrsplit(calledInCol, split="-", keep=2)]
                    table[, n_max:=length(unique(calledInCol)),
                        by=c("col_width", "rep")]

                    table <- table[,.(value=.chIPlabelRule(calledInCol,
                        data.table::first(n_max))),
                        by=c("V1", "col_width", "rep")]
                    table <- table[,.(value=.chIPReplicateLabelRule(value)),
                        by=c("V1", "col_width")]
                }
                else{
                    table[,n_max:=length(unique(calledInCol)),
                        by=c("col_width")]

                    table <- table[,.(value=.chIPlabelRule(calledInCol,
                        data.table::first(n_max))),
                        by=c("V1", "col_width")]
                }
            }
            else
            {
                table <- table[,.(value=aggregationFun(scoreCol)),
                    by=c("V1", "col_width")]
            }

            # one would need to construct a second table here for the neg labels

            # convert to sparse matrix
            table <- sparseMatrix(i=table$V1,
                j=as.integer(table$col_width), # what if non-numeric cell-id
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

        # overlap with ref. coordinates
        if(is.null(aggregationFun)) error("Aggregation function needs to be defined")

        # why an lapply 05.04.24
        #overlapTable <- overlapTable[,.("scoreCol"=lapply(.SD, aggregationFun, na.rm=TRUE)),
        #                              by=c("col_width", "V1"),
        #                             .SDcols=scoreCol]
        setnames(overlapTable, scoreCol, "scoreCol")
        overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
            by=c("col_width", "V1")]

        # convert to sparse matrix
        overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1,
            j=overlapTable$col_width,
            dims=c(nRefs, nColsWidth),
            x=overlapTable$scoreCol)

        colnames(overlapTable) <- levels
    }

    return(overlapTable)
}




#' Function to get, from the counts-per-match and for a motif of choice,
#' a matrix of the counts from the best motif match per peak (0 if none)
#' @param matchScore a data table containing motif ranges, motif_id, sample
#' and match score (output by `.getInsertProfile`)
#' @param peakRange a GRange object containing the peak ranges
#' @return a list of counts-per-peaks for motif X


.getPeakMatchScore <- function (matchScore,
    peakRange,
    seqCol="seqnames",
    ...
    ) {
    if (!is.null(peakRange$score)) peakRange$score <- NULL

    peaks <- as.data.table(peakRange)
    peaks$peakID <- seq_len(nrow(peaks))
    if ("strand" %in% colnames(matchScore)) matchScore$strand <- NULL
    res <- .genomicRangesMapping(peakRange, matchScore,
        byCols=c("motif_id", "sample"), seqNamesCol = seqCol,
        scoreCol="score", aggregationFun = max)
    lapply(res, \(x) as.matrix(x))

}

#' Function to get the background peak sets and strongest motif match scores in
#' each peak
#' @param se a SummarizedExperiment object that containing original fragment count
#' on peak level
#' @param genome BSgenome object
#' @param niterations number of background peaks to sample
#' @return a list containing background peak indices and a matrix containing the
#' strongest motif match score for peaks per sample

.getBackgroundPeaks <- function(se, genome, niterations) {
    se <- addGCBias(se,genome=genome)
    bg <- getBackgroundPeaks(object=se, niterations=niterations)
    #bgProfile <- .getInsertionBgProfiles(fragDt, peakRange)
}

#' get the motif activity score
#' @param backgroundPeaks a data.table containing background peak indices
#' @param peakMatchScore a list of matrix containing match score on peak level
#' for each motif
#' @param maxScorePeaks a matrix containing strongest motif match score for
#' peaks per sample
#' @return a list of matrix containing motif activity score on peak level per sample

.getMotifActivityScore <- function (backgroundPeaks,
    peakMatchScore,
    bgScore) {

    # compute the deviations for each motif
    activityScore <- lapply(peakMatchScore, function(score) {
        # find which peaks containing that motif
        idx <- which(rowSums(score)!=0)
        # compute background deviations
        bg_scores <- sapply(seq_len(ncol(backgroundPeaks)), function(x) {
            bg_i <- backgroundPeaks[idx,x]
            peak_set <- bgScore[bg_i,]
            # deviations
            bg_score <- colSums(peak_set)
        })

        motif_score <- colSums(score[idx,])
        bg_score <- rowMeans(bg_scores)
        bg_motif_score <- motif_score - rowMeans(bg_scores)
        z <- bg_motif_score/rowSds(bg_scores)

        dev_motif <- (motif_score-mean(motif_score))/mean(motif_score)
        dev_bg <- (bg_scores-colMeans(bg_scores))/colMeans(bg_scores)
        dev_z <- (dev_motif - rowMeans(dev_bg))/rowSds(dev_bg)
        # list(motif_score=motif_score, bg_motif_score=bg_motif_score, z=z,
        #     dev_motif=dev_motif, dev_bg=rowMeans(dev_bg), dev_z=dev_z,
        #     bg_score=bg_score)
        list(motif_score=motif_score, bg_scores=bg_scores)

    })
    return(activityScore)
}


computeMotifActivityScore <- function (se,
    atacFrag,
    motifRanges,
    species,
    minFrag = 30,
    maxFrag = 3000,
    genome,
    niterations=100,
    nullModel=FALSE,
    symmetric=TRUE,
    libNorm=FALSE,
    ...
    ) {

    # already resized peaks
    se <- filterPeaks(se, non_overlapping = TRUE)
    peakRange <- rowRanges(se)
    #seqlevelsStyle(peakRange) <- "NCBI"
    motifRanges <- .standardChromosomes(motifRanges, species = species)
    atacFrag <- lapply(atacFrag, function(dt) {
        gr <- dtToGr(dt)
        gr <- .standardChromosomes(gr, species = species)
        as.data.table(gr)
    })
    atacFrag <- lapply(atacFrag, \(x) {x$strand <- NULL; x$count <- NULL; x})
    # filter too short or too long fragments
    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)

    # match seqLevels
    res <- .matchSeqlevels(atacFrag, peakRange)
    atacFrag <- res$atacFrag
    peakRange <- res$ranges

    # calculate activity score
    fragDt <- rbindlist(atacFrag, idcol="sample")
    fragDt[,sample:=factor(sample)]
    rm(atacFrag)
    gc()
    frg <- copy(fragDt)
    matchScore <- .getInsertionProfiles(fragDt,
        motifRanges=motifRanges,
        nullModel=nullModel,
        symmetric=symmetric,
        libNorm=libNorm)
    fragDt <- frg
    rm(frg)
    peakMatchScore <- .getPeakMatchScore(matchScore, peakRange)
    backgroundPeaks <- .getBackgroundPeaks(se, genome=genome,
        niterations=niterations)
    mcols(peakRange)$id <- 1
    bgProfile <- .getInsertionBgProfiles(fragDt, peakRange,
        libNorm=libNorm, symmetric=symmetric)
    bgScore <- as.matrix(sparseMatrix(i=bgProfile$coord_id,
        j=as.integer(factor(bgProfile$sample)), x=as.numeric(bgProfile$score),
        dims=dim(se), dimnames=list(rownames(se),colnames(se))))


    # deviations, z-score, colSums of motifs and (mean deviations of) backgrounds
    activityScore <- .getMotifActivityScore(
        backgroundPeaks,
        peakMatchScore,
        bgScore)

    .m <- \(name) {
        res <- t(sapply(activityScore, \(x) x[[name]]))
        colnames(res) <- colnames(se)
        res
    }
    dt <- matchScore[,.(score=sum(score)), by=.(motif_id,sample)]
    ms <- reshape2::acast(dt, motif_id ~ sample, value.var = "score")

    # asy <- list(max_motif_score=.m("motif_score"), bg_score=.m("bg_score"),
    #     bg_motif_score=.m("bg_motif_score"), z=.m("z"),
    #     dev_motif=.m("dev_motif"), dev_z=.m("dev_z"), dev_bg=.m("dev_bg"),
    #     sum_motif_score=ms)
    # se <- SummarizedExperiment(assays = asy)

    return(list(activityScore=activityScore, sum_motif_score=ms))

}


