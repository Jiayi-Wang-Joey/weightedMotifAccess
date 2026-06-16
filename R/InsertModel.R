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
#' @param minFrag Minimum fragment length to be considered for counting insertions. Default is 30bp.
#' @param maxFrag Maximum fragment length to be considered for counting insertions. Default is
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
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importClassesFrom GenomicRanges GRanges
#' @author Emanuel Sonder
#' @export
getInsertionProfiles <- function(atacData,
                                 motifRanges,
                                 margin=200,
                                 shift=FALSE,
                                 minFrag=30, 
                                 maxFrag=3000,
                                 calcProfile=TRUE,
                                 profiles=NULL,
                                 symmetric=FALSE,
                                 stranded=FALSE,
                                 subSample=NULL,
                                 simplified=FALSE,
                                 BPPARAM=SerialParam()){

  # prep motif data
  seqLevelStyle <- GenomeInfoDb::seqlevelsStyle(motifRanges)
  motifData <- .processData(motifRanges, shift=FALSE, readAll=FALSE, 
                            seqLevelStyle=seqLevelStyle)

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
  atacFrag <- lapply(atacData, .processData, shift=shift, subSample=subSample, 
                     seqLevelStyle=seqLevelStyle, minFrag=minFrag, maxFrag=maxFrag, 
                     filterLength=TRUE)
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
  medZero <- function(x, len){floor(median(c(rep(0L,max(0,len-length(x))),x))) }
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
      atacProfiles <- rbindlist(atacProfiles, idcol="motif_id", use.names=TRUE, fill=TRUE, 
                                ignore.attr=TRUE)}
  }
  atacProfiles[,motif_id:=as.character(motif_id)]

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
      if("w_smooth" %in% colnames(profiles)) profiles[,w:=as.numeric(w_smooth)]
      atacInserts <- merge(atacInserts,
                           profiles[,c("rel_pos", "motif_id", "w"),with=FALSE],
                           by.x=c("motif_id","rel_pos"),
                           by.y=c("motif_id","rel_pos"), all.x=TRUE, all.y=FALSE, 
                           allow.cartesian=TRUE)
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
  setnames(motifScores, c("tot_count"), "inserts")}
  if("score" %in% colnames(motifScores)){
    setnames(motifScores, c("score", "chi2"), c("w_inserts", "dev"))}

  if(simplified){
    scoreCols <- intersect(c("inserts", "w_inserts", "dev"),
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
    names(res) <- c("insertion_weighted_counts", "insertion_profile")
  }

  return(res)
}

#' @name getInsertionProfiles
# todo add BPPARAM argument
#' @param atacData Named list of [GenomicRanges::GRanges-class], [data.table::data.table], data.frames or paths to .bed /. bam files
#' containing ATAC-seq fragment coordinates (i.e. chr/seqnames, start, end and optionally a strand column). List names will be used as sample names.
#' If a single object is provided and it contains a column named "sample", insertion counts will be computed for each sample separately.
#' @param peakRanges [GenomicRanges::GRanges-class] object containing coordinates of ATAC-seq peaks. Used for computing background weighted insertion counts.
#' @param motifRanges [GenomicRanges::GRanges-class] object containing coordinates of motif-matches.
#' Needs to contain a metadata column `motif_id` if insertion profiles should be computed for several motifs separately.
#' @param resize If TRUE, peaks will be resized to the provided width around their center. Default is FALSE.
#' @param width If `resize=TRUE`, width to which peaks will be resized.
#' @param minFrag Minimum fragment length to be considered for counting insertions. Default is 30bp.
#' @param maxFrag Maximum fragment length to be considered for counting insertions. Default is
#' @param shift If Tn5 insertion bias should be considered (only if strand column is provided).
#' @param ... Additional arguments passed to `getInsertionProfiles()`, e.g. for computing insertion profiles.
#' @author Emanuel Sonder
#' @import data.table
#' @import GenomicRanges 
#' @import SummarizedExperiment
#' @import Matrix
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @return A list containing the following elements: SummarizedExperiment object with weighted insertion counts for the motif matches (insertioncounts), a SummarizedExperiment object with weighted insertion counts for the peaks (bgcounts) and the insertion profile used for weighting (bgPeakProfile). If `calcProfile=TRUE` also a data.table with the computed insertion profiles per motif is returned (insertion_profiles).
#'
#' @export
getWeightedInsertions <- function(atacData, 
                                  peakRanges,
                                  motifRanges,
                                  genome,
                                  resize=FALSE,
                                  width=300,
                                  minFrag=30, 
                                  maxFrag=3000,
                                  shift=FALSE,
                                  ...){
                                    
  seqLevelStyle <- GenomeInfoDb::seqlevelsStyle(genome)
  if(resize){
    peakRanges <- .resizeRanges(peakRanges=peakRanges, width=width)
  }
  species <- genome@metadata$organism
  peakRanges <- .standardChromosomes(peakRanges, species=species,
                                     genome=genome, coerceToGenome=TRUE)
  motifRanges <- .standardChromosomes(motifRanges, species=species,
                                      genome=genome, coerceToGenome=TRUE)
  peakRanges <- sort(peakRanges)   
  peakRanges$peak_id <- 1:length(peakRanges)  
  GenomeInfoDb::seqlengths(peakRanges) <- GenomeInfoDb::seqlengths(genome)[GenomeInfoDb::seqlevels(peakRanges)]
  peakRanges <- GenomicRanges::trim(peakRanges) 
  if(sum(width(peakRanges)==0)>0 | resize & sum(width(peakRanges)<width)>0){
    warning("Some peaks have width 0 or smaller than the provided width after resizing. These will be removed.")
    peakRanges <- peakRanges[width(peakRanges)>0]
    if(resize) peakRanges <- peakRanges[width(peakRanges)>=width] 
  }       

  # get weighted insertion counts around motif matches
  args <- list(...)
  if("profiles" %in% names(args)) args$calcProfile <- FALSE
  insMot <-  getInsertionProfiles(atacData=atacData, motifRanges=motifRanges,
                                  minFrag=minFrag, maxFrag=maxFrag, shift=shift, ...)

  mmCounts <- insMot$insertion_weighted_counts[,.(w_inserts=sum(w_inserts)), 
                                                  by=.(motif_id, sample, motif_match_id, 
                                                       chr, start, end)]
  
  # aggregate per peak
  max0 <- function(x) {
    if (length(x) == 0) return(0)
    max(x, na.rm=TRUE)
  }
  peakCounts <- genomicRangesMapping(peakRanges, 
                                     assayTable=mmCounts,
                                     byCols=c("motif_id", "sample"), 
                                     scoreCol="w_inserts",
                                     aggregationFun=max0,
                                     type="any")
  samples <- colnames(peakCounts[[1]])
  samples <- samples[order(samples)]
  peakCounts <- lapply(peakCounts, function(x){x[,samples]})
  counts <- lapply(peakCounts, colSums)
  counts <- lapply(counts, Matrix::Matrix, nrow=1)
  counts <- Reduce("rbind", counts[-1], counts[[1]])
  rownames(counts) <- names(peakCounts)
  colnames(counts) <- samples
  
  # get weighted insertions of peaks
  # todo: switch to do.call
  args$calcProfile <- TRUE
  args$profiles <- NULL
  args <- c(list(atacData=atacData, motifRanges=peakRanges, 
                 minFrag=minFrag, maxFrag=maxFrag, shift=shift), args)
  insPeaks <-  do.call(getInsertionProfiles, args)
  bgCounts <- insPeaks$insertion_weighted_counts[,.(w_inserts=sum(w_inserts)), 
                                                    by=.(motif_id, sample, motif_match_id, 
                                                         chr, start, end)]
  zeroPeaks <- setdiff(1:length(peakRanges), bgCounts$motif_match_id)
  bgCounts <- rbind(bgCounts, data.table(motif_match_id=rep(zeroPeaks, each=ncol(counts)),
                                         sample=rep(colnames(counts), each=length(zeroPeaks)),
                                         w_inserts=0),fill=TRUE)
  bgCounts <- dcast(bgCounts, motif_match_id ~ sample, value.var="w_inserts", fill=0, 
                    fun.aggregate=sum, drop=FALSE)
  setorder(bgCounts, motif_match_id)
  bgCounts$motif_match_id <- NULL
  bgCounts <- Matrix::Matrix(as.matrix(bgCounts))
  bgCounts <- bgCounts[,samples]

  # ensure widths (resizing might have an effect)
  peakRanges <- .getGCContent(peakRanges, genome=genome)
  peakRanges$bias <- peakRanges$gc
  peakRanges$gc <- NULL
  bgCounts <- SummarizedExperiment(assays=list(weighted_counts=bgCounts),
                                   rowRanges=peakRanges,
                                   metadata=list(bgPeakProfile=insPeaks$insertion_profile))

  # add flbias to rowData of bgCounts
  if("subSample" %in% names(args)){
    subSample <- args[["subSample"]]
  } else {
    subSample <- NULL
  }
  atacData <- lapply(atacData, .processData, minFrag=minFrag, maxFrag=maxFrag, seqLevelStyle=seqLevelStyle, filterLength=TRUE, shift=shift, subSample=subSample)
  atacData <- rbindlist(atacData, idcol="sample")
  atacData[,frag_width:=as.integer(end-start+1L)]
  peakFragWidths <- genomicRangesMapping(peakRanges, 
                                         assayTable=atacData, 
                                         byCols="sample", 
                                         scoreCol="frag_width",
                                         aggregationFun=function(x){floor(median(x))}, # get mean fragment width per peak and sample, then get the mean across samples (to avoid bias from different sequencing depth across samples)
                                         type="any",
                                         shift=FALSE)
  rowRanges(bgCounts)$flbias <- rowMeans(peakFragWidths)

  counts <- SummarizedExperiment(assays=list(weighted_counts=counts),
                                 rowData=data.frame(motif_id=rownames(counts)),
                                 colData=data.frame(sample=colnames(counts)))

  insCounts <- list(insertioncounts=counts, 
                    bgcounts=bgCounts, 
                    insraw=insMot,
                    bgPeakProfile=insPeaks$insertion_profile)
  if(!is.null(insMot$insertion_profile)){
    insCounts$insertion_profiles <- insMot$insertion_profile
  }

  return(insCounts)
}
