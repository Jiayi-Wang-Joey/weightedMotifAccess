#' getInteractionsDeviations
#'
#' Get \code{\link[betterChromVAR]{betterChromVAR}} deviations for overlaps 
#' between motifs and a bait motif.
#' 
#' @param se An object inheriting `SummarizedExperiment`, containing a `counts`
#'   assay.
#' @param annotation A (sparse) matrix or object inheriting 
#'   `SummarizedExperiment` containing motif annotations for each row of `se`.
#' @param bait The name of the bait motif (corresponding to a column of 
#'   `annotation`), for which to search for interacting motifs.
#' @param minCount The minimum overlap size for a motif pair to be considered.
#' @param ... Passed to `betterChromVAR`
#'
#' @returns A `SummarizedExperiment` of deviations for the intersection of 
#'   all motifs with the bait motif.
#' @importFrom betterChromVAR betterChromVAR
#' @export
#' @examples
#' attach(betterChromVAR::getDummyData())
#' dev <- getInteractionsDeviations(counts, motifMatches, bait="motif1",
#'                                  minCount=5)
#' dev
getInteractionsDeviations <- function(se, annotation, bait, minCount=20, ...){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(nrow(se)==nrow(annotation))
  stopifnot(length(bait)==1 && bait %in% colnames(annotation))
  if(inherits(annotation, "SummarizedExperiment")){
    annotation <- assay(annotation)
  }
  moi2 <- annotation & annotation[,bait]
  moi2 <- moi2[,colSums(moi2)>=minCount,drop=FALSE]
  if(ncol(moi2)==0) stop("No interaction passing filtering.")
  ov <- colSums(moi2)
  dev <- betterChromVAR(se, moi2, ...)
  tot <- colSums(annotation[,colnames(moi2)])
  totBait <- as.integer(tot[bait])
  expected <- totBait*tot/nrow(annotation)
  rowData(dev)$overlap.percentOfBait <- round(100*(ov/totBait),2)
  jac <- .colOverlaps(annotation[,colnames(moi2)], annotation[,bait],
                      metric="jaccard")
  rowData(dev)$overlap.jaccard <- as.numeric(jac)
  rowData(dev)$overlap.log2Enr <- log2(ov/expected)
  rowData(dev)$isBait <- row.names(dev)==bait
  dev
}

#' .colOverlaps
#' 
#' Computes pairwise overlap metric between columns (taken from epiwraps)
#'
#' @param x A logical matrix.
#' @param y An optional nother logical matrix or vector. If NULL (default), 
#'   overlaps are computed between columns of `x`.
#' @param metric Either 'overlap', 'jaccard', or 'overlapCoef'.
#'
#' @returns A matrix of pairwise metric values.
.colOverlaps <- function(x, y=NULL, metric=c("overlap","jaccard","overlapCoef")){
  metric <- match.arg(metric)
  x <- as(x, "lMatrix")
  ys <- xs <- colSums(x)
  if(is.null(y)){
    y <- x
  }else{
    if(is.vector(y)) y <- as.matrix(y)
    y <- as(y, "lMatrix")
    ys <- colSums(y)
  }
  intersections <- t(x) %*% as.matrix(y)
  if(metric=="overlap"){
    out <- intersections
  }else if(metric=="overlapCoef"){
    minSize <- outer(xs, ys, "min")
    out <- intersections/minSize
  }else{
    unions <- outer(xs, ys, "+") - intersections
    out <- intersections/unions
  }
  if(is(out, "dgeMatrix")) out <- as.matrix(out)
  out
}


#' discoverMotifInteractions
#' 
#' Discover motifs interacting with a bait motif, based on the betterChromVAR
#' package.
#'
#' @param dev A deviations SummarizedExperiment produced by 
#'   \code{\link{getInteractionsDeviations}}
#' @param group The column of `colData(dev)` containing the experimental groups
#'   of interest. To column will be coerced to a factor, and all comparisons 
#'   will be made relative to the first level of the factor.
#' @param covar Optional columns of `colData(dev)` containing covariates to 
#'   correct for.
#' @param global Logical; whether to test the contrasts globally, rather than
#'   individually. Has no effect when there is a single contrast 
#'   (i.e. two-group comparison).
#' @param weights The weights to use, either 'fixed' (uses the sqrt number of 
#'   sites), 'trend' (default; uses the trend over the number of sites), or 
#'   'none' (no weights; not recommended). We recommend 'trend' if the number
#'   of motifs is sufficiently large (e.g. >100).
#' @param useAssay Which assay to use. We recommend setting 'adjZ' (default),
#'   which uses z-scores but corrects the scale of the bait's z-scores to its
#'   expectation for the intersection size.
#' @author Pierre-Luc Germain
#'
#' @returns A data.frame.
#' @importFrom limma lmFit makeContrasts eBayes topTable contrasts.fit
#' @importFrom stats model.matrix
#' @export
#' @examples
#' attach(betterChromVAR::getDummyData())
#' dev <- getInteractionsDeviations(counts, motifMatches, bait="motif1",
#'                                  minCount=5)
#' dev$group <- rep(LETTERS[1:2], each=5)
#' discoverMotifInteractions(dev, "group")
discoverMotifInteractions <- function(dev, group, covar=c(), global=FALSE,
                                      weights=c("trend","fixed","none"),
                                      useAssay=c("adjZ","deviations","z")){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  stopifnot(sum(rowData(dev)$isBait)==1)
  useAssay <- match.arg(useAssay)
  bait <- row.names(dev)[which(rowData(dev)$isBait)]
  wBait <- which(row.names(dev)==bait)
  N <- rowData(dev)$N
  if(useAssay=="adjZ"){
    useAssay <- "z"
    # shrink the bait's z to expectation for the interaction size
    nf <- sqrt(N/N[wBait])
    e1 <- outer(nf[-wBait], assay(dev,useAssay)[bait,])
  }else{
    e1 <- matrix(assay(dev,useAssay)[bait,], nrow=nrow(dev)-1, ncol=ncol(dev))
  }
  e2 <- assay(dev,useAssay)[-wBait,]
  cd <- as.data.frame(rbind(colData(dev),colData(dev)))
  cd$isInt <- rep(c(FALSE,TRUE), each=ncol(dev))
  cd$combined <- factor(paste(cd[[group]], cd$isInt, sep="_"))
  cd[[group]] <- factor(cd[[group]])
  formula <- paste(c("~0+combined", covar), collapse="+")
  mm <- model.matrix(as.formula(formula), data=cd)
  colnames(mm) <- gsub("combined","",colnames(mm))
  e <- cbind(e1,e2)
  if(weights=="fixed"){
    # use fixed weights
    motif_weight <- sqrt(rowData(dev)$N)
    w <- matrix(motif_weight[-wBait], nrow=nrow(e), ncol=ncol(e))
    w[, seq_len(ncol(dev))] <- motif_weight[wBait]
    fit <- lmFit(e, mm, weights=w)
  }else{
    fit <- lmFit(e, mm)
  }
  refG <- levels(cd[[group]])[1]
  grs <- levels(cd[[group]])[-1]
  names(grs) <- paste0(grs,"-",refG)
  contrast_strings <- sapply(grs, function(g) {
    paste0("(", g, "_TRUE - ", g, "_FALSE) - (", 
           refG, "_TRUE - ", refG, "_FALSE)")
  })
  cont_matrix <- makeContrasts(contrasts=contrast_strings, levels=mm)
  if(weights=="trend"){
    # let limma figure out the prior trend
    motif_weight <- sqrt(rowData(dev)$N[-wBait])
    fit <- eBayes(contrasts.fit(fit, cont_matrix), trend=motif_weight)
  }else{
    fit <- eBayes(contrasts.fit(fit, cont_matrix))
  }
  if(!isTRUE(global)){
    res <- dplyr::bind_rows(lapply(setNames(names(grs),names(grs)), \(co){
      res <- as.data.frame(topTable(fit, coef=co, number=Inf))
      colnames(res)[1] <- "difference"
      cbind(motif=row.names(res), res)
    }), .id="contrast")
  }else{
    res <- as.data.frame(topTable(fit, coef=names(grs), number=Inf))
    res <- cbind(motif=row.names(res), res)
  }
  res$B <- res$AveExpr <- NULL
  cols <- c("N",grep("^overlap", colnames(rowData(dev)), value=TRUE))
  d <- as.data.frame(rowData(dev)[,cols])
  colnames(d)[1] <- "overlap.N"
  cbind(res[,c("contrast","motif")], d[res$motif,],
        res[,setdiff(colnames(res), c("contrast", "motif"))])
}

#' exploreMotifInteraction
#' 
#' We compare deviations in sites that harbor a single of the two motifs, or
#' the pair of motifs at given distances from each other.
#'
#' @param counts An object inheriting from `RangedSummarizedExperiment`, with
#'   a `bias` column in its `rowData`.
#' @param motifs A named pair of motifs, in a format recognized by 
#'   \code{\link[motifmatchr]{matchMotif}}.
#' @param genome A BSgenome or FaFile object.
#' @param nDistQuantiles The number of distance quantiles between the motif 
#'   matches (up to `maxDist`)
#' @param maxDist The maximum distance between matches.
#' @param ... Any further argument passed to 
#'   \code{\link[betterChromVAR]{betterChromVAR}}.
#' @author Pierre-Luc Germain
#'
#' @returns A SummarizedExperiment object of deviations for different types of
#'   interactions between the motifs.
#' 
#' @details
#' Motif pairs are split according to the distance between their centers up to
#' `maxDist` (default 300) into `nDistQuantiles` (default 4) similarly-populated
#' quantiles. Deviations are reported for these different bins, as well as sites
#' harboring only a single motif.
#' @export
#' @importFrom betterChromVAR computeDeviationsAnalytic normalizeDevsForSize
#' @importFrom betterChromVAR betterChromVAR
#' @importFrom motifmatchr matchMotifs
exploreMotifInteraction <- function(counts, motifs, genome, bg=NULL,
                                    nDistQuantiles=5L, maxDist=500L, ...){
  stopifnot(length(motifs)==2)
  stopifnot(inherits(counts, "RangedSummarizedExperiment"))
  moi <- motifmatchr::matchMotifs(motifs, counts, genome=genome)
  moi2 <- motifCoOccurence(motifs, list(names(motifs)), rowRanges(counts),
                           genome = genome, nDistQuantiles=nDistQuantiles,
                           maxDist=maxDist)
  breaks <- round(attr(moi2, "breaks")[[1]][-1])
  bins <- c(paste0("d<=", breaks[1]), paste0(breaks[-length(breaks)],
                                             "<d<=", breaks[-1]))
  moi2flat <- Reduce(cbind, moi2)
  colnames(moi2flat) <- paste(colnames(moi2[[1]]),
                              rep(names(moi2), each=ncol(moi2[[1]])), sep=".")
  moiSinglet <- (assay(moi)-(Reduce("+", lapply(moi2,Matrix::rowSums))>0))>0
  colnames(moiSinglet) <- paste0(colnames(moiSinglet),"+,",
                                 rev(colnames(moiSinglet)),"-")
  moi <- cbind(assay(moi), moiSinglet, moi2flat)
  for(b in seq_along(bins)){
    colnames(moi) <- gsub(paste0("\\.bin",b), paste0("\n", bins[b]),
                          colnames(moi))
  }
  if(is.null(bg)){
    dev <- betterChromVAR(counts, moi, ...)
  }else{
    dev <- computeDeviationsAnalytic(counts, bg, moi, ...)
  }
  normalizeDevsForSize(dev)
}

