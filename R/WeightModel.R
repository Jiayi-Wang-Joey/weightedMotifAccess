#source("~/weightedMotifAccess/R/utils.R")
#source("~/weightedMotifAccess/R/peakWeight.R")
#' @title Define fragment type
#'
#' @description To classify fragments into nucleosome-free, mononucleosome,
#' etc.. according to their lengths
#' @param atacFrag a list of data tables containing the fragments information
#' @param cuts the length interval for fragments classification
#' @return a list of fragments data tables
#' @author Jiayi Wang
#' @noRd
.getType <- function(atacFrag, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(atacFrag, function(dt) {
        dt[,width:=end-start+1]
        # check if max(cuts) covers max(width)
        if (max(dt$width) > cuts[length(cuts)])
            cuts[length(cuts)+1] <- max(dt$width)

        dt[,type:=as.numeric(cut(width, breaks = cuts))]
        if (length(unique(dt$type))==1) {
            warnings("All fragments fell into 1 type")
        } else if(length(unique(dt$type))>=6) {
            stop("Too many types!")
        }
        for (i in unique(dt$type)) {
            dt[, paste0("type_",i) := as.integer(type == i)]
        }
        dt

    })
    res
}

#' @title Internal Fragment Binning
#'
#' @description
#' To assign each fragment to a GC- and FL-bin based on its GC-content and
#' length by quantile binning
#'
#' @param atacFrag a list of data tables containing fragments information of
#' each sample
#' @param nWidthBins and @param nGCBins the number of FL and GC bins,
#' respectively
#' @param genome the reference genome
#' @param BPPARAM Bioconductor parallel backend (defaults to SerialParam)
#' @author Jiayi Wang
#' @import BiocParallel
#' @importFrom data.table rbindlist as.data.table
#' @noRd
.getBins <- function(atacFrag,
    nWidthBins = 30,
    nGCBins = 10,
    genome,
    BPPARAM = BiocParallel::bpparam()) {

    print("start binning")
    fragDts <- BiocParallel::bplapply(atacFrag, BPPARAM = BPPARAM, \(dt) {
        gr <- .dtToGr(dt)
        gr <- .getGCContent(gr, genome = genome)
        dt <- as.data.table(gr)
        setnames(dt, "fragWidth", "width")
        dt
    })

    idCol <- "sample"
    fragDt <- rbindlist(fragDts, idcol=idCol)
    fragDt[, count := 1L]
    fragDt <- fragDt[!is.na(gc)]
    rm(fragDts)

    widthIntervals <- unique(quantile(fragDt$width,
        probs = seq(0,1,by=1/nWidthBins)))
    GCIntervals <-  unique(quantile(fragDt$gc,
        probs = seq(0,1,by=1/nGCBins)))


    fragDt[,widthBin:=cut(width,
        breaks=widthIntervals,
        include.lowest=TRUE, labels=FALSE)]
    fragDt[,GCBin:=cut(gc,
        breaks=GCIntervals,
        include.lowest=TRUE, labels=FALSE)]
    print("end binning")

    fragDt
}

# .getBins <- function(atacFrag,
#     nWidthBins = 30,
#     nGCBins = 10,
#     genome,
#     BPPARAM=BiocParallel::bpparam()) {
#     print("start binning")
#     fragDts <- BiocParallel::bplapply(atacFrag, BPPARAM=BPPARAM, \(dt){
#         gr <- dtToGr(dt)
#         gr <- .getGCContent(gr, genome = genome)
#         dt <- as.data.table(gr)
#         dt
#     })
#
#     idCol <- "sample"
#     fragDt <- rbindlist(fragDts, idcol=idCol)
#
#     rm(fragDts)
#
#     widthIntervals <- unique(quantile(fragDt$width,
#         probs = seq(0,1,by=1/nWidthBins)))
#     GCIntervals <-  unique(quantile(fragDt$gc,
#         probs = seq(0,1,by=1/nGCBins)))
#
#
#     fragDt[,widthBin:=cut(width,
#         breaks=widthIntervals,
#         include.lowest=TRUE, labels=FALSE)]
#     fragDt[,GCBin:=cut(gc,
#         breaks=GCIntervals,
#         include.lowest=TRUE, labels=FALSE)]
#     print("end binning")
#
#     fragDt
# }



#' @title Fragment Weighting for Bias Correction
#' @description
#' Internal function to compute weights based on Fragment Length and GC content
#' bins to align sample frequencies.
#'
#' @param atacFrag A \code{list} or \code{CompressedSplitData.table} of fragments.
#' @param genome A \code{BSgenome} object or character string.
#' @param smooth Logical; whether to apply 2D smoothing to weights.
#' @param nWidthBins,nGCBins Integer; number of bins for FL and GC.
#' @param aRange Numeric; bandwidth for smoothing.
#' @param ... Additional arguments passed to \code{fields::smooth.2d}.
#' @author Jiayi Wang
#' @import data.table
#' @import fields
#' @noRd

.weightFragments <- function (atacFrag,
    genome,
    smooth,
    nWidthBins=30,
    nGCBins=10,
    aRange=0,
    ...) {

    fragDt <- .getBins(atacFrag, genome = genome,
        nWidthBins = nWidthBins, nGCBins = nGCBins)
    fragDt[, bin := .GRP, by = .(widthBin, GCBin)]
    fragDt[,count_bin:=sum(count), by=.(sample, bin)]
    dt <- unique(fragDt[, .(sample, count_bin, bin)], by = c("sample", "bin"))
    dt[,freq_bin:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
    dt[,mean_freq_bin:=mean(freq_bin, na.rm=TRUE),by=bin]
    dt[,weight:=mean_freq_bin/freq_bin]
    fragDt <- fragDt[, !duplicated(names(fragDt)), with = FALSE]
    fragDt <- merge(fragDt, dt[, .(bin, sample, weight)],
        by = c("bin", "sample"), all.x = TRUE)

    if (smooth) {
        fb <- fragDt[,c("sample","GCBin","widthBin","weight")]
        fb <- unique(fb)
        dts <- lapply(split(fb, fb$sample), function(dt) {
            dt[,logWeight:=log2(weight)]
            sm <- smooth.2d(dt$logWeight, x=cbind(dt$widthBin, dt$GCBin),
                surface=FALSE,
                nrow=length(unique(fragDt$widthBin)),
                ncol=length(unique(fragDt$GCBin)),
                aRange=aRange)
            dimnames(sm) <- list(levels(dt$widthBin), levels(factor(dt$GCBin)))
            sm
        })
        tbl <- reshape2::melt(dts)
        names(tbl) <- c("widthBin", "GCBin", "smooth", "sample")
        fragDt <- merge(fragDt,
            tbl,
            by = c("GCBin", "widthBin", "sample"),
            allow.cartesian=TRUE)
        fragDt[,weight:=2^smooth]

    }

    split(fragDt, fragDt$sample)
}

#' @title .weightPeaks
#'
#' @description
#' Weight each peak by normalizing arrays by cyclic loess
#'
#' @param counts peak-level accessibility count matrix
#' @param nSamples the number of peaks used for loess fitting
#' @param family.loess if "gaussian" fitting is by least-squares, and if
#' "symmetric" a re-descending M estimator is used with Tukey's biweight function
#' @param span loess parameter, controls the degree of smoothing
#' @importFrom affy normalize.loess
#' @author Jiayi Wang
#' @noRd

.weightPeaks <- function(counts,
    nSample=1e4,
    nBins=100,
    family.loess="symmetric",
    span=0.3) {

    if(nrow(counts) <= nSample) {
        res <- affy::normalize.loess(counts+1L, span=span,
            subset=1:nrow(counts), family.loess=family.loess)
    } else {
        allAvg <- log2(rowMeans(counts+1L))
        bins <- cut(allAvg,
            breaks = seq(min(allAvg), max(allAvg), (max(allAvg)-min(allAvg))/nBins),
            include.lowest = TRUE)
        idx <- unlist(sapply(seq_len(length(levels(bins))), \(i) {
            ids <- which(bins==levels(bins)[i])
            if (length(ids) > 1)
                sample(ids,
                    size = min(length(ids), round(nSample/nBins)))
            else if (length(ids)==1)
                ids
            else NULL
        }))
        res <- affy::normalize.loess(counts+1L, span=span,
            subset=idx, family.loess=family.loess)
        }

}



.getOverlapCounts <- function(peakRanges,
    atacFrag,
    fragWeight = FALSE,
    peakWeight = FALSE,
    insertWeight =FALSE,
    cuts=c(0,120,300,500),
    genome,
    overlap = c("any", "start", "end", "within", "equal"),
    smooth = FALSE,
    aRange = 1.5,
    nWidthBins = 30,
    nGCBins = 10,
    motifRanges = NULL,
    profiles = NULL,
    ...
) {

    peaks <- data.table::as.data.table(peakRanges)
    peaks$peakID <- seq_len(nrow(peaks))

    if (fragWeight) {
        atacFrag <- .weightFragments(atacFrag,
            genome = genome,
            smooth = smooth,
            nGCBins = nGCBins,
            nWidthBins = nWidthBins,
            aRange = aRange)
    }

    atacFrag <- .getType(atacFrag, cuts = cuts)
    fragCounts <- lapply(atacFrag, function(frag) {
        types <- names(frag)[grepl("^type_", names(frag))]
        if (!fragWeight) frag[, count := 1L]
        if (fragWeight) {
            frag[,count:=weight*count]
            frag[,(types) := lapply(.SD, function(x) x*weight),
                .SDcols = types]
        }
        fragGR <- .dtToGr(frag)
        hits <- findOverlaps(fragGR, peakRanges, type = overlap)
        overlaps <- cbind(frag[queryHits(hits),],
            peaks[subjectHits(hits), c("peakID")])
        tmp <- overlaps[, c(list(counts = sum(count, na.rm = TRUE)),
            list(mean_width = mean(width, na.rm = TRUE)),
            list(median_width = median(width, na.rm = TRUE)),
            lapply(.SD, sum, na.rm = TRUE)),
            by = peakID, .SDcols = c(types)]

        res <- data.table(peakID = seq_len(nrow(peaks)))
        res <- merge(res, tmp, all =TRUE)
        res[is.na(res)] <- 0
        res

    })

    cols <- names(fragCounts[[1]])[grepl("^type_|counts|width",
        names(fragCounts[[1]]))]
    allCounts <- lapply(cols, function(x) {
        lst <- lapply(fragCounts, function(.) data.frame(.)[,x])
        mat <- do.call(cbind, lst)
        colnames(mat) <- names(fragCounts)
        mat
    })
    names(allCounts) <- cols
    if (peakWeight) {
        allCounts[["counts"]] <- .weightPeaks(allCounts[["counts"]])
        allCounts[["type_1"]] <- .weightPeaks(allCounts[["type_1"]])
    }
    return(SummarizedExperiment(assays = allCounts, rowRanges = peakRanges))
}


#' @title Count (Weighted) Fragments in Peak Regions
#'
#' @description
#' Count the number of fragments in each peak region, with or without
#' fragment-level and peak-level weighting (normalization).
#'
#' @param files A character vector of paths to fragment files (e.g., BED or TSV).
#' @param atacFrag A \code{list} of \code{data.table} objects or a \code{GRangesList}.
#' @param ranges A \code{GRanges} or \code{data.table} object containing peak regions.
#' @param genome A \code{BSgenome} object or string (e.g., "hg38") for GC content.
#' @param species Character; species name (e.g., "human") for chromosome filtering.
#' @param fragWeight Logical; whether to apply fragment-level bias correction.
#' @param peakWeight Logical; whether to apply cyclic loess normalization on counts.
#' @param resize Logical; whether to resize peak ranges to a fixed width.
#' @param width Integer; width to resize peaks to if \code{resize=TRUE}.
#' @param nWidthBins,nGCBins Integer; number of bins for fragment weighting.
#' @param cuts Numeric vector; fragment length thresholds for classification.
#' @param minFrag,maxFrag Integer; fragment length filters.
#' @param smooth Logical; whether to apply smoothing on fragment weights.
#' @param aRange Numeric; bandwidth for smoothing.
#' @param ... Additional arguments passed to internal weighting functions.
#' @author Jiayi Wang
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' containing assays for different fragment types.
#' @import SummarizedExperiment
#' @importFrom fields smooth.2d
#' @importFrom data.table as.data.table
#' @importFrom GenomicRanges findOverlaps GPos resize GRanges
#' @importFrom Rsamtools FaFile
#'
#' @export
#'
getWeightedCounts <- function(
        files,
    atacFrag,
    ranges,
    genome,
    species,
    fragWeight = FALSE,
    peakWeight = FALSE,
    resize = FALSE,
    width = 300,
    nWidthBins = 30,
    nGCBins = 10,
    cuts = c(0, 120, 300, 500),
    minFrag = 30,
    maxFrag = 3000,
    smooth = FALSE,
    aRange = 1.5,
    coerceToGenome = TRUE,
    ...
) {
    if (!inherits(genome, c("BSgenome", "FaFile"))) {
        stop(
            "The 'genome' argument must be either a BSgenome object or an Rsamtools::FaFile."
        )
    }

    if (is.null(atacFrag)) {
        atacFrag <- .importFragments(files)
    }

    .sanityCheck(atacFrag, ranges)

    ranges <- .standardChromosomes(
        ranges,
        species = species,
        genome = genome,
        coerceToGenome = coerceToGenome
    )

    atacFrag <- lapply(atacFrag, function(dt) {
        gr <- .dtToGr(dt)
        gr <- .standardChromosomes(
            gr,
            species = species,
            genome = genome,
            coerceToGenome = coerceToGenome
        )
        as.data.table(gr)
    })

    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)

    res <- .matchSeqlevels(atacFrag, ranges)
    atacFrag <- res$atacFrag
    ranges <- res$ranges
    seqlengths(ranges) <- seqlengths(genome)[seqlevels(ranges)]

    if (resize) {
        ranges <- .resizeRanges(peakRanges = ranges, width = width)
    }

    asy <- .getOverlapCounts(
        peakRanges = ranges,
        atacFrag = atacFrag,
        fragWeight = fragWeight,
        cuts = cuts,
        genome = genome,
        smooth = smooth,
        aRange = aRange,
        species = species,
        nWidthBins = nWidthBins,
        nGCBins = nGCBins,
        peakWeight = peakWeight,
        ...
    )

    asy
}




