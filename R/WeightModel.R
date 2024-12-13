source("~/WeightInsertModels/R/utils.R")
source("~/WeightInsertModels/R/peakWeight.R")
#' @description To classify fragments into nucleosome-free, mononucleosome, 
#' etc.. according to their lengths
#' @param atacFrag a list of data tables containing the fragments information
#' @param cuts the length interval for fragments classification
#' @return a list of fragments data tables 

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

#' .getBins
#' 
#' @description
#' To assign each fragment to a GC- and FL-bin based on its GC-content and 
#' length by quantile binning
#' @param atacFrag a list of data tables containing fragments information of 
#' each sample
#' @param nWidthBins and @param nGCBins the number of FL and GC bins, 
#' respectively
#' @param genome the reference genome

.getBins <- function(atacFrag, 
    nWidthBins = 30, 
    nGCBins = 10, 
    genome,
    BPPARAM=BiocParallel::SerialParam()) {
    print("start binning")
    fragDts <- BiocParallel::bplapply(atacFrag, BPPARAM=BPPARAM, \(dt){
        dt[,width:=end-start+1]
        gr <- dtToGr(dt)
        gr <- .getGCContent(gr, genome = genome)
        dt <- as.data.table(gr)
        dt
    })
    fragDt <- rbindlist(fragDts, idcol = 
            if ("barcode" %in% colnames(fragDts[[1]])) NULL else "sample")
    
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

#' moderateBinFrequencies
#'
#' @param bins A vector of FL/GC bin labels for each element of `counts`
#' @param samples A vector of sample labels for each element of `counts`
#' @param counts A vector of counts per sample/bin
#'
#' @return A vector of moderated frequencies for each element of `counts`
#' @author Pierre-Luc
moderateBinFrequencies <- function (bins, samples, counts) {
    totPerSamp <- tapply(counts, samples, sum)
    binFreq <- counts/totPerSamp[samples]
    mv <- aggregate(binFreq, by=list(bin=bins), FUN=function(x){
        c(mu=mean(x,na.rm=TRUE), v=var(x, na.rm=TRUE))
    })
    mv$mu <- mv$x[,1]
    mv$v <- mv$x[,2]
    mv$alpha <- ((1-mv$mu)/mv$v - 1/mv$mu)*mv$mu^2
    mv$beta <- mv$alpha*(1/mv$mu-1)
    (counts + mv$alpha[bins])/(totPerSamp[samples] + mv$alpha[bins] + mv$beta[bins])
}  


#' .weightFragments
#' 
#' @description
#' Give a weight to each fragment such that the fragment frequency of each sample
#' aligns with the mean frequency.
#' @param atacFrag a list of data tables containing fragments information of 
#' each sample
#' @param genome the reference genome
#' @param nWidthBins and @param nGCBins the number of FL and GC bins, 
#' respectively
#' @param aRange bandwidth in smoothing, controls the degree of smoothing
#' @param moderating whether moderated fragment frequencies is used 
.weightFragments <- function (atacFrag, 
    genome,
    smooth,
    nWidthBins,
    nGCBins,
    aRange,
    moderating,
    ...) {
    print("start weighting")
    fragDt <- .getBins(atacFrag, genome = genome, 
        nWidthBins = nWidthBins, nGCBins = nGCBins)
    fragDt[, bin := .GRP, by = .(widthBin, GCBin)]
    fragDt[,count_bin:=sum(count), by=.(sample, bin)]
    dt <- unique(fragDt[, .(sample, count_bin, bin)], by = c("sample", "bin"))
    
    if (moderating) {
        dt$freq_bin <- moderateBinFrequencies(dt$bin, dt$sample, dt$count_bin)
    } else {
        dt[,freq_bin:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
    }
    dt[,mean_freq_bin:=mean(freq_bin, na.rm=TRUE),by=bin]
    dt[,weight:=mean_freq_bin/freq_bin]
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
    print("end weighting")
    if ("barcode" %in% colnames(fragDt)) {
        split(fragDt, fragDt$barcode)
    } else {
        split(fragDt, fragDt$sample)
    }
    #return(split(fragDt, fragDt$sample))
}

#' .weightPeaks
#'
#' @description
#' Weight each peak by normalizing arrays by loess
#' @param counts peak-level accessibility count matrix
#' @param nSamples the number of peaks used for loess fitting
#' @param family.loess if "gaussian" fitting is by least-squares, and if 
#' "symmetric" a re-descending M estimator is used with Tukey's biweight function
#' @param span loess parameter, controls the degree of smoothing
.weightPeaks <- function(counts, 
    nSample=1e4, 
    nBins=100,
    family.loess="symmetric",
    span=0.3,
    peakWeight) {
    if (peakWeight=="loess") {
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
    } else if (peakWeight=="binnedLoess") {
        e <- log2(1L+cpm(calcNormFactors(DGEList(counts))))
        offset <- binnedLoess(e)    
        res <- 2^(e-offset)
    } else if(peakWeight=="monotonicSmoothed") {
        e <- log2(1L+cpm(calcNormFactors(DGEList(counts))))
        offset <- monotonicSmoothed(e)    
        res <- 2^(e-offset)
    } else if (peakWeight=="smoothedBinLm") {
        e <- log2(1L+cpm(calcNormFactors(DGEList(counts))))
        offset <- smoothedBinLm(e)    
        res <- 2^(e-offset)
    }
    
    res
}
    

#' @description
#' count the number of fragments in each peak region
#' 
#' @param peakRanges a GRange or data.table object that contains the peak ranges
#' @param atacFrag a list of data.table that contains the fragment ranges, 
#' each data.table represents a sample
#' @param mode whether to 
#' @return a list of count table; by all, nucleosome-free, mono...
.getOverlapCounts <- function(peakRanges, 
    atacFrag,
    fragWeight = FALSE,
    cuts,
    genome,
    overlap = c("any", "start", "end", "within", "equal"),
    smooth = FALSE, 
    aRange,
    nWidthBins = 30,
    nGCBins = 10,
    peakWeight = c("none", "loess", "binnedLoess"),
    moderating = FALSE,
    singleCell = FALSE,
    ...
) {
    
    peakWeight <- match.arg(peakWeight, 
        choices = c("none", "loess", "binnedLoess", 
            "monotonicSmoothed", "smoothedBinLm"))
    peaks <- data.table::as.data.table(peakRanges)
    peaks$peakID <- seq_len(nrow(peaks))
    
    if (fragWeight) {
        atacFrag <- .weightFragments(atacFrag, 
            genome = genome, 
            smooth = smooth,
            nGCBins = nGCBins,
            nWidthBins = nWidthBins,
            aRange = aRange,
            moderating = moderating,
            singleCell = singleCell)
    }
    
    atacFrag <- .getType(atacFrag, cuts = cuts)
    
    fragCounts <- lapply(atacFrag, function(frag) {
        types <- names(frag)[grepl("^type_", names(frag))]
        if (fragWeight) {
            frag[,count:=weight*count] # remove
            frag[,(types) := lapply(.SD, function(x) x*weight), 
                .SDcols = types]
        } 
        fragGR <- dtToGr(frag)
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
    
    cols <- names(fragCounts[[1]])[grepl("^type_|counts|_width", 
        names(fragCounts[[1]]))]
    allCounts <- lapply(cols, function(x) {
        lst <- lapply(fragCounts, function(.) data.frame(.)[,x])
        mat <- do.call(cbind, lst)
        colnames(mat) <- names(fragCounts)
        mat
    })
    names(allCounts) <- cols
    if (peakWeight != "none") {
        allCounts[["counts"]] <- .weightPeaks(allCounts[["counts"]], 
            peakWeight=peakWeight)
        allCounts[["type_1"]] <- .weightPeaks(allCounts[["type_1"]], 
            peakWeight=peakWeight)
    }
    allCounts
}

    
getCounts <- function (files,
    atacFrag,
    ranges, # motif matches or peaks, set a parameter to define
    genome,
    species,
    fragWeight = FALSE,
    paired = TRUE,
    resize = TRUE, 
    width = 300,
    nWidthBins = 30,
    nGCBins = 10,
    cuts = c(0,120,300,500),
    minFrag = 30,
    maxFrag = 3000,
    smooth = FALSE,
    aRange = 0,
    peakWeight=c("none", "loess"),
    moderating = FALSE,
    singleCell = FALSE,
    ...) {
    
    peakWeight <- match.arg(peakWeight, 
        choices = c("none", "loess", "binnedLoess", 
            "monotonicSmoothed", "smoothedBinLm"))
    
    # get fragments ranges
    if (is.null(atacFrag)) atacFrag <- .importFragments(files)
    
    # sanity check
    .sanityCheck(atacFrag, ranges)
    
    # standard chromosomes
    ranges <- .standardChromosomes(ranges, species = species)
    atacFrag <- lapply(atacFrag, function(dt) {
        gr <- dtToGr(dt)
        gr <- .standardChromosomes(gr, species = species)
        as.data.table(gr)
    })
    
    # filter too short or too long fragments
    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)
    
    # match seqLevels
    res <- .matchSeqlevels(atacFrag, ranges)
    atacFrag <- res$atacFrag
    ranges <- res$ranges
    
    
    if (resize) ranges <- .resizeRanges(peakRanges = ranges, width = width)
    
    asy <- .getOverlapCounts(peakRanges = ranges, 
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
            moderating = moderating,
            singleCell = singleCell)

    if (singleCell) {
        cd <- rbindlist(lapply(names(atacFrag), \(x) 
            data.table(sample_id=atacFrag[[x]]$sample[1],
                group_id=atacFrag[[x]]$motif[1],
                barcode=x)))
        cd <- cd[match(cd$barcode, colnames(asy[[1]]))]
        SummarizedExperiment(assays = asy, 
            rowRanges = ranges,
            colData = DataFrame(cd), checkDimnames=FALSE)
    } else {
        SummarizedExperiment(assays = asy, rowRanges = ranges)  
    }
    
}    
    
    
    


