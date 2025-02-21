source("~/WeightInsertModels/R/utils.R")
source("~/WeightInsertModels/R/WeightModel.R")

importFrags <- function(samplePath,
    chrSelect=NULL,
    maxFragSize=2000,
    verbose=TRUE){
    if(!is.null(chrSelect)){
        fragDt <- .importTabix(samplePath=samplePath, chrSelect=chrSelect)
    }
    else{
        fragDt <- fread(samplePath)
    }
    fragDt[,width:=end-start+1]
    fragDt <- subset(fragDt, width<maxFragSize)
    return(fragDt)
}

.importTabix <- function(samplePath,
    chrSelect){
    which <-  GRanges(chrSelect, IRanges(1,10^8))
    tf <- TabixFile(samplePath)
    frags <- rtracklayer::import(con=tf, format="bed", which=which)
    elementMetadata(frags) <- NULL
    frags <- as.data.table(frags)
    return(frags)
}



#' @param atacFrag a list of data.table containing 
getAverageIntervals <- \(samplePaths, genome, 
    nWidthBins=30, nGCBins=10,
    BPPARAM=BiocParallel::MulticoreParam(workers = 4),
    chunks = TRUE,
    chromsomes = paste0("chr", c(1:22, "X", "Y"))) {
    
    res <- BiocParallel::bplapply(samplePaths, BPPARAM=BPPARAM, \(path) {
        if (chunks) {
            chr_res <- lapply(chromsomes, \(chr) {
                frag <- importFrags(path, chrSelect = chr)
                #frag[,width:=end-start+1]
                gr <- dtToGr(frag)
                gr <- .getGCContent(gr, genome = genome)
                frag <- as.data.table(gr)
                
                # intervals
                widthIntervals <- unique(quantile(frag$width, 
                    probs = seq(0,1,by=1/(nWidthBins))))
                GCIntervals <-  unique(quantile(frag$gc, 
                    probs = seq(0,1,by=1/nGCBins)))
                
                list(width=widthIntervals, gc=GCIntervals)
            })
    
            intervals <- list(width=rowMeans(sapply(chr_res, \(.) .$width)),
                gc=rowMeans(sapply(chr_res, \(.) .$gc)))
            list(width=widthIntervals, gc=GCIntervals)
        } else {
            frag <- importFrags(path, chrSelect = chromsomes)
            frag[,width:=end-start+1]
            gr <- dtToGr(frag)
            gr <- .getGCContent(gr, genome = genome)
            frag <- as.data.table(gr)
            
            # intervals
            widthIntervals <- unique(quantile(frag$width, 
                probs = seq(0,1,by=1/(nWidthBins-2))))
            GCIntervals <-  unique(quantile(frag$gc, 
                probs = seq(0,1,by=1/nGCBins-2)))
            list(width=widthIntervals, gc=GCIntervals)
        }
        
        
    })

    intervals <- list(width=rowMeans(sapply(res, \(.) .$width)),
        gc=rowMeans(sapply(res, \(.) .$gc)))
    intervals$width <- c(0, intervals$width, Inf)
    intervals$gc <- c(0, intervals$gc, 1)
    list(width=widthIntervals, gc=GCIntervals)
}


getAverageBinFreq <- \(samplePaths, genome, 
    intervals, chunk = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4),
    chromosomes = paste0("chr", c(1:22, "X", "Y"))) {
    
    freq <- BiocParallel::bplapply(samplePaths, BPPARAM = BPPARAM, \(path) {
        .processFrags <- \(frag) {
            gr <- dtToGr(frag)
            gr <- .getGCContent(gr, genome = genome)
            frag <- as.data.table(gr)
            frag[, widthBin := cut(width, 
                breaks = intervals$width, 
                include.lowest = TRUE, labels = FALSE)] 
            frag[, GCBin := cut(gc, 
                breaks = intervals$gc, 
                include.lowest = TRUE, labels = FALSE)] 
            frag[, .(count = .N), by = .(widthBin, GCBin)]
        }
        
        if (chunk) {
            dt <- lapply(chromosomes, \(chr) {
                frag <- importFrags(path, chrSelect = chr)
                .processFrags(frag)
            }) |> rbindlist()
        } else {
            frag <- importFrags(path, chrSelect = chromosomes)
            dt <- .processFrags(frag)
        }
        dt[, count_bin := sum(count), by = .(widthBin, GCBin)]
        dt[, bin := .GRP, by = .(widthBin, GCBin)]
        dt
    })
    
    freqDt <- rbindlist(freq, idcol = "sample")
    freqDt[, freq_bin := (count_bin + 1L) / (sum(count_bin) + 1L), by = sample]
    freqDt[, mean_freq_bin := mean(freq_bin), by = bin]
    freqDt[, weight := mean_freq_bin / freq_bin]
    freqDt
}

##################### new ##########################

source("utils.R")
source("WeightModel.R")

importFrags <- function(samplePath,
    chrSelect=NULL,
    maxFragSize=2000,
    verbose=TRUE){
    if(!is.null(chrSelect) & file.exists(paste0(samplePath, ".tbi"))){
        fragDt <- .importTabix(samplePath=samplePath, chrSelect=chrSelect)
    }
    else{
        fragDt <- fread(samplePath)
        setnames(fragDt, colnames(fragDt)[1:3],
            c("seqnames", "start", "end"))
        if(!is.null(chrSelect)) fragDt <- subset(fragDt, seqnames %in% chrSelect)
    }
    fragDt[,width:=end-start+1]
    fragDt <- subset(fragDt, width<maxFragSize)
    return(fragDt)
}

.importTabix <- function(samplePath,
    chrSelect){
    which <-  GRanges(chrSelect, IRanges(1,10^8))
    tf <- TabixFile(samplePath)
    frags <- rtracklayer::import(con=tf, format="bed", which=which)
    elementMetadata(frags) <- NULL
    frags <- as.data.table(frags)
    return(frags)
}



#' @param atacFrag a list of data.table containing 
getAverageIntervals <- \(samplePaths, genome, 
    nWidthBins=30, nGCBins=10,
    BPPARAM=BiocParallel::MulticoreParam(workers = 4),
    chunks = TRUE,
    chromosomes = paste0("chr", c(1:22, "X", "Y")),
    extend=FALSE) {
    
    res <- BiocParallel::bplapply(samplePaths, BPPARAM=BPPARAM, \(path) {
        if (chunks) {
            chr_res <- lapply(chromosomes, \(chr) {
                frag <- importFrags(path, chrSelect = chr)
                #frag[,width:=end-start+1]
                gr <- dtToGr(frag)
                gr <- .getGCContent(gr, genome = genome)
                frag <- as.data.table(gr)
                
                # intervals
                widthIntervals <- unique(quantile(frag$width, 
                    probs=seq(0,1, length.out=nWidthBins+1)))
                GCIntervals <-  unique(quantile(frag$gc, 
                    probs=seq(0,1,length.out=nGCBins+1)))
                
                list(width=widthIntervals, gc=GCIntervals)
            })
            
            intervals <- list(width=rowMedians(sapply(chr_res, \(.) .$width)),
                gc=rowMedians(sapply(chr_res, \(.) .$gc)))
            intervals
        } else {
            frag <- importFrags(path, chrSelect=chromosomes)
            frag[,width:=end-start+1]
            gr <- dtToGr(frag)
            gr <- .getGCContent(gr, genome = genome)
            frag <- as.data.table(gr)
            
            # intervals
            widthIntervals <- unique(quantile(frag$width, 
                probs = seq(0,1,length.out=nWidthBins+1)))
            GCIntervals <-  unique(quantile(frag$gc, 
                probs = seq(0,1,length.out=nGCBins+1)))
            list(width=widthIntervals, gc=GCIntervals)
        }
    })
    
    intervals <- list(width=rowMeans(sapply(res, \(.) .$width)),
        gc=rowMeans(sapply(res, \(.) .$gc)))
    
    if(extend){
        intervals$width <- c(0, intervals$width, Inf)
        if(tail(intervals$gc, 1)<1) intervals$gc <- c(intervals$gc, 1)
        if(head(intervals$gc, 1)>0) intervals$gc <- c(0, intervals$gc)}
    
    
    widthBinsDt <- data.table(lb_width=intervals$width)
    widthBinsDt <- widthBinsDt[,ub_width:=data.table::shift(lb_width, 
        n=1, type="lead")]
    widthBinsDt <- widthBinsDt[complete.cases(widthBinsDt),]
    widthBinsDt$width_bin_id <- 1:nrow(widthBinsDt)
    
    gcBinsDt <- data.table(lb_gc=intervals$gc)
    gcBinsDt <- gcBinsDt[,ub_gc:=data.table::shift(lb_gc, n=1, type="lead")]
    gcBinsDt <- gcBinsDt[complete.cases(gcBinsDt),]
    gcBinsDt$gc_bin_id <- 1:nrow(gcBinsDt)
    
    # get all combinations of GC and width bins (2D-bins)
    binsDt <- data.table(expand.grid(gcBinsDt$gc_bin_id, 
        widthBinsDt$width_bin_id))
    colnames(binsDt) <- c("gc_bin_id", "width_bin_id")
    binsDt <- merge(binsDt, gcBinsDt, by.x=c("gc_bin_id"), by.y=c("gc_bin_id"))
    binsDt <- merge(binsDt, widthBinsDt, by.x=c("width_bin_id"), 
        by.y=c("width_bin_id"))
    binsDt[,bin_id:=paste(gc_bin_id, width_bin_id, sep="_")]
    binsDt
}


getAverageBinFreq <- \(samplePaths, genome, 
    intervals, chunk = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4),
    chromosomes = paste0("chr", c(1:22, "X", "Y"))) {
    
    freq <- BiocParallel::bplapply(samplePaths, BPPARAM = BPPARAM, \(path) {
        if (chunk) {
            fragCount <- lapply(chromosomes, \(chr) {
                frag <- importFrags(path, chrSelect = chromosomes)
                frag$gc_content <- .getGCContent(frag, genome=genome, 
                    contentOnly=TRUE)
                frag <- .mapFragToBin(frag, intervals)
                fragCount <- frag[,.(count_bin=.N),by=bin_id]
            }) |> rbindlist()
            fragCount <- fragCount[,.(count_bin=.N),by=bin_id]
        } else {
            frag <- importFrags(path, chrSelect = chromosomes)
            frag$gc_content <- .getGCContent(frag, genome=genome, 
                contentOnly=TRUE)
            frag <- .mapFragToBin(frag, intervals)
            fragCount <- frag[,.(count_bin=.N),by=bin_id]
        }
        
        fragCount <- merge(intervals, 
            fragCount, by.x="bin_id", 
            by.y="bin_id", all.y=TRUE)
        fragCount[,count_bin:=fifelse(is.na(count_bin), 0, count_bin)]
    })
    
    # compute bin weights
    freqDt <- rbindlist(freq, idcol = "sample")
    freqDt[, freq_bin := (count_bin + 1L) / (sum(count_bin) + 1L), by = sample]
    freqDt[, mean_freq_bin := mean(freq_bin), by = bin_id]
    freqDt[, weight := mean_freq_bin / freq_bin]
    
    return(freqDt)
}

# TODO: Number of merges can be reduced
.mapFragToBin <- function(fragDt, bins, addWeight=FALSE){
    
    fragDt[,gc_1:=gc_content]
    fragDt[,gc_2:=gc_content]
    
    # get GC bin
    gcBins <- bins[,.(lb_gc=first(lb_gc),
        ub_gc=first(ub_gc)), by=gc_bin_id]
    setkey(gcBins, lb_gc, ub_gc)
    fragDt <- foverlaps(fragDt,
        gcBins[,c("lb_gc", "ub_gc", 
            "gc_bin_id"), 
            with=FALSE],
        by.x=c("gc_1", "gc_2"),
        by.y=c("lb_gc", "ub_gc"), 
        mult="first")
    fragDt$gc_1 <- NULL
    fragDt$gc_2 <- NULL
    
    # get fragment length bin
    widthBins <- bins[,.(lb_width=first(lb_width),
        ub_width=first(ub_width)), by=width_bin_id]
    setkey(widthBins, lb_width, ub_width)
    fragDt[,width_1:=width]
    fragDt[,width_2:=width]
    fragDt <- foverlaps(fragDt,
        widthBins[,c("lb_width", "ub_width", 
            "width_bin_id"), 
            with=FALSE],
        by.x=c("width_1", "width_2"),
        by.y=c("lb_width", "ub_width"),
        mult="first")
    fragDt$width_1 <- NULL
    fragDt$width_2 <- NULL
    fragDt[,bin_id:=paste(gc_bin_id, width_bin_id, sep="_")]
    
    if(addWeight){
        fragDt <- merge(fragDt, bins[,c("bin_id", "weight"), with=FALSE], 
            by.x=c("bin_id"), 
            by.y=c("bin_id"), all.x=TRUE)
    }
    
    return(fragDt)
}

getWeightCounts <- function(samplePaths, # named list
    frequencies,
    regions,
    genome,
    chromosomes = paste0("chr", c(1:22, "X", "Y")),
    chunk=TRUE){
    
    res <- BiocParallel::bpmapply(function(path, sampleName, frequencies, regions, genome, 
        chromosomes, chunk=TRUE) {
        
        bins <- subset(frequencies, sample==sampleName)
        if (chunk) {
            weightedCounts <- lapply(chromosomes, \(chr) {
                fragDt <- importFrags(path, chrSelect=chr)
                fragDt$gc_content <- .getGCContent(fragDt, genome=genome, 
                    contentOnly=TRUE)
                fragDt$sample <- sampleName
                weightedFrags <- .mapFragToBin(fragDt, bins, addWeight=TRUE)
                
                # could subset regions to reduce computations
                weightedCounts <- .genomicRangesMapping(regions,
                    weightedFrags, 
                    byCols="sample",
                    scoreCol="weight",
                    seqNamesCol="seqnames",
                    aggregationFun=sum)
            })
            weightedCounts <- Reduce("cbind",weightedCounts[-1], weightedCounts[[1]])
            weightedCounts <- Matrix(rowSums(weightedCounts), ncol=1)
        }
        else{
            fragDt <- importFrags(path, chrSelect=chromosomes)
            fragDt$gc_content <- .getGCContent(fragDt, genome=genome, 
                contentOnly=TRUE)
            wFragDt <- .mapFragToBin(fragDt, bins, addWeight=TRUE)
            wFragDt$sample <- sampleName
            weightedCounts <- .genomicRangesMapping(regions,
                wFragDt, 
                byCols="sample",
                seqNamesCol="seqnames",
                scoreCol="weight", 
                aggregationFun=sum)
        }
    }, samplePaths, names(samplePaths), 
        MoreArgs=list(frequencies=frequencies,
            regions=regions,
            genome=genome, 
            chunk=chunk, 
            chromosomes=chromosomes), 
        BPPARAM=BPPARAM)
    res <- Reduce("cbind", res[-1], res[[1]])
    
    return(res)
}

###### test ######
samples <- list.files("/mnt/plger/plger/DTFAB/fullFrags/GATA1/seq_files/", 
    full.names = TRUE)
intervals <- getAverageIntervals(samplePaths = samples, 
    genome=BSgenome.Hsapiens.UCSC.hg38)
frequencies <- getAverageBinFreq(samples, 
    genome = BSgenome.Hsapiens.UCSC.hg38, intervals)

