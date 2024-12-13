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





