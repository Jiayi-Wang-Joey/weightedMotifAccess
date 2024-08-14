#' @description
#' Resize the peaks
#' 
#' @param peakRanges: a GRanges object of peak ranges
#' @param width: the re-defined size of each peak
#' @param fix: the fixed point for resizing
#' @return a GRange object with resized ranges

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

#' @description To calculate the GC content of each fragment or peak
#' @param gr: a GRanges object 
#' @param genome: a BSgenome object, the corresponding genome 
#' @return a GRanges objects with an additional metadata column gc that contains
#' GC content

.getGCContent <- function(gr, genome) {
    # Sanity check
    if (!class(gr) == "GRanges") {
        stop("peakRanges must be a GRanges object")
    }
    seqs <- Biostrings::getSeq(x = genome, gr)
    mcols(gr)$gc <- letterFrequency(seqs, "GC",as.prob=TRUE)[,1]
    gr
} 

#' @description remove fragments that are too short or too long
#' @param atacFrag: a list of data tables containing the fragments information
#' @param min and @param max the minimum and maximum limit of fragment length
#' @return a list of data tables of filtered fragments 
.filterFrags <- function(atacFrag, min = 30, max = 2000) {
    res <- lapply(atacFrag, function(frag) {
        frag[,width:=end-start+1]
        frag[width>=min & width<=max,]
        frag
    })
    res
}

#' @description match the chromosomes between fragments and peaks/motifs
#' @param atacFrag: a list of data tables containing the fragments information
#' @param ranges: a genomic object of peaks/motifs

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

#' @description
#' Convert a data table to GenomicRange object

dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end"){
    setnames(dt, seqCol, "seqnames", skip_absent = TRUE)
    gr <- GRanges(seqnames=dt[["seqnames"]], ranges=IRanges(start=dt[[startCol]], 
        end=dt[[endCol]]))
    mcols(gr)$count <- dt$count
    if (!is.null(dt$sample)) {
        mcols(gr)$sample <- dt$sample
    }
    
    if (!is.null(dt$motif)) {
        mcols(gr)$motif <- dt$motif
    }
    
    if (!is.null(dt$barcode)) {
        mcols(gr)$barcode <- dt$barcode
    }
    
    if(startCol==endCol)
    {
        gr <- GPos(seqnames=dt[["seqnames"]], pos=dt[[startCol]])
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

#' @description
#' Sanity check to ensure the input arguments have the correct classes
#' 
#' @param atacFrag: a list of data.tables that contain the ranges of fragments
#' @param peakRanges: a GRange object that contains the ranges of peaks
#' @param motifRanges:  a GRange object that contains the ranges of motifs, 
#' check metadata columns
#' check seqnames to factor in datatable
.sanityCheck <- function(atacFrag, 
    ranges) {
    lapply(atacFrag, function(frag) {
        if (!is.data.table(frag)) {
            stop("Each element in atacFrag should be a data.table")
        }
    })
    
    if (!class(ranges)=="GRanges") {
        stop("ranges should be a GRanges object")
    }
    
    
}

