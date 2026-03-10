#' Get chromVAR bin-bin background selection probabilities
#'
#' @param object A chromVAR counts object
#' @param bias Vector of GCBias (default from object)
#' @param w Standard deviation of the Gaussian kernel)
#' @param bs Number of bins per dimension (total bins = bs^2)
#'
#' @return a list with the slots `peak2bin` (which bin each peak belongs to) and
#'   `binBinProbs` (the probability of a peak from a given bin being selected
#'   as background for another).
#' @importFrom stats cov
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom matrixStats colMins colMaxs
getBinProbMatrix <- function(x, bias=NULL, w=0.1, bs=50){
    if (inherits(x, "SummarizedExperiment") ||
            inherits(x, "SingleCellExperiment")) {
        if(is.null(bias)) bias <- rowData(x)$bias
        x <- rowSums(assay(x))
    }else if(is.null(bias)){
        stopifnot("`bias` not provided, and not found in the object.")
    }

    # Mahalanobis transformation
    intensity <- log10(x)
    norm_mat <- matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
    chol_cov_mat <- chol(cov(norm_mat))
    transMat <- t(forwardsolve(t(chol_cov_mat), t(norm_mat)))

    # Calculate min/max for the bin boundaries
    minCoords <- colMins(transMat)
    maxCoords <- colMaxs(transMat)
    range_coords <- maxCoords - minCoords

    # Map peak coordinates to bins
    idx1 <- round((transMat[, 1] - minCoords[1]) /
            range_coords[1] * (bs - 1)) + 1
    idx2 <- round((transMat[, 2] - minCoords[2]) /
            range_coords[2] * (bs - 1)) + 1
    # Ensure indices fall within [1, bs]
    idx1 <- pmax(1, pmin(bs, idx1))
    idx2 <- pmax(1, pmin(bs, idx2))

    # Linearize index
    peak_to_bin <- idx1 + (idx2 - 1) * bs
    bin_density <- tabulate(peak_to_bin, nbins = bs^2)

    # bin center grid (for distance calculation)
    bins1 <- seq(minCoords[1], maxCoords[1], length.out = bs)
    bins2 <- seq(minCoords[2], maxCoords[2], length.out = bs)
    bin_data <- expand.grid(bins1, bins2)

    # bin-to-bin probability matrix
    bin_dist <- as.matrix(dist(bin_data))
    W <- dnorm(bin_dist, 0, w)

    normalizer <- as.vector(W %*% bin_density)
    # Avoid division by zero for empty regions
    normalizer[normalizer == 0] <- 1
    binBinProbs <- W/normalizer

    return(list(
        peak2bin = peak_to_bin,
        binBinProbs = binBinProbs
    ))
}


#' Analytic version of chromVAR's computeDeviations
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', or a count (sparse) matrix.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay.
#' @param bias Per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(object)$bias`.
#' @param w Standard deviation of the Gaussian kernel. Values close to zero will
#'   effectively mean that only peaks from the same bins are used as background,
#'   which is suboptimal if the bin is sparsely populated. High values (e.g. >1)
#'   will lead to homogeneous sampling, which will fail to correct for bias.
#'   A value of 0.05 will be highly similar to the original chromVAR, and values
#'   below 0.2 are recommended.
#' @param bs Number of bins per dimension (total bins = bs^2)
#' @param expectation Optional vector of length equal to `nrow(object)`
#'   giving the expected counts. If NULL, defaults to mean counts.
#' @author Pierre-Luc Germain
#'
#' @return A SummarizedExperiment containing the adjusted deviations and
#'   z-scores for each motif/sample.
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom Matrix crossprod sparseMatrix
computeDeviationsAnalytic <- function(object, annotations, w=0.1, bs=50,
    bias=NULL, expectation=NULL) {

    # Check input validity

    stopifnot(nrow(object) == nrow(annotations))
    stopifnot(is.null(expectation) || length(expectation)==nrow(object))

    motifCD <- NULL
    if( inherits(annotations, "SummarizedExperiment") ){
        motifCD <- colData(annotations)
        annotations <- assay(annotations)
    }

    if(max(annotations) > 1 || min(annotations)<0)
        warning("`annotations` should be either binary or weights from 0 to 1.")

    if( inherits(object, "SummarizedExperiment") ||
            inherits(object, "SingleCellExperiment") ){
        if(is.null(bias)) bias <- rowData(object)$bias
        counts <- assay(object, "counts")
    }else{
        object <- SummarizedExperiment(list(counts=object))
        counts <- assay(object)
    }

    stopifnot(!is.null(bias) && length(bias)==nrow(object))

    if(!is(counts, "matrix") && !inherits(counts, "sparseMatrix"))
        stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
            " or a (sparse) matrix of counts.")

    if(is.null(expectation)){
        expectation <- rowMeans(counts)
        if(any(expectation==0))
            stop("The object contains rows with no counts; please filter them out.")
    }else{
        if(any(expecation==0))
            stop("Some peaks have an expectation of zero; please remove them.")
    }

    # get background bins (B)
    background <- getBinProbMatrix(expectation, bias = bias, w = w, bs = bs)
    bin_map <- background$peak2bin
    binBinProbs <- background$binBinProbs

    # sparse mapping from peaks to bins
    bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation),
        dims=c(nrow(binBinProbs), length(expectation)))

    # bin-level expectations and variances (B x S)
    E_us <- as.matrix(binBinProbs %*% as.matrix(bin2peakMat %*% counts))
    E2_us <- as.matrix(binBinProbs %*% as.matrix(bin2peakMat %*% (counts^2)))
    V_us <- E2_us - (E_us^2)

    # motif-containing peaks per bin (M x B)
    motif_bin_counts <- t(annotations) %*% t(bin2peakMat)

    # motif-level background stats (M x S)
    motif_bg_exp <- as.matrix(motif_bin_counts %*% E_us)
    motif_bg_sd <- sqrt(pmax(0, as.matrix(motif_bin_counts %*% V_us)))

    # observed motif sums (M x S)
    observed_motif_sum <- Matrix::crossprod(annotations, counts)

    # deviation = (Obs - bgExpect) / bgExpect
    deviations <- (observed_motif_sum - motif_bg_exp) / motif_bg_exp

    # Z-score = (Obs - bgExpect) / bgExpect
    z_scores <- (observed_motif_sum - motif_bg_exp) / motif_bg_sd

    SummarizedExperiment(
        assays = list(deviations = deviations, z = z_scores),
        colData = colData(object),
        rowData = motifCD,
        metadata = metadata(object)
    )
}
