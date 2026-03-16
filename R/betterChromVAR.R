#' getBinProbMatrix
#' 
#' Get chromVAR-like bin-bin background selection probabilities
#' 
#' @param object A chromVAR counts object
#' @param bias Vector of GCBias (default from object)
#' @param w Standard deviation of the Gaussian kernel)
#' @param bs Number of bins per dimension (total bins = bs^2)
#' @param pseudo Optional pseudocount to be added. This should not be needed 
#'   with standard workflows.
#' 
#' @return a list with the slots `peak2bin` (which bin each peak belongs to) 
#'   and `binBinProbs` (the probability of a peak from a given bin being 
#'   selected as background for another).
#' @references
#'   Schep, A.N., Wu, B., Buenrostro, J.D., Greenleaf W.J. chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data Nature Methods (2017)
#'   
#' @importFrom stats cov
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom matrixStats colMins colMaxs
#' @export
getBinProbMatrix <- function(x, bias=NULL, w=0.1, bs=50, pseudo=0){
  if (inherits(x, "SummarizedExperiment") || 
      inherits(x, "SingleCellExperiment")) {
    if(is.null(bias)) bias <- rowData(x)$bias
    x <- rowSums(assay(x))
  }else if(is.null(bias)){
    stopifnot("`bias` not provided, and not found in the object.")
  }
  
  # Mahalanobis transformation
  norm_mat <- cbind(log10(x+pseudo), bias)
  cv <- cov(norm_mat)
  diag(cv) <- diag(cv) + pseudo/1000
  transMat <- t(forwardsolve(t(chol(cv)), t(norm_mat)))
  
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
  bin_dist <- dist(bin_data)
  W <- dnorm(as.matrix(bin_dist), 0, w)
  
  normalizer <- as.vector(W %*% bin_density)
  # Avoid division by zero for empty regions
  normalizer[normalizer < 1e-9] <- 1
  binBinProbs <- W/normalizer
  
  tt <- table(binBinProbs>0)
  if( (tt["TRUE"]/sum(tt)) < 0.2 )
    binBinProbs <- as(binBinProbs, "sparseMatrix")
    
  return(list(
    peak2bin = peak_to_bin,
    #binBinDist = bin_dist,
    binBinProbs = binBinProbs
  ))
}


#' betterChromVAR
#' 
#' A fast, analytic implementation of chromVAR's computeDeviations, 
#' additionally enabling balanced expectations and bias shrinkage.
#' 
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', or a count (sparse) matrix.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param grouping An optional factor or vector coercible to a factor indicating
#'   the groupings of the columns of `object`. This is optionally used to 1) 
#'   compute the base expectation such that rare cell types are given as much 
#'   weight as abundant ones, and 2) apply shrinkage (if `shrinkage!="none"`) 
#'   on a per-grouping fashion. In single-cell data, the grouping can for 
#'   instance be the interaction of samples and cell types.
#' @param bias Per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(object)$bias`.
#' @param w Standard deviation of the Gaussian kernel. Values close to zero will
#'   effectively mean that only peaks from the same bins are used as background,
#'   which is suboptimal if the bin is sparsely populated. High values (e.g. >1)
#'   will lead to homogeneous sampling, which will fail to correct for bias.
#'   A value of 0.05 will be highly similar to the original chromVAR, and values
#'   below 0.2 are recommended.
#' @param bs Number of bins per dimension (total bins = `bs^2`).
#' @param sigma Sigma parameter for the 2D smoothing. Ignored unless 
#'   `shrinkage="smooth"`.
#' @param expectation Optional vector of length equal to `nrow(object)` 
#'   giving the expected counts. If NULL, defaults to mean counts (eventually
#'   grouped, see `grouping`).
#' @param nthreads Either an integer scalar indicating the number of threads to
#'   use, or a BiocParallelParam object. This is only used for subsets of the 
#'   steps.
#' @author Pierre-Luc Germain
#' @references
#'   Schep, A.N., Wu, B., Buenrostro, J.D., Greenleaf W.J. chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data Nature Methods (2017)
#' 
#' @return A SummarizedExperiment containing the adjusted deviations and 
#'   z-scores for each motif/sample.
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom Matrix crossprod sparseMatrix kronecker Diagonal cbind2
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpnworkers
#' @export
betterChromVAR <- function(object, annotations, grouping=NULL, bias=NULL, 
                           expectation=NULL, verbose=FALSE, bs=50, sigma=1,
                           nthreads=NULL, w=0.05,
                           shrinkage=c("none", "average", "smooth")){
  
  # Check input validity
  shrinkage <- match.arg(shrinkage)
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
    expectation <- .get_expectation(counts, grouping)
  }
  if(any(expectation==0)){
      stop("Some peaks have an expectation of zero, most likely because they ",
           "have zero counts. Please remove them.")
  }

  if(is.null(grouping)) grouping <- rep(factor("all"), ncol(object))
  grouping <- factor(grouping)
  stopifnot(length(grouping)==ncol(object))
  ngroups <- length(levels(grouping))
                    
  if(is.null(nthreads)){
    BPPARAM <- SerialParam(progress=(verbose && ngroups>1))
  }else if(is.integer(nthreads) && length(nthreads)==1 && nthreads>0L){
    BPPARAM <- MulticoreParam(nthreads, progress=verbose)
  }else if(!inherits(nthreads, "BiocParallelParam")){
    stop("`nthreads` should either be a positive integer, or a ",
         "BiocParallelParam object.")
  }
  
  if(verbose) message("Preparing bias bins")
  
  # get background bins (B)
  background <- getBinProbMatrix(expectation, bias = bias, w = w, bs = bs)
  bin_map <- background$peak2bin
  binBinProbs <- background$binBinProbs
  
  # sparse mapping from peaks to bins
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))

  # motif-containing peaks per bin (M x B)
  motif_bin_counts <- t(annotations) %*% t(bin2peakMat)
  binCounts <- NULL
  
  if(shrinkage != "none"){
    binCounts <- bin2peakMat %*% counts
    cs <- colSums(binCounts)
    
    if(verbose) message("Applying shrinkage")
    il <- split(seq_len(ncol(binCounts)), grouping)
    binCounts <- Reduce(cbind2, bplapply(il, BPPARAM=BPPARAM, \(i){
      binCounts2 <- binCounts[,i]
      cs2 <- cs[i]
      if(shrinkage=="average"){
        # method of moment shrinkage towards per-bin average across cells
        binCounts2 <- shrinkColumnProps(binCounts2)
      }else if(shrinkage=="smooth"){
        # method of moment shrinkage towards cell's 2D-smoothed proportions
        stopifnot(sigma>0)
        G <- .diagKernalMatrix(sqrt(nrow(binCounts2)), sigma=sigma)
        # Create the 2D Kronecker Smoothing Matrix
        G <- Matrix::kronecker(G, G)
        binCounts2 <- shrinkColumnProps(binCounts2,
                                        .fastColNorm(G %*% binCounts2))
      }
      binCounts2 %*% Diagonal(x=cs2)
    }))
    binCounts <- binCounts[,order(unlist(il))]
  }
  
  if(verbose) message("Computing deviations")
  
  
  # bin-level expectations and variances (B x S)
  i <- seq_len(ncol(counts))
  if((nW <- BiocParallel::bpnworkers(BPPARAM))>1 & ncol(counts>5000)){
    chunks <- split(i, cut(i, nW, labels=FALSE))
    res <- bplapply(chunks, BPPARAM=BPPARAM, \(i){
      if(shrinkage=="none"){
        binCounts2 <- bin2peakMat %*% counts[,i]
      }else{
        binCounts2 <- binCounts[,i]
      }
      .getDeviations(binBinProbs, annotations, motif_bin_counts, binCounts2,
                     counts[,i])
    })
    res <- list(deviations = Reduce(cbind2, lapply(res, \(x) x$deviations)),
                z = Reduce(cbind2, lapply(res, \(x) x$z)))
  }else{
    if(is.null(binCounts)) binCounts <- bin2peakMat %*% counts
    res <- .getDeviations(binBinProbs, annotations, motif_bin_counts,
                          binCounts, counts)
  }
  
  SummarizedExperiment(
    assays = res,
    colData = colData(object),
    rowData = motifCD,
    metadata = metadata(object)
  )
}

.getDeviations <- function(binBinProbs, annotations, motif_bin_counts,
                           binCounts, counts){
  # bin-level expectations and variances (B x S)
  E <- binBinProbs %*% binCounts
  V <- (binBinProbs %*% binCounts^2) - (E^2)
  
  # motif-level background stats (M x S)
  motif_bg_exp <- motif_bin_counts %*% E
  motif_bg_sd <- motif_bin_counts %*% V
  if(is(motif_bg_sd, "sparseMatrix")){
    motif_bg_sd@x <- sqrt(pmax(0, motif_bg_sd@x))
    motif_bg_sd <- drop0(motif_bg_sd)
  }else{
    motif_bg_sd <- sqrt(pmax(0, as.matrix(motif_bg_sd)))
  }
  
  # observed motif sums (M x S)
  observed_motif_sum <- Matrix::crossprod(annotations, counts)
  
  # deviation = (Obs - bgExpect) / bgExpect; z = (Obs-exp)/sdExpect
  deviations <- observed_motif_sum - motif_bg_exp
  z_scores <- deviations / motif_bg_sd
  deviations <- deviations / motif_bg_exp
  
  list(deviations=deviations, z=z_scores)
}

.fastColNorm <- function(x, cs=colSums(x)){
  x %*% Diagonal(x = 1/cs)
}
.fastColAgg <- function(x, by){
  by <- factor(by)
  stopifnot(length(by)==ncol(x))
  mapMat <- sparseMatrix(
    i = as.integer(by),
    j = seq_along(by),
    dims = c(length(levels(by)), ncol(x)),
    dimnames = list(levels(by), colnames(x))
  )
  tcrossprod(x, mapMat)
}

# Create a 1D Gaussian kernel matrix
.diagKernalMatrix <- function(n, sigma, sparsify=TRUE){
  x <- seq_len(n)
  dist_mat <- outer(x, x, "-")
  G <- exp(-(dist_mat^2) / (2 * sigma^2))
  G <- G/rowSums(G)
  if(sparsify && n>10){
    G[which(G<median(G))] <- 0
    G <- as(G, "sparseMatrix")
  }
  G
}

.get_expectation <- function(counts, grouping=NULL){
  if(is.null(grouping)) return(rowMeans(counts))
  grouping <- factor(grouping)
  stopifnot(length(grouping)==ncol(counts))
  # compute expectation based on an average of group averages
  agcnt <- .fastColAgg(counts, grouping)
  cs <- colSums(agcnt)
  agcnt <- .fastColNorm(agcnt, cs=cs)*median(cs)
  rowMeans(agcnt)
}


#' shrinkColumnProps
#' 
#' Empirical Bayes shrinkage of a matrix of counts towards a prior proportion 
#' (by default the mean across columns).
#'
#' @param x A matrix of counts, with features as rows and samples as columns.
#' @param shrinkTo A vector (of length equal to `nrow(x)`) of sampling 
#' probabilities to shrink towards, or a matrix (of the same dimensions as `x`) 
#' of such probabilities. If omitted, will be the weighted mean of the columns'
#' relative frequencies.
#' @param var.theo Logical; whether to use theoretical (i.e. binomial) variances
#'  of the proportions, rather than the observed (weighted) variance.
#'
#' @returns A matrix of the same dimensions as `x` representing the shrunk 
#'  column-wise proportions.
#' @export
#'
#' @examples
#' # generate a matrix of 5 sampling (with different total counts) of 20 
#' # features based on the same base frequency :
#' baseFreq <- abs(rnorm(20))
#' baseFreq <- baseFreq/sum(baseFreq)
#' mat <- sapply(c(10,20,30,40,50), \(tot){
#'   rpois(length(baseFreq), baseFreq*tot)
#' })
#' # apply shrinkage and confirm that shrunk proportions are better correlated
#' shrunk_mat <- shrinkColumnProps(mat)
#' mean(cor(shrunk_mat))>mean(cor(mat))
shrinkColumnProps <- function(x, shrinkTo=NULL, var.theo=FALSE) {
  cs <- colSums(x)
  bigTotal <- sum(cs)
  p <- .fastColNorm(x, cs=cs)
  
  if(is.null(shrinkTo)){
    mu <- rowSums(x) / bigTotal
  }else{
    mu <- shrinkTo
  }
  
  if(var.theo){
    # theoretical (i.e. binomial) variances of the proportions
    if(!is.array(mu)){
      pos <- mu*bigTotal
    }else if(is.null(shrinkTo)){
      pos <- x
    }else{
      pos <- mu %*% Diagonal(x=cs)
    }
    v_weighted <- bigTotal*pos*(1-pos)
  }else{
    # weighted row variances
    v_weighted <- rowSums(sweep((p - mu)^2, 2, cs, "*")) / bigTotal
  }
  
  # Estimate M (precision parameter)
  # This formula adjusts the observed weighted variance 
  # by subtracting the expected binomial noise
  M <- (mu * (1 - mu) - v_weighted) /
    (v_weighted - (mu * (1 - mu) / mean(cs)))
  
  # Stability Handling
  # If denominator is negative or zero, variance is too low to estimate M
  M[is.na(M) | !is.finite(M) | M <= 0] <- 1000 
  
  # Back-calculate Alpha and Beta Priors
  alpha <- mu * M
  beta <- (1 - mu) * M
  
  # Shrinkage (Posterior Mean)
  # x_ij is the count, n_j is the column sum
  # Result = (x_ij + alpha_i) / (n_j + alpha_i + beta_i)
  
  # Numerator: mat + alpha (vector added to each column)
  # Denominator: n_j (vector) + M (vector) -> requires a matrix
  den_mat <- sweep(matrix(M, nrow = nrow(x), ncol = ncol(x)), 2, cs, "+")
  
  (x + alpha) / den_mat
}



#' CVnorm: chromVAR-inspired ATAC-seq normalization
#' 
#' Corrects ATAC peak counts by removing the effects of technical biases 
#' (GC/accessibility) using the chromVAR background binning approach and an 
#' optional variance-based bias shrinkage (inspired from the qsmooth package) 
#' to preserve group biological signal.
#'
#' @param object A matrix of counts, or a SummarizedExperiment-like object with
#'   an assay named 'counts'.
#' @param bias Per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(object)$bias`.
#' @param grouping Optional grouping for the baseline expectation (prevents 
#'   bias toward more abundant groups).
#' @param smoothGrouping Optional grouping to determine correction strength. 
#'   If bias is consistent within these groups, correction is reduced.
#' @param bs Number of bins per dimension (total bins = `bs^2`).
#' @param w Standard deviation of the Gaussian kernel for bin smoothing.
#' @param Z Logical; whether to return standardized residuals (Z-scores) 
#'   instead of the (default) corrected counts.
#' @param enforceZeros Logical; whether to enforce that zero counts should 
#'   remain zeroes after correction (ignored if `Z=TRUE`).
#' @author Pierre-Luc Germain
#' @references
#'   - Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'     inferring transcription-factor-associated accessibility from 
#'     single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#'   - Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Corrado
#'     Bravo H (2018). “Smooth quantile normalization.” Biostatistics 19 (2),
#'     doi: 10.1093/biostatistics/kxx028
#'     
#' @return If `object` is a matrix, then a matrix of corrected counts of the 
#'   same dimensions. If `object` is a SummarizedExperiment-like object, then
#'   the object is returned with an extra "corrected" assay.
#' @export
CVnorm <- function(object, bias=NULL, grouping=NULL, smoothGrouping=grouping, 
                   bs=50, w=0.05, Z=FALSE, enforceZeros=TRUE){
  
  # input validity
  if (inherits(object, "SummarizedExperiment") || 
      inherits(object, "SingleCellExperiment")) {
    if(is.null(bias)) bias <- rowData(object)$bias
    counts <- assay(object, "counts")
  } else {
    counts <- object
  }
  stopifnot(!is.null(bias) && length(bias) == nrow(counts))
  stopifnot(is.null(grouping) || length(grouping)==ncol(counts))
  stopifnot(is.null(smoothGrouping) || length(smoothGrouping)==ncol(counts))
  
  # global profile
  expectation <- .get_expectation(counts, grouping)
  if(any(expectation==0)){
    stop("Some peaks have an expectation of zero, most likely because they ",
         "have zero counts. Please remove them.")
  }
  peak_p <- expectation / sum(expectation)
  
  # bias bins
  background <- getBinProbMatrix(expectation, bias=bias, w=w, bs=bs)
  bin_map <- background$peak2bin
  binBinProbs <- background$binBinProbs
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))
  
  # Bin-level observed vs expected
  cs <- colSums(counts)
  binCounts <- bin2peakMat %*% counts
  bin_p_expected <- as.numeric(bin2peakMat %*% peak_p)
  
  # raw bias factors (B x S)
  smooth_obs_p <- binBinProbs %*% .fastColNorm(binCounts, cs=cs)
  smooth_exp_p <- as.numeric(binBinProbs %*% bin_p_expected)
  bias_factor <- smooth_obs_p / smooth_exp_p
  
  # variance-based weighting (qsmooth logic)
  if(!is.null(smoothGrouping) && length(unique(smoothGrouping))>1){
    g <- factor(smoothGrouping)
    n_groups <- length(levels(g))
    
    log_R <- log(as.matrix(bias_factor))
    
    # SST (Total Sum of Squares per bin)
    SST <- rowSums((log_R - rowMeans(log_R))^2)
    
    # SSW (Within-group Sum of Squares)
    group_means <- .fastColAgg(log_R, g) %*% Diagonal(x=1/as.numeric(table(g)))
    # Compute means and expand back to sample dimensions
    group_means_mat <- group_means[, as.integer(g)]
    SSW <- rowSums((log_R - group_means_mat)^2)
    
    # SSB (Between-group Sum of Squares)
    SSB <- pmax(0, SST - SSW)
    
    # Mean Squares (accounting for degrees of freedom)
    MSW <- SSW / (ncol(counts) - n_groups)
    MSB <- SSB / (n_groups - 1)
    
    # Weight: MSW / (MSW + MSB)
    # If MSB is large, bias is consistent within groups -> weight -> 0
    # If MSW is large, bias is noisy/sample-specific -> weight -> 1
    weights <- MSW / (MSW + MSB + 1e-12)
    
    # Shrink the log-bias factors towards based on the weights
    log_R_weighted <- log_R * weights
    bias_factor <- exp(log_R_weighted)
  }
  
  # peak-level correction
  peak_bias_factors <- bias_factor[bin_map, ]
  expected_counts <- (peak_p %*% t(cs)) * peak_bias_factors
  
  if(isTRUE(Z)){
    # Pearson residuals
    out <- (counts - expected_counts) / sqrt(expected_counts)
  } else {
    # residuals + global component to keep values in "count-like" scale
    out <- (counts - expected_counts) + (peak_p %*% t(cs))
    out[which(out<0)] <- 0
    if(isTRUE(enforceZeros)) out[which(counts==0L)] <- 0
  }
  
  if(!inherits(object, "SummarizedExperiment")) return(out)
  
  assay(object, "corrected") <- out
  return(object)
}