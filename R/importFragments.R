#' @title Normalize fragment input
#'
#' @description
#' Coerces the various accepted forms of the \code{atac} argument into a named
#' \code{list} of \code{data.table}s with \code{seqnames}/\code{start}/\code{end}
#' columns, which is what the rest of the pipeline expects.
#'
#' @param atac A named \code{list} of \code{data.table}s/\code{data.frame}s, a
#'   named \code{list} of \code{GRanges} or a \code{GRangesList}, or a named
#'   character vector of file paths. For a single sample, the bare object (one
#'   \code{data.table}, \code{GRanges} or path) may be passed unwrapped.
#' @return A named \code{list} of \code{data.table}s.
#' @author Jiayi Wang
#' @importFrom data.table as.data.table fread setnames
#' @noRd
.importFragments <- function(atac) {
    if (is.character(atac)) atac <- as.list(atac)
    if (is(atac, "GRangesList")) atac <- as.list(atac)
    # Single-sample shorthand. This must stay ahead of the is.list() check
    # below: a data.table/data.frame is itself a list, so testing is.list()
    # first would treat one sample's columns as several samples.
    if (is(atac, "GRanges") || is.data.frame(atac)) atac <- list(atac)
    if (!is.list(atac)) {
        stop("'atac' should be a named list of data.tables or GRanges, a ",
            "GRangesList, or a named character vector of file paths.")
    }
    if (length(atac) == 0L) stop("'atac' is empty")

    if (is.null(names(atac))) names(atac) <- .fragNames(atac)

    lapply(atac, .importFragment)
}

# Derive sample names when the user did not supply any: file basenames without
# their (possibly compressed) extension, or sample1..N otherwise.
.fragNames <- function(atac) {
    if (all(vapply(atac, \(x) is.character(x) && length(x) == 1L, logical(1)))) {
        nms <- sub("\\.(gz|bgz)$", "", basename(unlist(atac)), ignore.case=TRUE)
        nms <- sub("\\.[^.]+$", "", nms)
        if (!anyDuplicated(nms)) return(nms)
    }
    paste0("sample", seq_along(atac))
}

# Extensions recognized as tabular fragment files, optionally compressed.
.fragTabularPattern <- "\\.(bed|tsv|txt)(\\.(gz|bgz))?$"

# Number of leading '#' comment lines, which some tools prepend to their
# fragment files. Left in place they confuse fread()'s header detection, which
# then warns that the file is invalid while still returning the right columns.
# gzfile() reads plain and gzipped files alike.
.nCommentLines <- function(path) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    n <- 0L
    while (length(l <- readLines(con, n=1L)) > 0L && startsWith(l, "#"))
        n <- n + 1L
    n
}

# Coerce a single element of 'atac' into a data.table of fragments.
.importFragment <- function(x) {
    if (is(x, "GRanges")) return(.checkFragCols(as.data.table(x)))
    if (is.data.frame(x)) return(.checkFragCols(as.data.table(x)))
    if (!(is.character(x) && length(x) == 1L)) {
        stop("Each element of 'atac' should be a data.table, a GRanges, or a ",
            "single file path; got an object of class '", class(x)[1], "'.")
    }

    f <- basename(x)
    if (grepl("\\.bam$", f, ignore.case=TRUE)) {
        .checkFragCols(.readBamFragments(x))
    } else if (grepl(.fragTabularPattern, f, ignore.case=TRUE)) {
        dt <- fread(x, skip=.nCommentLines(x), select=1:3,
            col.names=c("seqnames", "start", "end"))
        .checkFragCols(dt)
    } else if (grepl("\\.rds$", f, ignore.case=TRUE)) {
        # may hold either a GRanges or a data.table/data.frame
        .importFragment(readRDS(x))
    } else {
        stop("Cannot determine the format of '", x, "'. Fragment files should ",
            "be .bam, .bed, .tsv, .txt (optionally .gz) or .rds.")
    }
}

# Read fragment coordinates from the read pairs of an (indexed) bam file, and
# apply the usual ATAC shift.
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom Rsamtools ScanBamParam
.readBamFragments <- function(bamPath) {
    message("Importing bam file: ", basename(bamPath))
    param <- Rsamtools::ScanBamParam(what=c("pos", "qwidth", "isize"))
    readPairs <- GenomicAlignments::readGAlignmentPairs(bamPath, param=param)
    r1 <- GenomicAlignments::first(readPairs)
    r2 <- GenomicAlignments::second(readPairs)
    frags <- GRanges(seqnames(r1),
        IRanges(start=pmin(GenomicAlignments::start(r1),
                GenomicAlignments::start(r2)),
            end=pmax(GenomicAlignments::end(r1),
                GenomicAlignments::end(r2))))
    frags <- granges(frags, use.mcols=TRUE)

    # ATAC shift
    start(frags) <- start(frags) + 4L
    end(frags) <- end(frags) - 4L
    as.data.table(frags)
}

# Fragments need genomic coordinates; 'chr' is accepted as an alias for
# 'seqnames' (as in .processData()). Only the coordinates are kept: everything
# downstream recomputes width, strand is unused, and any count column would be
# overwritten by the weighting anyway. Keeping them around otherwise collides
# with the 'width' column that .dtToGr()/as.data.table() reintroduce.
.checkFragCols <- function(dt) {
    if (!("seqnames" %in% names(dt)) && "chr" %in% names(dt))
        setnames(dt, "chr", "seqnames")
    missing <- setdiff(c("seqnames", "start", "end"), names(dt))
    if (length(missing) > 0L) {
        stop("Fragments are missing the column(s): ",
            paste(missing, collapse=", "),
            ". They should contain seqnames (or chr), start and end.")
    }
    dt[, c("seqnames", "start", "end"), with=FALSE]
}
