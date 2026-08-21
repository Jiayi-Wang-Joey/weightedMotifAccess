# weightedMotifAccess 0.99.0

* Initial Bioconductor submission.
* `getWeightedCounts()`: the `files` and `atacFrag` arguments are replaced by a
  single `atac` argument, which accepts a named list of `data.table`s
  (or `data.frame`s), a named list of `GRanges` or a `GRangesList`, or a named
  character vector of file paths; for a single sample the bare object can be
  passed unwrapped. `.bed`, `.tsv` and `.txt` fragment files (optionally `.gz`-compressed)
  are now read as documented, and file-type detection is anchored on the
  extension and case-insensitive.
