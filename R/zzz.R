# Suppress R CMD check NOTEs for data.table non-standard evaluation variables
utils::globalVariables(c(
    # data.table list alias
    ".",
    # WeightModel.R
    "bin", "GCBin", "widthBin", "count", "count_bin",
    "freq_bin", "mean_freq_bin", "weight", "logWeight", "peakID",
    "type", "width", "sample",
    # utils.R
    "chr", "col_width", "col_depth", "V1"
))
