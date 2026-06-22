# Suppress R CMD check NOTEs for data.table non-standard evaluation variables
utils::globalVariables(c(
    # data.table list alias
    ".",
    # WeightModel.R
    "bin", "GCBin", "widthBin", "count", "count_bin",
    "freq_bin", "mean_freq_bin", "weight", "logWeight", "peakID",
    "type", "width", "sample",
    # utils.R
    "chr", "col_width", "col_depth", "V1",
    # InsertModel.R
    "rel_pos", "insert", "motif_center", "ml", "motif_id",
    "strand_insert", "w_inserts", "motif_match_id", "frag_width",
    "pos_count_global", "med_pos_count_global", "w", "w_smooth",
    "pos_count", "dev", "end_margin", "start_margin"
))
