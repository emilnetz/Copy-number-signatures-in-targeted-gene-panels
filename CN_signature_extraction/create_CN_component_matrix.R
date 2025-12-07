library(sigminer)

#read in CN seg file 
segTabs <- read.table("/data/gpfs-1/users/netze_c/work/projects/cnv_sigs/data/MH_combined_loh.tsv")


cn <- read_copynumber(segTabs, seg_cols = c("chromosome", "start", "end", "segVal"), 
    genome_measure = "wg", complement = TRUE, add_loh = TRUE, join_adj_seg = FALSE)

tally <- sig_tally(cn, method = "S")

CN48 <- tally$all_matrices$CN_48

saveRDS(CN48, "/data/gpfs-1/users/netze_c/work/projects/cnv_sigs/results/matrices/HRD_CN48_matrix.rds")