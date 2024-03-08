options(warn = 1)

library(data.table)
library(GenomicFeatures)

max_TSS_distance <- 200L # for an anchor to be considered TSS-proximal


#
#  read the manually annotated loops
#

loops <- fread("data/muszka/long_and_short_range_loops_D_mel_manual.tsv", header = T, sep = ",")
# remove old annotations
loops[, x := NULL]
loops[, y := NULL]
loops[, `own_loop_no` := NULL]
loops[, `old_16-18h_WE` := NULL]
loops[, `in_CNS_L3` := NULL]
loops[, `in_WB_adult` := NULL]

# change colum name to more concise ones
for (n in grep("^Dmel_", colnames(loops), value=TRUE))
  setnames(loops, n, sub("^Dmel_OregonR_WT_", "Dmel_", sub("^Dmel_w1118_WT_", "Dmel_", n)))
setnames(loops, "Dmel_larva_brain", "Dmel_larval_brain")
setnames(loops, "dm6_HiC_WE_16-18h_Schuettengruber2014_1000", "Dmel_16-18h_whole_embryo")

setcolorder(loops, c(
  "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "loop_id",
  "Dmel_6-8h_Neuroblasts", "Dmel_6-8h_Neurons", "Dmel_6-8h_Glia",
  "Dmel_10-12h_Neuroblasts", "Dmel_10-12h_Neurons", "Dmel_10-12h_Glia",
  "Dmel_14-16h_Neuroblasts", "Dmel_14-16h_Neurons", "Dmel_14-16h_Glia",
  "Dmel_16-18h_whole_embryo", "Dmel_larval_brain", "Dmel_adult_brain"
  ))


#
#  read the long-range loops from Mohana et al. 2023
#

loops_long <- fread("data/muszka/long_range_loops_D_mel.tsv", header = T, sep = "\t")

# merge 'loops_long' with 'loops' into 'loops_new'
dt <- loops_long[, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color", "loop_id")]
setnames(dt, "color", "color_long")
setnames(dt, "loop_id", "loop_id_long")
loops_new <- merge(loops, dt, by = c("chr1", "x1", "x2", "chr2", "y1", "y2"), all=TRUE)
stopifnot(nrow(loops_new) == nrow(loops))

# copy the published loop_ids to the new data
loops_new[, loop_id := loop_id_long]
loops_new[, loop_id_long := NULL]
loops_new[, color := ifelse(is.na(color_long), "0,0,255", color_long)]
loops_new[, color_long := NULL]

# sort the loops
loops_new[, xm := (x1 + x2) / 2]
loops_new[, ym := (y1 + y2) / 2]
setkey(loops_new, chr1, xm, chr2, ym)
loops_new[, xm := NULL]
loops_new[, ym := NULL]

# add new loop_ids
max_existing_loop_id <- max(as.integer(sub("^L", "", loops_new$loop_id)), na.rm=TRUE)
sel <- is.na(loops_new$loop_id)
loops_new$loop_id[sel] <- paste0("L", seq(from=max_existing_loop_id+1, length=sum(sel)))


#
#  read the transcript database
#

# read gene id to symbol mapping
gtf <- rtracklayer::import("data/muszka/dmel-all-r6.36.gtf.gz")
gene_symbols <- unique(as.data.table(elementMetadata(gtf)[, c("gene_id", "gene_symbol")]))

# all FlyBase transcripts
txdb <- loadDb("data/muszka/dmel-all-r6.36.sqlite")
tx <- transcripts(txdb, columns=c("gene_id", "tx_name", "TXTYPE"))
tx$gene_id <- sapply(tx$gene_id, paste)
tx$gene_symbol <- gene_symbols$gene_symbol[match(tx$gene_id, gene_symbols$gene_id)]

# take only the coding transcripts
coding_tx <- tx[tx$TXTYPE == "mRNA"]

# convert transcripts to TSSes
tss <- resize(coding_tx, width=1, fix='start')


#
#  find the nearest TSS for all loop anchors
#

annotate_nearest_gene <- function(loops)
{
  gr <- with(loops, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))

  nearest_tss <- as.data.table(nearest(gr, tss, select = "all"))

  nearest_tss[, distance := distance(gr[queryHits], tss[subjectHits])]
  nearest_tss[, gene_id := tss[subjectHits]$gene_id]
  nearest_tss[, gene_symbol := tss[subjectHits]$gene_symbol]
  nearest_tss <- nearest_tss[, list(distance = min(distance)), by = list(queryHits, gene_id, gene_symbol)]

  nearest_tss_aggr <- nearest_tss[, list(
    anchor_TSS_proximal = as.integer(min(distance) <= max_TSS_distance),
    anchor_distance_to_TSS = paste(unique(distance), collapse=","),
    anchor_nearest_gene_id = paste(unique(gene_id), collapse=","),
    anchor_nearest_gene_symbol = paste(unique(gene_symbol), collapse=",")
  ), by = queryHits]

  stopifnot(all.equal(nearest_tss_aggr$queryHits, seq_len(nrow(loops))))
  loops$anchor_TSS_proximal <- nearest_tss_aggr$anchor_TSS_proximal
  loops[, anchor_type := ifelse(anchor_TSS_proximal, "P", "I")]
  loops$anchor_distance_to_TSS <- nearest_tss_aggr$anchor_distance_to_TSS
  loops$anchor_nearest_gene_id <- nearest_tss_aggr$anchor_nearest_gene_id
  loops$anchor_nearest_gene_symbol <- nearest_tss_aggr$anchor_nearest_gene_symbol

  anchor_distance_fun <- function(anchor, anchor_midpoint)
  {
    return(anchor_midpoint[anchor == 2] - anchor_midpoint[anchor == 1])
  }
  anchor_distance_aggr <- loops[, list(loop_anchor_distance = anchor_distance_fun(anchor, anchor_midpoint)),
    by = list(loop_id)]
  loops$loop_anchor_distance <- with(anchor_distance_aggr, loop_anchor_distance[match(loops$loop_id, loop_id)])

  return(loops)
}

loops_new_long <- with(loops_new, rbind(
  data.table(seq = seq_along(loop_id), loop_id, anchor = 1L, anchor_chr = chr1, anchor_start = x1, anchor_end = x2),
  data.table(seq = seq_along(loop_id), loop_id, anchor = 2L, anchor_chr = chr2, anchor_start = y1, anchor_end = y2)
  ))
setkey(loops_new_long, seq)
loops_new_long[, seq := NULL]
loops_new_long[, anchor_midpoint := as.integer((anchor_start + anchor_end) / 2)]

loops_new_long_annotated <- annotate_nearest_gene(loops_new_long)


#
#  save the loops with extra annotations added
#

write.table(loops_new, file = "data/muszka/long_and_short_range_loops_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

write.table(loops_new_long_annotated, file = "data/muszka/long_and_short_range_loops_D_mel_annotated.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_proximal <- with(loops_new_long_annotated, data.frame(gene_id =
  unique(do.call(c, strsplit(anchor_nearest_gene_id[anchor_TSS_proximal == 1], ",")))))
write.table(genes_proximal, file = "data/muszka/long_and_short_range_loops_genes_proximal_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_distal <- with(loops_new_long_annotated, data.frame(gene_id =
  unique(do.call(c, strsplit(anchor_nearest_gene_id[anchor_TSS_proximal == 0], ",")))))
write.table(genes_distal, file = "data/muszka/long_and_short_range_loops_genes_distal_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)


#
#  save the gene universe for completeness
#

genes_universe <- data.frame(gene_id = unique(coding_tx$gene_id))
write.table(genes_universe, file = "data/muszka/genes_universe_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

coding_tx_anno <- coding_tx
names(coding_tx_anno) <- with(coding_tx_anno, paste(gene_id, tx_name, sep = "_"))
rtracklayer::export.bed(coding_tx_anno, "data/muszka/genes_universe_D_mel.bed")


#
#  statistics
#

message(nrow(genes_proximal), " genes attributed to promoter-proximal loop anchors (found within +/- ", max_TSS_distance, " bp of a promoter anchor).")
message(nrow(genes_distal), " genes attributed as the closest to promoter-distal loop anchors")
message(nrow(genes_universe), " genes in the universe")
