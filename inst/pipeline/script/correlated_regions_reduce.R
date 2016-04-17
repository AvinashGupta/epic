
suppressPackageStartupMessages(library(GetoptLong))

cutoff = 0.01
gap1 = "bp(1000)"
gap2 = "bp(1000)"
GetoptLong("config=s", "configuration R script",
	       "cutoff=f", "cutoff for filter cr",
	       "gap1=s", "gap for neg_cr",
	       "gap2=s", "gap for pos_cr")

library(epic)
load_config(config)


files = dir(qq("@{OUTPUT_DIR}/rds"), pattern = "^cr_filtered_fdr_.*\\.rds$")
if(length(files) == 1) {
	cutoff = gsub("^cr_filtered_fdr_(.*)\\.rds$", "\\1", files[1])
	cr_filtered = readRDS(files[1])
} else {
	cr_filtered = readRDS(qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))
}

reduce_cr_gap_test(cr_filtered)

neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]

neg_cr_reduced = reduce_cr(neg_cr, EXPR, TXDB, gap = eval(parse(text = gap1)), mc.cores = 4)
pos_cr_reduced = reduce_cr(pos_cr, EXPR, TXDB, gap = eval(parse(text = gap2)), mc.cores = 4)

cr_reduced = c(neg_cr_reduced, pos_cr_reduced)
cr_reduced = cr_reduced[order(as.vector(seqnames(cr_reduced)), cr_reduced$gene_id, start(cr_reduced))]

if(length(unique(attr(cr_filtered, "factor"))) > 1) {
	cr_reduced = add_subtype_specificity(cr_reduced)
}

saveRDS(cr_reduced, file = qq("@{OUTPUT_DIR}/rds/cr_reduced_fdr_@{cutoff}.rds"))
