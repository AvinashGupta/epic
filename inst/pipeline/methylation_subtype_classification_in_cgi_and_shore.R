
# this script depends on `differential_methylation_in_cgi_and_shore.R`

suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script",
           "n_class=i", "number of classes expected")

library(epic)
load_config(config)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

cat("subtype classification based on cgi\n")
gr_cgi = readRDS(qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi.rds}"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_cgi, n_class = n_class, pct_cutoff = 0.2, corr_cutoff = 0.6, k = 500, ha = ha)
dev.off()


cat("subtype classification based on cgi shores\n")
gr_shore = readRDS(qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi_shore.rds}"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_shore, n_class = n_class, pct_cutoff = 0.2, corr_cutoff = 0.6, k = 500, ha = ha)
dev.of()


cat("subtype classification based on other parts in genome\n")
gr_complement = readRDS(qq("@{OUTPUT_DIR/rds/mean_meth_1kb_neither_cgi_nor_shore.rds}"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_complement, n_class = n_class, pct_cutoff = 0.02, corr_cutoff = 0.8, k = 1000, ha = ha)
dev.of()
