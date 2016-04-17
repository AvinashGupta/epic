
suppressPackageStartupMessages(library(GetoptLong))
cutoff = 0.01
GetoptLong("config=s", "configuration R script",
		   "cutoff=f", "cutoff for filter cr",
	       "chr=s", "chromosomes")

library(epic)
load_config(config)



files = dir(qq("@{OUTPUT_DIR}/rds"), pattern = "^cr_filtered_fdr_.*\\.rds$")
if(length(files) == 1) {
	cutoff = gsub("^cr_filtered_fdr_(.*)\\.rds$", "\\1", files[1])
	cr_filtered = readRDS(files[1])
} else {
	cr_filtered = readRDS(qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))
}


neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]

n1 = length(neg_cr) * 5
n2 = length(pos_cr) * 5
w1 = sum(width(neg_cr))
w2 = sum(width(pos_cr))

pdf(qq("@{OUTPUT_DIR}/cr_number.pdf"), width = 4, height = 6)
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
barplot(c("neg_cr" = n1, "pos_cr" = n2), beside = TRUE, col = c("green", "red"), 
    ylab = "#CpG", axes = FALSE)
axis(side = 2, at = c(0, 2, 4, 6, 8)*100000, labels = c("0K", "200K", "400K", "600K", "800K"))
barplot(c("neg_cr" = w1, "pos_cr" = w2), beside = TRUE, col = c("green", "red"), 
    ylab = "sum(width(cr))", axes = FALSE)
axis(side = 2, at = c(0, 10, 20, 30)*1000000, labels = c("0MB", "10MB", "20MB", "30MB"))
dev.off()


pdf(qq("@{OUTPUT_DIR}/hilbert_sig.pdf"), width = 8, height = 8)
cr_hilbert(cr = cr_filtered, txdb = TXDB, merge_chr = TRUE)
cr_hilbert(cr = cr_filtered, txdb = TXDB, merge_chr = FALSE)
dev.off()


pdf(qq("@{OUTPUT_DIR}/hilbert_all.pdf"), width = 18, height = 12)
cr_hilbert(template = paste0(OUTPUT_DIR, "/rds/@{chr}_cr.rds"), txdb = TXDB, merge_chr = TRUE)
cr_hilbert(template = paste0(OUTPUT_DIR, "/rds/@{chr}_cr.rds"), txdb = TXDB, merge_chr = FALSE)
dev.off()


pdf(qq("@{OUTPUT_DIR}/cr_overlap.pdf"), width = 8, height = 5)
cr_overlap_to_genomic_features(cr_filtered, GENOMIC_FEATURE_LIST, chromosome = CHROMOSOME)
dev.off()


## how cr enriched at tss
pdf(qq("@{OUTPUT_DIR}/cr_tss.pdf"), width = 10, height = 4)
cr_enriched_at_tss(cr_filtered, TXDB)
dev.off()
