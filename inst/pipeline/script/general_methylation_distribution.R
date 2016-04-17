
suppressPackageStartupMessages(library(GetoptLong))
shore_extend = 2000
GetoptLong("config=s", "configuration R script",
	       "shore_extend=i", "base pairs that extend from CGI")

library(epic)
load_config(config)

sample_id = rownames(SAMPLE)
n_sample = length(sample_id)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

cat("general methylation distribution whole genome-wide\n")
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"), width = 0.16*n_sample + 2, height = 10)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, ha = ha)
dev.off()

cat("general methylation distribution in cpg islands\n")
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"), width = 0.16*n_sample + 2, height = 10)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01, ha = ha)
dev.off()


cat("general methylation distribution in cgi shores\n")
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"), width = 0.16*n_sample + 2, height = 10)
extended_cgi = GENOMIC_FEATURE_LIST$cgi
start(extended_cgi) = start(extended_cgi) - shore_extend
end(extended_cgi) = end(extended_cgi) + shore_extend
shore = setdiff(extended_cgi, GENOMIC_FEATURE_LIST$cgi)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = shore, p = 0.01, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = shore, p = 0.01, ha = ha)
dev.off()

cat("general methylation distribution in neither cgi nor shores\n")
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf"), width = 0.16*n_sample + 2, height = 10)
chromInfo = getChromInfoFromUCSC(GENOME)
chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
complement = setdiff(chromGr, extended_cgi)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = complement, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = complement, ha = ha)
dev.off()
