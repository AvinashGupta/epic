
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script")

library(epic)
load_config(config)

sample_id = rownames(SAMPLE)
n_sample = length(sample_id)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

gr_list = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = GENOMIC_FEATURE_LIST)
gr_name = names(gr_list)

pdf(qq("@{OUTPUT_DIR}/heatmap_diff_methylation_in_genomic_features.pdf"), width = 16, height = 16)
for(i in seq_along(gr_list)) {
	qqcat("making heatmap for @{gr_name[i]}\n")
	heatmap_diff_methylation_in_genomic_features(gr_list[[i]], annotation = SAMPLE$type, 
		annotation_color = COLOR$class, title = gr_name[i], txdb = TXDB, gf_list = GENOMIC_FEATURE_LIST, ha = ha)
}
dev.off()
