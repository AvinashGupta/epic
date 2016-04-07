
suppressPackageStartupMessages(library(GetoptLong))

cutoff = 0.01
gtf = NULL
peak = NULL
GetoptLong("config=s", "configuration R script",
	       "cutoff=f", "cutoff for filter cr",
	       "chr=s", "chromosome",
	       "gtf=s", "gtf file, it is used to extract gene symbol",
	       "peak=s", "peak name")

library(epic)
load_config(config)


files = dir(qq("@{OUTPUT_DIR}/rds"), pattern = "^cr_filtered_fdr_.*\.rds$")
if(length(files) == 1) {
	cutoff = gsub("^cr_filtered_fdr_(.*)\.rds$", "\\1", files[1])
	cr_filtered = readRDS(files[1])
} else {
	cr_filtered = readRDS(qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))
}

if(!chr %in% unique(as.character(seqnames(cr_filtered)))) {
	stop("'chr' does not exist in 'cr_filtered'.")
}

if(is.null(peak)) {
	peak_list = NULL
} else {
	peak_list = get_peak_list(peak)
}

tx_list = transcriptsBy(TXDB, by = "gene")
cr_subset = cr_filtered[seqnames(cr_filtered) == chr]
if(is.null(gtf)) {
	for(gi in unique(cr_subset$gene_id)) {
		pdf(qq("@{OUTPUT}/gviz/gviz_@{chr}_@{gi}.pdf"), width = 16, height = 12)
	    cr_gviz(cr_filtered, gi, expr, txdb, tx_list = tx_list[[gi]]$tx_name, 
	    	gf_list = GENOMIC_FEATURE_LIST[intersect(c("cgi", "tfbs", "enhancer"), names(GENOMIC_FEATURE_LIST))], 
	    	hm_list = peak_list)
	    dev.off()
	}
} else {
	gn = extract_field_from_gencode(gtf, level = "gene", primary_key = "gene_id", field = "gene_name")
	for(gi in unique(cr_subset$gene_id)) {
	    pdf(qq("@{OUTPUT}/gviz/gviz_@{chr}_@{gi}_@{gn[gi]}.pdf"), width = 16, height = 12)
	    cr_gviz(cr_filtered, gi, expr, txdb, tx_list = tx_list[[gi]]$tx_name, 
	    	gf_list = GENOMIC_FEATURE_LIST[intersect(c("cgi", "tfbs", "enhancer"), names(GENOMIC_FEATURE_LIST))], 
	    	hm_list = peak_list, symbol = gn[gi])
	    dev.off()
	}
}

