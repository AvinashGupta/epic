
suppressPackageStartupMessages(library(GetoptLong))

cutoff = NULL
GetoptLong("config=s", "configuration R script",
	       "cutoff=f", "cutoff for filter cr")

library(epic)
load_config(config)

if(is.null(cutoff)) {
	cutoff = CR_CUTOFF
}

cr_filtered = filter_correlated_regions(chromosome = CHROMOSOME, template = qq("`OUTPUT_DIR`/rds/@{chr}_cr.rds", code.pattern = "`CODE`"), cutoff = cutoff)

cv = rowIQRs(expr)/rowMedians(expr)
names(cv) = rownames(expr)
cr_filtered$expr_cv = cv[cr_filtered$gene_id]

if(length(unique(attr(cr_filtered, "factor"))) > 1) {
	cr_filtered = add_subtype_specificity(cr_filtered)
}

saveRDS(cr_filtered, file = qq("@{OUTPUT_FOLDER}/rds/cr_filtered_fdr_@{cutoff}.rds"))

