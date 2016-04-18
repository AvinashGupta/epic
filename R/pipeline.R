 


pipeline_by_qsub = function(config_file, prefix = "", email = NULL, enforce = FALSE) {

	load_config(config_file, validate = FALSE)

	if(length(unique(SAMPLE$class)) > 1) {
		x = pipeline_step("Rscript -e 'epic::epic()' differential_methylation_in_cgi_and_shore", 
			output = c(qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi.rds}"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi_shore.rds}"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_neither_cgi_nor_shore.rds}"),
				       qq("@{OUTPUT_DIR}/heatmap_diff_methylation_1kb_window.pdf"),
				       qq("@{OUTPUT}/genome_diff_1kb_window_correlation.pdf")),
			name = qq("@{prefix}differential_methylation_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce)

		pipeline_step("Rscript -e 'epic::epic()' differential_methylation_in_genomic_features",
			output = c(qq("@{OUTPUT_DIR}/heatmap_diff_methylation_in_genomic_features.pdf")),
			name = qq("@{prefix}differential_methylation_in_genomic_features"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce)

		pipeline_step("Rscript -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore",
			output = c(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf")),
			name = qq("@{prefix}methylation_subtype_classification_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = x)
	} else {
		pipeline_step("Rscript -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore",
			output = c(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi.rds}"), 
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi_shore.rds}"), 
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_neither_cgi_nor_shore.rds}")),
			name = qq("@{prefix}methylation_subtype_classification_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce)
	}

	pipeline_step("Rscript -e 'epic::epic()' general_methylation_distribution",
		output = c(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf")),
		name = qq("@{prefix}general_methylation_distribution"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce)

	#####################################
	### CR
	#####################################

	if(is.null(EXPR) || is.null(TXDB)) {
		return(NULL)
	}

	dependency = NULL
	for(chr in CHROMOSOME) {
		x = pipeline_step(qq("Rscript -e 'epic::epic()' correlated_regions --chr @{chr}"),
			output = qq("@{OUTPUT_DIR}/rds/@{chr}_cr.rds"),
			name = qq("@{prefix}correlated_regions_@{chr}"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce)
		dependency = c(dependency, x)
	}
	dependency = paste(dependency, collapse = ",")

	x = pipeline_step("Rscript -e 'epic::epic()' correlated_regions_filter",
		output = qq("@{OUTPUT_FOLDER}/rds/cr_filtered_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_filter"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = dependency)

	pipeline_step("Rscript -e 'epic::epic()' correlated_regions_reduce",
		output = qq("@{OUTPUT_DIR}/rds/cr_reduced_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_reduce"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = x)

	pipeline_step("Rscript -e 'epic::epic()' correlated_regions_downstream",
		output = c(qq("@{OUTPUT_DIR}/cr_number.pdf"),
			       qq("@{OUTPUT_DIR}/hilbert_sig.pdf"),
			       qq("@{OUTPUT_DIR}/hilbert_all.pdf"),
			       qq("@{OUTPUT_DIR}/cr_overlap.pdf"),
			       qq("@{OUTPUT_DIR}/cr_tss.pdf")),
		name = qq("@{prefix}correlated_regions_downstream"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = x)

	for(chr in CHROMOSOME) {
		pipeline_step(qq("Rscript -e 'epic::epic()' correlated_regions_gviz --chr @{chr}"),
			output = NULL,
			name = qq("@{prefix}correlated_regions_gviz"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = x)
	}

	if(is.null(MARKS)) {
		return(NULL)
	}

	for(mk in MARKS) {
		for(which in c("pos", "neg")) {
			pipeline_step("Rscript -e 'epic::epic()' correlated_regions_enriched --peak @{mk} --which @{which}",
				output = NULL,
				name = qq("@{prefix}correlated_regions_enriched_@{mk}_@{which}"),
				walltime = "10:00:00",
				mem = "10G",
				nodes = 1,
				email = email,
				enforce = enforce,
				dependency = x)
		}
	}
}


pipeline_step = function(..., output = NULL, name, walltime = "1:00:00", mem = "1G", nodes = 1, 
	email = NULL, dependency = NULL, enforce = FALSE) {
	
	if(is.null(output)) {
		enforce = TRUE
	}
	if(!enforce && all(file.exists(output))) {
		return(NULL)
	}
	cmd = unlist(list(...))
	cmd = paste(cmd, collapse = "\n")
	script = qq("#!/bin/sh
#PBS -j oe
#PBS -o @{OUTPUT_DIR}/temp/
@{ifelse(is.null(dependency), '', paste0('#PBS -W depend=afterok:', dependency))}
#PBS -N @{name}
#PBS -l walltime=@{walltime}
#PBS -l mem=@{mem}
#PBS -l nodes=@{nodes}
@{ifelse(is.null(email), '', paste0('#PBS -M', email))}

@{cmd}
")
	temp_file = tempfile(tmpdir = OUTPUT_DIR, fileext = ".sh")
	writeLines(script, temp_file)

	con = pipe(qq("qsub @{temp_file}"))
	x = scan(con, what = "character")
	return(x[1])
}


# x = pipeline_step("ls", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de")
# x = pipeline_step("ls -l", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de", dependency = x)

