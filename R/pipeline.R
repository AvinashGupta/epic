 

# == title
# run pipeline through qsub system
#
# == param
# -config_file path of configuration script
# -prefix prefix of the job name
# -email email
# -enforce enforce run all the steps
# -Rscript_binary path of Rscript binary
#
# == details
# Automatically run 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
pipeline_by_qsub = function(config_file, prefix = "", email = NULL, enforce = FALSE, Rscript_binary = "Rscript") {

	OUTPUT_DIR = NULL
	SAMPLE = NULL
	EXPR = NULL
	TXDB = NULL
	CHROMOSOME = NULL
	MARKS = NULL

	load_config(config_file, validate = FALSE)

	tmpdir = paste0(OUTPUT_DIR, "/temp")

	if(length(unique(SAMPLE$class)) > 1) {
		x = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' differential_methylation_in_cgi_and_shore --config '@{config_file}'"), 
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
			enforce = enforce,
			tmpdir = tmpdir)

		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' differential_methylation_in_genomic_features --config '@{config_file}'"),
			output = c(qq("@{OUTPUT_DIR}/heatmap_diff_methylation_in_genomic_features.pdf")),
			name = qq("@{prefix}differential_methylation_in_genomic_features"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir)

		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore --config '@{config_file}'"),
			output = c(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf")),
			name = qq("@{prefix}methylation_subtype_classification_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = x,
			tmpdir = tmpdir)
	} else {
		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore --config '@{config_file}'"),
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
			enforce = enforce,
			tmpdir = tmpdir)
	}

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' general_methylation_distribution --config '@{config_file}'"),
		output = c(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf")),
		name = qq("@{prefix}general_methylation_distribution"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		tmpdir = tmpdir)

	#####################################
	### CR
	#####################################

	if(is.null(EXPR) || is.null(TXDB)) {
		return(NULL)
	}

	dependency = NULL
	for(chr in CHROMOSOME) {
		x = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions --config '@{config_file}' --chr @{chr}"),
			output = qq("@{OUTPUT_DIR}/rds/@{chr}_cr.rds"),
			name = qq("@{prefix}correlated_regions_@{chr}"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir)
		dependency = c(dependency, x)
	}
	dependency = paste(dependency, collapse = ",")

	pid_cr_filter = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_filter --config '@{config_file}'"),
		output = qq("@{OUTPUT_FOLDER}/rds/cr_filtered_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_filter"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = dependency,
		tmpdir = tmpdir)

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_reduce --config '@{config_file}'"),
		output = qq("@{OUTPUT_DIR}/rds/cr_reduced_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_reduce"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = pid_cr_filter,
		tmpdir = tmpdir)

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_downstream --config '@{config_file}'"),
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
		dependency = pid_cr_filter,
		tmpdir = tmpdir)

	for(chr in CHROMOSOME) {
		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_gviz --config '@{config_file}' --chr @{chr}"),
			output = NULL,
			name = qq("@{prefix}correlated_regions_gviz"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = pid_cr_filter,
			tmpdir = tmpdir)
	}

	if(is.null(MARKS)) {
		return(NULL)
	}

	for(mk in MARKS) {
		for(which in c("pos", "neg")) {
			pipeline_step("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_enriched --config '@{config_file}' --peak @{mk} --which @{which}",
				output = NULL,
				name = qq("@{prefix}correlated_regions_enriched_@{mk}_@{which}"),
				walltime = "10:00:00",
				mem = "10G",
				nodes = 1,
				email = email,
				enforce = enforce,
				dependency = pid_cr_filter,
				tmpdir = tmpdir)
		}
	}
}


pipeline_step = function(..., output = NULL, name, walltime = "1:00:00", mem = "1G", nodes = 1, 
	email = NULL, dependency = NULL, enforce = FALSE, tmpdir) {
	
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
	temp_file = tempfile(tmpdir = tmpdir, fileext = ".sh")
	writeLines(script, temp_file)

	x = system(qq("qsub @{temp_file}"), intern = TRUE)
	return(x[1])
}


# x = pipeline_step("ls", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de")
# x = pipeline_step("ls -l", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de", dependency = x)

