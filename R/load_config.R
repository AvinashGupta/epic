

# == title
# Load and validate configuration file
#
# == param
# -config_file path of configuration file
# -export_env environment where to export variables
# -validate whether do validation
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
load_config = function(config_file, export_env = parent.frame(), validate = TRUE) {
	cat(
"Configuration file should provide following variables:

SAMPLE: a data frame, row names must be sample id and there must be a
  'class' column which corresponds to classes of samples. There can also
  be other annotation columns.

COLOR: a list of color settings corresponding to annotation column in 
  'SAMPLE'. The value of each list element must be either a named vector
  or a color mapping function. 'COLOR$class' must be defined or random color
  will be assigned. Names of other color settings should be same as
  corresponding columns in 'SAMPLE'.

TXDB (optional): a 'TxDB' object.

EXPR (optional): a matrix which contains expression values. Column names 
  should be sample id and row names should be gene ids. Note gene ids in the 
  expression matrix should be same type as genes in 'TXDB'.

CHROMOSOME: a vector of chromosome names.

GENOME: abbreviation of species.

OUTPUT_DIR: path of output directory. Several sub directories will be created.

GENOMIC_FEATURE_LIST: a list of genomic features as GRanges objects. There
  must be a element named 'cgi'.

MARKS (optional): a vector of ChIP-Seq markers.

methylation_hooks() must be defined.
chipseq_hooks() is optional unless you want to do integrative analysis.

")
	SAMPLE = NULL
	COLOR = NULL
	TXDB = NULL
	EXPR = NULL
	CHROMOSOME = NULL
	GENOME = NULL
	OUTPUT_DIR = NULL
	GENOMIC_FEATURE_LIST = NULL
	MARKS = NULL

	sys.source(config_file, envir = environment())

	if(!validate) {
		assign("SAMPLE", SAMPLE, envir = export_env)
		assign("COLOR", COLOR, envir = export_env)
		assign("TXDB", TXDB, envir = export_env)
		assign("EXPR", EXPR, envir = export_env)
		assign("CHROMOSOME", CHROMOSOME, envir = export_env)
		assign("GENOME", GENOME, envir = export_env)
		assign("OUTPUT_DIR", OUTPUT_DIR, envir = export_env)
		assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
		assign("MARKS", MARKS, envir = export_env)
	}

	# test SAMPLE
	if(is.null(SAMPLE$class)) {
		SAMPLE$class = rep("sample", nrow(SAMPLE))
		warning("There should be a 'class' column in 'SAMPLE'.")
	}
	if(is.null(rownames(SAMPLE))) {
		stop("'SAMPLE must have row names.")
	}

	if(is.null(COLOR$class)) {
		class = unique(COLOR$class)
		COLOR$class = structure(rand_color(length(class)), names = class)
		warning("'COLOR$class' should be defined")
	}

	cn = intersect(names(COLOR), colnames(SAMPLE))
	COLOR = COLOR[cn]
	SAMPLE = SAMPLE[, cn, drop = FALSE]

	sample_id = rownames(SAMPLE)

	# test methylation_hooks()
	cat("Randomly pick one chromosome.\n")
	chr = sample(CHROMOSOME, 1)
	methylation_hooks$set(chr)
	meth = methylation_hooks$meth(row_index = 1:5)
	cov = methylation_hooks$coverage(row_index = 1:5)
	methylation_hooks$site()[1:2]
	methylation_hooks$GRanges()[1:2]

	if(length(intersect(sample_id, colnames(meth))) == 0) {
		stop("Cannot match column names in methylation data to sample ids in 'SAMPLE'.")
	}
	if(length(intersect(sample_id, colnames(cov))) == 0) {
		stop("Cannot match column names in coverage data to sample ids in 'SAMPLE'.")
	}

	if(is.null(EXPR) || is.null(TXDB)) {
		# test txdb and expr
		if(length(intersect(sample_id, colnames(EXPR))) == 0) {
			stop("Cannot match column names in 'EXPR' to sample ids in 'SAMPLE'.")
		}
		if(length(intersect(colnames(meth), colnames(EXPR))) == 0) {
			stop("Cannot match column names in 'EXPR' to column names in methylation data.")
		}

		genes = genes(TXDB)
		if(length(intersect(names(genes), rownames(EXPR))) == 0) {
			stop("Cannot match genes in 'EXPR' to 'TXDB'.")
		}

		# check chromosome and species
		if(length(intersect(CHROMOSOME, getChromInfoFromUCSC(GENOME)[, 1])) == 0) {
			stop("Cannot match 'GENOME' to 'CHROMOSOME'.")
		}

		chr = as.vector(seqnames(genes))
		names(chr) = names(genes)
		EXPR = EXPR[chr[rownames(EXPR)] %in% CHROMOSOME, , drop = FALSE]
	}

	# initialize the folder structure
	dir.create(OUTPUT_DIR, showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/gviz"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/rds"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/temp"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/enriched_heatmap"), showWarnings = FALSE)

	# genomic features
	if(is.null(GENOMIC_FEATURE_LIST$cgi)) {
		stop("'GENOMIC_FEATURE_LIST' must have an element 'cgi'.")
	}
	gf_name = names(GENOMIC_FEATURE_LIST)
	for(i in seq_along(GENOMIC_FEATURE_LIST)) {
		if(length(setdiff(as.character(seqnames(GENOMIC_FEATURE_LIST[[i]])), CHROMOSOME))) {
			stop(qq("GENOMIC_FEATURE_LIST$@{gf_name[i]} should only contain chromosomes in 'CHROMOSOME'"))
		}
	}

	# validate chipseq data input
	for(mk in MARKS) {
		sample_id = chipseq_hooks$sample_id(mk)
		if(length(sample_id) == 0) {
			stop(qq("No ChIP-Seq sample detected for mark @{mk}."))
		}
		if(length(intersect(rownames(SAMPLE), sample_id)) == 0) {
			stop(qq("No ChIP-Seq sample in 'SAMPLE' for mark @{mark}."))
		}
		if(length(setdiff(sample_id, rownames(SAMPLE)))) {
			stop(qq("There are ChIP-Seq samples which are not in 'SAMPLE'. Please modify the 'sample_id' hook."))
		}
		sid = sample(sample_id, 1)
		qqcat("random pick one sample: @{sid}\n")
		peak_list = get_peak_list(mk, sid)
		for(i in seq_along(peak_list)) {
			if(length(setdiff(as.character(seqnames(peak_list[[i]])), CHROMOSOME))) {
				stop(qq("peaks should only contain chromosomes in 'CHROMOSOME', mark @{mk}. Please modify the 'peak' hook."))
			}
			if(!"density" %in% colnames(mcols(peak_list[[i]]))) {
				stop(qq("Signal of peaks should be in the 'density' column. Please modify the 'peak' hook."))
			}
		}
	}


	assign("SAMPLE", SAMPLE, envir = export_env)
	assign("COLOR", COLOR, envir = export_env)
	assign("TXDB", TXDB, envir = export_env)
	assign("EXPR", EXPR, envir = export_env)
	assign("CHROMOSOME", CHROMOSOME, envir = export_env)
	assign("GENOME", GENOME, envir = export_env)
	assign("OUTPUT_DIR", OUTPUT_DIR, envir = export_env)
	assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
	assign("MARKS", MARKS, envir = export_env)

	cat("\nvalidation passed and global variables are imported.\n")
}
