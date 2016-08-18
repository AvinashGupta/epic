

# == title
# Hook functions to extract methylation data
#
# == param
# -... Arguments for the parameters, see "details" section
# -RESET reset to default values
# -READ.ONLY whether only return read-only options
# -LOCAL switch local mode
#
# == detail
# Methylation from whole genome bisulfite sequencing is always huge and it does not
# make sense to read them all into the memory (imaging there are 20M CpG sites on single strand in human genome). 
# This hook sets how to read the methylation
# data and how to return methylation data (e.g. CpG coverage, methylation rate...).
#
# All downstream functions which analyze methylation data needs this hook to be already set.
# 
# There are following hooks:
#
# -get_data how to get the object which contains methylation data. The function accepts a single chromosome name and 
#      returns an object which is used as the first argument in other hook functions
# -meth how to extract methylation rate. The function should have three arguments:
#       the object returned from ``set()``, index of rows and index of columns. Normally,
#       the first argument (``obj``) can be ignored when calling this hook. Note the methylation
#       matrix should have column names. The methylation rate should between 0 and 1.
# -raw how to extract raw methylation value, same setting as ``meth``. This hook is optional.
# -site the function should return a vector of positions of CpG sites
# -coverage how to extract CpG coverage, same setting as ``meth``.
# -GRanges howt to extract CpG sites as a `GenomicRanges::GRanges` object.
#
# Following two hooks can be used if above hooks are set:
#
# -set select chromosome as the current chromosome
# -sample_id a vector of sample ids which contains methylation data
#
# Note: positions of CpG sites in a chromosome should be sorted.
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
methylation_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
methylation_hooks = setGlobalOptions(
	set = list(.value = NULL, .class = "function"),
	get_data = list(.value = function(chr) stop("you need to define `get_data`"),
	                             .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	meth  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `meth` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	raw  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `raw` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	site         = list(.value = function(obj)  stop("you need to define `site` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 3
	                             	}),
	coverage     = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `coverage` hook"),
		                       .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	GRanges      = list(.value =function(obj)  stop("you need to define `GRanges` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	obj = NULL,
	sample_id = list(.value = NULL, .class = "character")
)

methylation_hooks$set = function(chr) {

	if(!is.null(attr(methylation_hooks$obj, "chr"))) {
		if(attr(methylation_hooks$obj, "chr") == chr) {
	        qqcat("[@{chr}] @{chr} is already set.\n")
	        return(invisible(NULL))
	    }
	}
    
    obj = methylation_hooks$get_data(chr)
    attr(obj, "chr") = chr

    methylation_hooks$obj = obj

    methylation_hooks$sample_id = colnames(methylation_hooks$meth())

    return(invisible(NULL))
}

# .obj_is_set = function() {
# 	!is.null(methylation_hooks$obj)
# }

# == title
# Hook functions to extract ChIP-Seq peak data
#
# == param
# -... Arguments for the parameters, see "details" section
# -RESET reset to default values
# -READ.ONLY whether only return read-only options
# -LOCAL switch local mode
#
# == details
# This hook defines how to get sample ids for a specific marks and how to get peak regions
# by given mark type and sample id:
#
# -sample_id how to extract sample ids. The argument is the name of the mark and it returns a vector of sample ids.
# -peak how to get peak regions. The argument is the name of the mark and a single sample id. The function returns a `GenomicRanges::GRanges` object.
# -chromHMM how to get chromHMM data. The argument is a single sample id and the function returns a `GenomicRanges::GRanges` object. Note
#           the first column in meta columns should be the chromatin states.
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
chipseq_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
chipseq_hooks = setGlobalOptions(
	sample_id = list(.value = function(mark) stop("you need to define `sample_id` hook"),
		             .class = "function",
		             .validate = function(f) length(as.list(f)) == 2),
	peak = list(.value = function(mark, sid) stop("you need to define `peak` hook"),
		        .class = "function",
		        .validate = function(f) length(as.list(f)) == 3),
	chromHMM = list(.value = function(sid) stop("you need to define `chromHMM` hook"),
		            .class = "function",
		            .validate = function(f) length(as.list(f)) == 2)
)

# == title
# Get a list of peak regions
#
# == param
# -mark mark type
# -sample_id a vector of sample ids
#
# == details
# It works after `chipseq_hooks` is set.
#
# == value
# A list of `GenomicRanges::GRanges` objects.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_peak_list = function(mark, sample_id = chipseq_hooks$sample_id(mark)) {
    peak_list = lapply(sample_id, function(sid) chipseq_hooks$peak(mark, sid))
    names(peak_list) = sample_id
    peak_list
}
