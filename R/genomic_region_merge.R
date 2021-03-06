
# == title
# Segmentation by a logical vector
#
# == param
# -l a logical vector
#
# == details
# The logical vector will be segmented according to their values.
# It returns intervals for continuous `TRUE` values.
#
# == values
# A data frame in which the first column is the index of start sites in the original vector and
# the second column is the index of end sites.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# l = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
# logical_segment(l)
logical_segment = function(l) {
	 w = which(l)

	 # no cluster
	 if(! length(w)) {
	 	return(data.frame(start_index = numeric(0), end_index = numeric(0)))
	 }

	 d = diff(w) > 1

	 # only one cluster
	 if(! any(d)) {
	 	return(data.frame(start_index = w[1], end_index = w[length(w)]))
	 } 

	 # more than one clusters
	 # second step: dynamic merge
	 #   if distance between two clusters are smaller than short cluster, merge the two clusters
	 di = which(d)
	 end = w[di]
	 start = w[di+1]
	 start = c(w[1], start)
	 end = c(end, w[length(w)])

	return(data.frame(start_index = start, end_index = end))
}

# == title
# Mark the numbers to be number of base pairs
#
# == param
# -x a numeric vector. It will be convert to integers by `base::as.integer`.
#
# == details
# It adds a new ``bp`` class to the vector.
#
# == value
# A same vector as input
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# bp(10)
# bp(10.1)
bp = function(x) {
	x = as.integer(x)
	class(x) = c(class(x), "bp")
	x
}

# == title
# Mark that the numbers represent number of kilo bases
#
# == param
# -x a numeric vector.
#
# == details
# The input values are multiplied by 1000 and send to `bp`.
#
# == value
# A numeric vector measured in bp
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# kb(10)
# kb(10.01)
kb = function(x) {
	bp(x*1000)
}

# == title
# Mark that the numbers represent number of mega bases
#
# == param
# -x a numeric vector.
#
# == details
# The input values are multiplied by 1000000 and send to `bp`.
#
# == value
# A numeric vector measured in bp
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# mb(10)
# mb(10.01)
mb = function(x) {
	bp(x*1000000)
}

# == title
# Print bp class objects
#
# == param
# -x `bp` class object
# -... other arguments
# 
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
print.bp = function(x, ...) {
	x = paste0(x, " bp")
	print(x, quote = FALSE)
}


# == title
# Merge genomic regions
#
# == param
# -gr a `GenomicRanges::GRanges` object
# -max_gap maximum gap to merge, measured in base pairs. Only work if ``gap`` is the ratio.
# -gap a numeric value means to extend each region by ``gap`` times of its width before merging. If ``gap`` represents
#      number of base pairs, use `bp`, `kb` or `mb` to wrap it. If ``gap`` represents absolute number of base pairs,
#      the functionality is same as `GenomicRanges::reduce`.
# -add_n_column whether to add a column which represents number of regions merged, used internally.
#
# == details
# `GenomicRanges::reduce` only merges regions with fixed gap width, but sometimes it is not reasonable to set gap
# to a same width for all regions. Assuming we have a list of differentially methylated regions (DMRs) and we want to reduce
# the number of DMRs by merging neighouring DMRs. DMRs distribute differently in different places in the genome, e.g. DMRs are dense
# and short in CpG-rich regions (e.g. CpG islands) while long in CpG-sparse regions (e.g. gene bodies and intergenic regions),
# thus the merging should be applied based to the width of every DMR itself. `reduce2` can merge regions by the width of every region itself.
# This type of merging is dynamic because after each iteration of merging, width and number of regions change which results in changing of extension
# of some regions. The whole merging will proceed iteratively unless there is no new merging or the gap between regions reaches ``max_gap``.
#
# If there are numeric meta columns, corresponding values will be summed up for the merged regions. There will
# be a new column ``.__n__.`` added which represents number of regions that are merged.
#
# == value
# a `GenomicRanges::GRanges` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,12), 
#     end = c(2, 5, 10, 20)))
# reduce2(gr, gap = bp(2))
# reduce2(gr, gap = 0.1)
#
reduce2 = function(gr, max_gap = 1000, gap = bp(1000), add_n_column = TRUE) {
	if(length(gr) %in% c(0, 1)) {
		return(gr)
	}

	if(add_n_column) gr$.__n__. = rep(1, length(gr))

	if(inherits(gap, "bp")) {
		gr_reduced = reduce(gr, min.gapwidth = gap, with.revmap = TRUE)
		revmap = mcols(gr_reduced)[, 1]
		if(is.null(mcols(gr))) {
			mcols(gr_reduced) = NULL	
		} else {
			mt = mcols(gr)
			l = sapply(mt, class) == "numeric"
			if(sum(l) == 0) {
				mcols(gr_reduced) = NULL
			} else {
				mat = as.matrix(mt[, l, drop = FALSE])
				cn = colnames(mat)
				mat = do.call("rbind", lapply(revmap, function(ind) colSums(mat[ind, , drop = FALSE])))
				colnames(mat) = cn
				mcols(gr_reduced) = mat
			}
		}
		return(gr_reduced)
	} else {

		if(round(gap) == gap && gap > 50) {
			qqcat("Extend each region by @{gap} times of its width. If the gap value represents absolute base pairs, use `bp()` to wrap it.")
		}
		n = length(gr)
		gr_extend = gr
		gr_width = width(gr)
		start(gr_extend) = start(gr) - round(gr_width*gap)
		end(gr_extend) = end(gr) + round(gr_width*gap)

		gr_reduced = reduce(gr_extend, with.revmap = TRUE)
		revmap = mcols(gr_reduced)[, 1]
		s = start(gr)
		e = end(gr)
		os = sapply(revmap, function(ind) min(s[ind]))
		oe = sapply(revmap, function(ind) max(e[ind]))
		start(gr_reduced) = os
		end(gr_reduced) = oe

		if(is.null(mcols(gr))) {
			mcols(gr_reduced) = NULL	
		} else {
			mt = mcols(gr)
			l = sapply(mt, class) == "numeric"
			if(sum(l) == 0) {
				mcols(gr_reduced) = NULL
			} else {
				mat = as.matrix(mt[, l, drop = FALSE])
				cn = colnames(mat)
				mat = do.call("rbind", lapply(revmap, function(ind) colSums(mat[ind, , drop = FALSE])))
				colnames(mat) = cn
				mcols(gr_reduced) = mat
			}
		}

		n2 = length(gr_reduced)

		if(n2 == 1) {
			return(gr_reduced)
		} else if(n == n2) {
			return(gr_reduced)
		} else {
			qqcat("  regions have been reduced from @{n} to @{n2}...\n")
			gr = reduce2(gr_reduced, max_gap = max_gap, gap = gap, add_n_column = FALSE)
			return(gr)
		}
	}
}
