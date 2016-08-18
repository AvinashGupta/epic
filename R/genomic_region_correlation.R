###############################################################
# correlation between two sets of genomic regions
# or enrichment of one type of genomic region in the other
###############################################################


# == title
# Correlation between two sets of genomic regions
#
# == param
# -gr_list_1 a list of `GenomicRanges::GRanges` objects, should be a named list, e.g. low methylated regions in different samples.
# -gr_list_2 a list of `GenomicRanges::GRanges` objects, should be a named list, e.g. a list of different genomic features.
# -background a `GenomicRanges::GRanges` object. The correlation is only looked in background regions.
# -chromosome a vector of chromosome names
# -species species, used for random shuffling genomic regions
# -nperm number of random shufflings. If it is set to 0, no random shuffling will be performed.
# -mc.cores number of cores for parallel calculation
# -stat_fun method to calculate correlations. There are some pre-defined functions:
#           `genomic_corr_reldist`, `genomic_corr_absdist` measure how two sets of genomic regions are close; `genomic_corr_jaccard`,
#           `genomic_corr_nintersect`, `genomic_corr_pintersect`, `genomic_corr_sintersect` measures how two sets of genomic regions are overlapped.
#           The self-defined function should accept at least two arguments which are two GRanges object.
#           The third argument is ``...`` which is passed from the main function. The function
#           should only return a numeric value.
# -... pass to ``stat_fun``
# -bedtools_binary random shuffling is perfomed by ``bedtools``. If ``bedtools`` is not in ``PATH``, the path of ``bedtools`` can be set here.
# -tmpdir dir for temporary files
#
# == details
# The correlation between two sets of genomic regions basically means how much the first type of genomic regions
# are overlapped or close to the second type of genomic regions.
#
# The significance of the correlation is calculated by random shuffling the regions. 
# In random shuffling, regions in ``gr_list_1`` will be shuffled. So if you want to shuffle ``gr_list_2``,
# just switch the first two arguments.
#
# Pleast note random shuffling is done by bedtools, so bedtools should be installed and exists in ``PATH``
# and should support ``-i -g -incl`` options.
#
# == value
# A list containing following elements:
#
# -stat statistic value
# -fold_change stat/E(stat), stat divided by expected value which is generated from random shuffling
# -p.value p-value for over correlated. So, 1 - p.value is the significance for being less correlated
# -stat_random_mean mean value of stat in random shuffling
# -stat_random_sd standard deviation in random shuffling
#
# If ``perm`` is set to 0 or 1, ``fold_change``, ``p.value``, ``stat_random_mean`` and ``stat_random_sd`` are all ``NULL``.
#
# == seealso
# `genomic_corr_reldist`, `genomic_corr_jaccard`, `genomic_corr_absdist`, `genomic_corr_nintersect`, 
# `genomic_corr_pintersect`, `genomic_corr_sintersect`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
# gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
# genomic_regions_correlation(gr1, gr2, nperm = 0)
# genomic_regions_correlation(list(gr1 = gr1), list(gr2 = gr2), nperm = 0)
#
genomic_regions_correlation = function(gr_list_1, gr_list_2, background = NULL,
	chromosome = paste0("chr", 1:22), species = "hg19",
	nperm = 0, mc.cores = 1, stat_fun = genomic_corr_jaccard, ..., 
	bedtools_binary = Sys.which("bedtools"), tmpdir = tempdir()) {
	
	qq.options(LOCAL = TRUE)
	qq.options(code.pattern = "@\\{CODE\\}")
	
	## check input
	if(inherits(gr_list_1, "GRanges")) {
		gr_name_1 = deparse(substitute(gr_list_1))
		gr_list_1 = list(gr_list_1)
		names(gr_list_1) = gr_name_1
	}
	if(inherits(gr_list_2, "GRanges")) {
		gr_name_2 = deparse(substitute(gr_list_2))
		gr_list_2 = list(gr_list_2)
		names(gr_list_2) = gr_name_2
	}
	
	if(is.null(names(gr_list_1))) {
		stop("`gr_list_1` should have names.\n")
	}
	if(is.null(names(gr_list_2))) {
		stop("`gr_list_2` should have names.\n")
	}

	if(nperm > 1) {
		if(!file.exists(bedtools_binary)) {
			stop("Cannot find binary file for bedtools.")
		}
	}

	message("set strand to * and merge potential overlapped regions in gr_list_1.")
	gr_list_1 = lapply(gr_list_1, function(gr) {
		strand(gr) = "*"
		reduce(sort(gr))
	})
	message("set strand to * and merge potential overlapped regions in gr_list_2.")
	gr_list_2 = lapply(gr_list_2, function(gr) {
		strand(gr) = "*"
		reduce(sort(gr))
	})
	
	# limit in chromosomes
	message("subset regions in selected chromosomes.")
	gr_list_1 = lapply(gr_list_1, function(gr) gr[ seqnames(gr) %in% chromosome])
	gr_list_2 = lapply(gr_list_2, function(gr) gr[ seqnames(gr) %in% chromosome])

	gr_name_1 = names(gr_list_1)
	gr_name_2 = names(gr_list_2)
	
	# limit in background
	if(!is.null(background)) {
		strand(background) = "*"
		background = background[ seqnames(background) %in% chromosome ]
		background = reduce(background)

		message("overlaping `gr_list_1` to background")
		gr_list_1 = lapply(gr_list_1, function(gr) {
			intersect(gr, background)
		})
		message("overlaping `gr_list_2` to background")
		gr_list_2 = lapply(gr_list_2, function(gr) {
			intersect(gr, background)
		})

		# since `background` will be send to bedtoos, here we need a data frame
		background_df = as.data.frame(background)[1:3]
	}

	chr_len_df = getChromInfoFromUCSC(species)
	chr_len_df = chr_len_df[chr_len_df[[1]] %in% chromosome, , drop = FALSE]  # needed for bedtools shuffle

	# prepare values that will be returned
	foldChange = matrix(0, nrow = length(gr_list_2), ncol = length(gr_list_1))
	rownames(foldChange) = names(gr_list_2)
	colnames(foldChange) = names(gr_list_1)
	p = foldChange
	stat = foldChange
	stat_random_mean = stat
	stat_random_sd = stat

	# cache chr_len_df and background files
	chr_len_df_tmp = tempfile(tmpdir = tmpdir)
	write.table(chr_len_df, file = chr_len_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)
	if(!is.null(background)) {
		background_df_tmp = tempfile(tmpdir = tmpdir)
		write.table(background_df, file = background_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)
	}

	# loop in gr_list_1
	for(i in seq_along(gr_list_1)) {

		stat_random = matrix(0, nrow = length(gr_list_2), ncol = nperm)
		
		# stat
		for(j in seq_along(gr_list_2)) {
			message(qq("calculating correlation between @{gr_name_1[i]} and @{gr_name_2[[j]]}"))
			stat[j, i] = do.call("stat_fun", list(gr_list_1[[i]], gr_list_2[[j]], ...))
		}

		if(nperm > 1) {
			# random shuffle gr_list_1
			# cache gr_list_1
			gr_list_1_df = as.data.frame(gr_list_1[[i]])[1:3]
			gr_list_1_df_tmp = tempfile()
			write.table(gr_list_1_df, file = gr_list_1_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)

			res = mclapply(seq_len(nperm), mc.cores = mc.cores, function(k) {
				
				if(is.null(background)) {
					gr_random = systemdf(qq("'@{bedtools_binary}' shuffle -i '@{gr_list_1_df_tmp}' -g '@{chr_len_df_tmp}'"))
				} else {
					gr_random = systemdf(qq("'@{bedtools_binary}' shuffle -i '@{gr_list_1_df_tmp}' -g '@{chr_len_df_tmp}' -incl '@{background_df_tmp}'"))
				}

				gr_random = GRanges(seqnames = gr_random[[1]], ranges = IRanges(gr_random[[2]], gr_random[[3]]))

				# x contains stat for every gr_list_2
				x = numeric(length(gr_list_2))
				for(j in seq_along(gr_list_2)) {
					
					message(qq("calculating correlation between random_@{gr_name_1[i]} and @{gr_name_2[[j]]}, @{k}/@{nperm}"))

					x[j] = do.call("stat_fun", list(gr_random, gr_list_2[[j]], ...))
				}

				return(x)
			})

			# `res` is a list, convert to a matrix
			for(k in seq_along(res)) {
				stat_random[, k] = res[[k]]
			}

			file.remove(gr_list_1_df_tmp)

			stat_random_mean[, i] = rowMeans(stat_random)
			stat_random_sd[, i] = rowSds(stat_random)

			foldChange[, i] = stat[, i]/stat_random_mean[, i]
			p[, i] = sapply(seq_len(length(gr_list_2)), function(i) sum(stat_random[i, ] > stat[i])/nperm)
		} else {
			stat_random_mean[, i]  = NA
			stat_random_sd[, i] = NA
			foldChange[, i] = NA
			p[, i] = NA
		}
	}

	file.remove(chr_len_df_tmp)
	if(!is.null(background)) {
		file.remove(background_df_tmp)
	}

	if(nperm <= 1) {
		foldChange = NULL
		p = NULL
		stat_random_mean = NULL
		stat_random_sd = NULL
	}

	res = list(stat = stat,
		       fold_change = foldChange, 
		       p.value = p,
		       stat_random_mean = stat_random_mean,
		       stat_random_sd = stat_random_sd)
	
	qq.options(LOCAL = FALSE)
	return2(res)
}

# == title
# Relative distance between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
#
# == details
# For regions in ``query`` and ``reference``, they are all degenerated as single points
# which are the middle points of regions. For each middle point in ``query``, it looks 
# for two nearest points in ``reference`` on its left and right. The statistic is defined as the ratio of the distance
# to the nearest neighbour point to the distance of two neighbour points. If ``query`` and ``reference`` are not correlated at all,
# It is expected that the ratio follows a uniform distribution. So final statisitic are the KS-statistic
# between the real distribution of rations to the uniform distribution.
#
# == reference
# Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
# 
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_reldist(gr1, gr2)
#
genomic_corr_reldist = function(query, reference) {
	# GRanges for mid-points
	query_mid = ceiling( (start(query) + end(query))/2 )
	reference_mid = ceiling( (start(reference) + end(reference))/2 )
	
	# we don't need strand information here
	query_mid_gr = GRanges(seqnames = seqnames(query),
		                   ranges = IRanges(start = query_mid,
		                   	                end = query_mid))
	reference_mid_gr = GRanges(seqnames = seqnames(reference),
		                   ranges = IRanges(start = reference_mid,
		                   	                end = reference_mid))
	
	reference_mid_gr = sort(reference_mid_gr)
	
	# look for nearest position
	mtch = nearest(query_mid_gr, reference_mid_gr) # index in reference
	l = !is.na(mtch)        # if there is no site in reference on some chromosome, the value would be NA
	query_mid_gr = query_mid_gr[l]   # remove these regions in query
	m2 = reference_mid_gr[ mtch[l] ]
	mtch = mtch[l]
	
	# look for another point in reference
	ind = ifelse(start(query_mid_gr) > start(m2), mtch + 1, mtch - 1)  # ind contains index in reference
	l = ind >= 1 & ind <= length(reference)
	m3 = reference_mid_gr[ ind[l] ]
	query_mid_gr = query_mid_gr[l]
	m2 = m2[l]

	suppressWarnings(d <- distance(query_mid_gr, m2) / distance(m2, m3))
	d = d[!is.na(d) & !is.infinite(d)]  # this is not necessary, just double check
	
	stat = 0
	try(stat <- (integrate(ecdf(d), lower = 0, upper = 0.5)$value - 0.25) / 0.25, silent = TRUE)
	
	return(stat)
}

# == title
# Jaccard coefficient between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
# -background background regions that should be only looked in, a `GenomicRanges::GRanges` object
#
# == details
# Jaccard coefficient is defined as the total length of intersection divided by total
# length of union of two sets of genomic regions.
#
# You can set the background when calculating Jaccard coefficient. For example,
# if the interest is the Jaccard coefficient between CpG sites in ``query`` and in ``reference``
# ``background`` can be set with a `GenomicRanges::GRanges` object which contains positions of CpG sites.
#
# Be careful with the ``strand`` in your `GenomicRanges::GRanges` object!
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_jaccard(gr1, gr2)
#
genomic_corr_jaccard = function(query, reference, background = NULL) {
	if(is.null(background)) {
		res = sum(width(intersect(query, reference))) / sum(as.numeric(width(union(query, reference))))
	} else {
		gr1 = intersect(query, reference)
		gr1 = intersect(gr1, background)

		gr2 = union(query, reference)
		gr2 = intersect(gr2, background)
		res = sum(width(gr1)) / sum(width(gr2))
	}
	return(res)
}

# == title
# Absolute distance between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
# -method function in which input is a vector of distance and output is a scalar
# -... pass to ``method``
#
# == details
# For regions in ``query`` and ``reference``, they are all degenerated as single points
# which are the middle points of corresponding regions. For each middle point in ``query``, it looks 
# for the nearest point in ``reference``. Assuming the distance vector is ``d``, the final statistic is ``method(d)``.
#
# == reference
# Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_absdist(gr1, gr2)
#
genomic_corr_absdist = function(query, reference, method = mean, ...) {
	# GRanges for mid-points
	query_mid = ceiling( (start(query) + end(query))/2 )
	reference_mid = ceiling( (start(reference) + end(reference))/2 )
	
	query_mid_gr = GRanges(seqnames = seqnames(query),
		                   ranges = IRanges(start = query_mid,
		                   	                end = query_mid))
	reference_mid_gr = GRanges(seqnames = seqnames(reference),
		                   ranges = IRanges(start = reference_mid,
		                   	                end = reference_mid))
	
	# look for nearest position
	suppressWarnings(mtch <- distanceToNearest(query_mid_gr, reference_mid_gr))
	dst = mtch@elementMetadata[, 1]
	dst = dst[!is.na(dst)]
	stat = method(dst, ...)
	
	return(stat)
}

# == title
# Intersections between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
# -... pass to `GenomicRanges::countOverlaps`
#
# == details
# It calculates number of regions in ``query`` that overlap with ``reference``.
#
# Please note this value is not equal to the number of intersections betweenn two sets of regions,
# because one region in ``query`` may overlap with more than one
# regions in ``reference``
#
# Be careful with the ``strand`` in your GRanges object!
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_nintersect(gr1, gr2)
genomic_corr_nintersect = function(query, reference, ...) {
	x = countOverlaps(query, reference, ...)
	res = sum(x > 0)
	return(res)
}

# == title
# Intersection between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
# -... pass to `percentOverlaps`
#
# == details
# For each region in ``query``, it calculates the percent that is covered by ``reference``.
#
# The returned value is percent which is how much ``query`` is covered by ``reference``. 
#
# Be careful with the ``strand`` in your GRanges object!
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_pintersect(gr1, gr2)
genomic_corr_pintersect = function(query, reference, ...) {
	x = percentOverlaps(query, reference, ...)
	
	w = width(query)
	res = sum(x*w)/sum(w)
	
	return(res)
}

# == title
# Intersection between two sets of genomic regions
#
# == param
# -query genomic region 1, a `GenomicRanges::GRanges` object
# -reference genomic region 2, a `GenomicRanges::GRanges` object
# -background subset of regions that should be only looked into, a `GenomicRanges::GRanges` object
#
# == details
# It calculates the total length of overlapped regions in ``query``.
#
# Be careful with the ``strand`` in your GRanges object!
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_sintersect(gr1, gr2)
genomic_corr_sintersect = function(query, reference, background = NULL) {
	x = intersect(query, reference)
	if(!is.null(background)) {
		x = intersect(x, background)
	}
	res = sum(width(x))
	return(res)
}
