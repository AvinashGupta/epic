### functions related to QC
## the functions make related plots and output some tables

### mat_density_plot for coverage and methylation, add mean and median lines on the heatmap


# == title
# Basic qc plot for bisulfite sequencing data
#
# == param
# -sample_id a vector of sample ids
# -chromosome a vector of chromosome names
#
# == detail
# For each sample id, it will produce five plots:
#
# - mean/median CpG coverage per chromosome
# - histogram of CpG coverage
# - methylation per chromosome 
# - histogram of methylation
# - mean Methylation for each CpG coverage 
#
# == value
# A list of corresponding statistics
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
wgbs_qcplot = function(sample_id, chromosome = paste0("chr", 1:22), background = NULL) {

	# coverage and methylation per chromosome
	data = rep(list(list(cov = NULL, meth = NULL, strand = NULL, cov_count = NULL)), length(sample_id))
	names(data) = sample_id
	
	for(chr in chromosome) {
		
		methylation_hooks$set(chr)
		if(!is.null(background)) {
			mtch = as.matrix(findOverlaps(methylation_hooks$GRanges(), background))
		}
		for(sid in sample_id) {

			cv = methylation_hooks$coverage(col_index = sid)
			ind = which(cv != 0)
			data[[sid]]$cov_count[[chr]] = c("zero" = length(cv) - length(ind), "all" = length(cv))
			if(!is.null(background)) {
				ind = intersect(ind, mtch[, 1])
			}
			if(length(ind) > 10000) {
				ind = sort(sample(ind, 10000))
			}
			cv = cv[ind]

			mh = as.vector(methylation_hooks$meth(row_index = ind, col_index = sid))
			
			strd = rep("*", length(mh))
			
			data[[sid]]$cov[[chr]] = cv
			data[[sid]]$meth[[chr]] = mh
			data[[sid]]$strand[[chr]] = strd
		}
	}

	for(sid in sample_id) {

		cov = data[[sid]]$cov
		meth = data[[sid]]$meth
		strand = data[[sid]]$strand
		cov_count = data[[sid]]$cov_count

		par(mfrow = c(2, 3))
		
		# mean coverage per chromosome
		cpg_coverage_mean = sapply(cov, mean)
		cpg_coverage_median = sapply(cov, median)
		plot(c(0, length(cpg_coverage_mean)), c(0, max(c(cpg_coverage_mean, cpg_coverage_median))), axes = FALSE, ann = FALSE, type="n")
		for(i in seq_along(cpg_coverage_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_coverage_mean[i], cpg_coverage_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_coverage_median[i], cpg_coverage_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_coverage_mean)-0.5, labels = names(cpg_coverage_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("Coverage per chromosome (@{sid})"), ylab = "mean and median CpG coverage")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		
		# coverage distribution
		cov_count = do.call("cbind", cov_count)
		zero_cov_rate = sprintf("%.2f", sum(cov_count[1, ])/sum(cov_count[2, ])*100)
		if(all(unique(unlist(strand)) %in% c("+", "-"))) {
			x = unlist(cov)
			q99 = quantile(unlist(cov), 0.99)
			y = unlist(strand)
			x1 = x[y == "+"]
			x2 = x[y == "-"]
			ta = table(x)
			ta1 = table(x1)
			ta2 = table(x2)
			xlim = range(c(as.numeric(names(ta)), as.numeric(names(ta1)), as.numeric(names(ta2))))
			ylim = range(ta, ta1, ta2)
			plot(as.numeric(names(ta)), ta, xlim = xlim, ylim = ylim, main = qq("histogram of CpG coverage (@{sid})\n@{zero_cov_rate}% have zero coverage"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
			par(new = TRUE)
			plot(as.numeric(names(ta1))+0.2, ta1, xlim = xlim, ylim = ylim, type = "h", col = "red", log = "x", axes = FALSE, ann = FALSE)
			par(new = TRUE)
			plot(as.numeric(names(ta2))+0.4, ta2, xlim = xlim, ylim = ylim, type = "h", col = "blue", log = "x", axes = FALSE, ann = FALSE)
			#axis(side = 2)
			legend("topright", lty = 1, col = c("black", "red", "blue"), legend = c("strand *", "strand +", "strand -"))
			par(new = FALSE)
		} else {
			ta = table(unlist(cov))
			q99 = quantile(unlist(cov), 0.99)
			plot(as.numeric(names(ta)), ta, main = qq("histogram of CpG coverage (@{sid})\n@{zero_cov_rate}% have zero coverage"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			abline(v = q99, lty = 2, col = "blue"); text(q99, 0, "q99", adj = c(0, 0))
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
		}

		# ## barplot of zero-coverage and non-zero-coverage
		# plot(c(0, length(cov_count)), c(0, max(unlist(cov_count))), axes = FALSE, ann = FALSE, type = "n")
		# for(i in seq_along(cov_count)) {
		# 	rect(i-1, 0, i, cov_count[[i]]["zero"], col = "orange")
		# 	rect(i-1, cov_count[[i]]["zero"], i, cov_count[[i]]["all"], col = "blue")
		# 	abline(v = i, lty = 3, col = "grey")
		# }
		# abline(v = 0, lty = 3, col = "grey")
		# par(las = 3)
		# axis(side = 1, at = seq_along(cov_count) - 0.5, labels = names(cov_count))
		# axis(side = 2)
		# box()
		# par(las = 0)
		# title(main = qq("%zeor/non_zero CpG coverage (@{sid})"), ylab = "number of CpG sites")
		# legend("topright", pch = 15, col = c("orange", "blue"), legend = c("zero", "non-zero"))

		# mean methylation per chromosome
		cpg_methyrate_mean = sapply(meth, mean)
		cpg_methyrate_median = sapply(meth, median)
		plot(c(0, length(cpg_methyrate_mean)), c(0, 1), axes = FALSE, ann = FALSE, type = "n")
		for(i in seq_along(cpg_methyrate_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_methyrate_mean[i], cpg_methyrate_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_methyrate_median[i], cpg_methyrate_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_methyrate_mean) - 0.5, labels = names(cpg_methyrate_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("methylation per chromosome (@{sid})"), ylab = "mean and median methylation")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		# distribution of methylation on all chromosomes
		hist(unlist(meth), main = qq("histogram of methylation (@{sid})"), xlab = "methylation")

		
		# methylation to coverage
		coverage2methyrate = tapply(unlist(meth), unlist(cov), mean)
		plot(as.numeric(names(coverage2methyrate)), coverage2methyrate, ylim = c(0, 1), pch=16, log = "x", cex = 0.8, xlab = "CpG coverage", ylab = "mean methylation", main = qq("Mean Methylation for each CpG coverage (@{sid})"))
		coverage2methyrate = tapply(unlist(meth), unlist(cov), median)
		points(as.numeric(names(coverage2methyrate)), coverage2methyrate, pch=16, cex = 0.8, col = "red")
		legend("bottomleft", pch = 16, col = c("black", "red"), legend = c("mean", "median"))
		abline(v = q99, lty = 2, col = "blue"); text(q99, 0, "q99 of cov", adj = c(0, 0))
		par(mfrow = c(1, 1))
	}
	
	return2(data, invisible = TRUE)
}

# == title
# Plot coverage and methylation for a single sample
#
# == param
# -sid a single sample id
# -chromosome a vector of chromosome names
# -species species
# -nw number of windows
# -... pass to `gtrellis::gtrellis_layout`
#
# == details
# The whole genome is segented by ``nw`` windows and mean methylation and mean CpG coverage
# are visualized as two tracks.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
plot_coverage_and_methylation_on_genome = function(sid, chromosome = paste0("chr", 1:22), 
	species = "hg19", nw = 10000, ...) {

	w = round(read.chromInfo(species = species)$chr.len["chr1"]/nw)
	flag = 0
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"), transparency = 0.8)
	for(chr in chromosome) {
		methylation_hooks$set(chr)

		meth = methylation_hooks$meth(col_index = sid)[,1]
		cov = methylation_hooks$coverage(col_index = sid)[,1]; cov = log10(cov+1)
		site = methylation_hooks$site()
		gr = methylation_hooks$GRanges()
		chr_len = read.chromInfo(species = species)$chr.len[chr]
		chr_gr = GRanges(seqnames = chr, ranges = IRanges(1, chr_len))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(mtch[,2], mtch[,1], function(i) mean(meth[i]))
		cov = tapply(mtch[,2], mtch[,1], function(i) mean(cov[i]))
		
		if(flag == 0) {
			gtrellis_layout(category = chromosome, species = species, 
				n_track = 2, track_ylim = c(0, quantile(cov, 0.95), 0, 1), 
				track_ylab = c("log10(coverage)", "methylation"),
				add_name_track = TRUE, add_ideogram_track = TRUE, ...)
			flag = 1
		}

		add_track(gr2, track = 2, category = chr, panel.fun = function(gr) {
			x = (start(gr) + end(gr))/2
			y = cov
			grid.points(x, y, pch = ".", gp = gpar(col = "#FF000010"))
		})
		add_track(gr2, track = 3, category = chr, panel.fun = function(gr) {
			x = (start(gr) + end(gr))/2
			y = meth
			grid.points(x, y, pch = ".", gp = gpar(col = col_fun(y)))
		})
	}
}




# == title
# Plot methylation for multiple samples
#
# == param
# -sample_id a vector of sample ids
# -annotation annotation of samples (e.g. subtypes)
# -chromosome a vector of chromosome names
# -species species
# -nw number of windows
# -... pass to `gtrellis::gtrellis_layout`
#
# == details
# The whole genome is segented by ``nw`` windows. Methylation in different classes are visualized as separated tracks.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
plot_multiple_samples_methylation_on_genome = function(sample_id, annotation, 
	chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000, ...) {
	
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"), transparency = 0.8)
	type = unique(annotation)
		
	n = table(annotation)[type]
	ty = numeric(2*length(n))
	ty[seq_along(n)*2-1] = 0.5
	ty[seq_along(n)*2] = n + 0.5

	gtrellis_layout(category = chromosome, species = species, 
		n_track = length(type), track_ylab = type, track_ylim = ty, track_height = n,
		add_name_track = TRUE, add_ideogram_track = TRUE, ...)

	w = round(read.chromInfo(species = species)$chr.len["chr1"]/nw)
	for(chr in chromosome) {
		methylation_hooks$set(chr)

		meth = methylation_hooks$meth(col_index = sample_id)
		site = methylation_hooks$site()
		gr = methylation_hooks$GRanges()
		chr_len = read.chromInfo(species = species)$chr.len[chr]
		chr_gr = GRanges(seqnames = chr, ranges = IRanges(1, chr_len))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth[i, , drop = FALSE]))

		meth = do.call("rbind", meth)

		for(i in seq_along(type)) {
			sid = sample_id[annotation == type[i]]
			m = meth[, sid, drop = FALSE]
			add_track(gr2, track = i+1, category = chr, panel.fun = function(gr) {
				x = (start(gr2) + end(gr2))/2
				for(i in seq_along(sid)) {
					y = rep(i, length(x)) + (runif(length(x))-0.5)*0.8
					grid.points(x, y, pch = ".", gp = gpar(col = col_fun(m[, i])))
				}
			})
		}
	}
}

# -x a data frame or a list
.mat_dist = function(x, anno, col, reorder_column = TRUE, od = seq_along(anno), ha = NULL, title = NULL, ...) {

	if(reorder_column) {
		if(inherits(x, "list")) {
			od = order(factor(anno, levels = unique(anno), ordered = TRUE), sapply(x, median, na.rm = TRUE))
		} else {
			od = order(factor(anno, levels = unique(anno), ordered = TRUE), apply(x, 2, median, na.rm = TRUE))
		}
	} 
	if(is.null(ha)) ha = HeatmapAnnotation(df = data.frame(type = anno), col = list(type = col))

	densityHeatmap(x, anno = ha, title = title, column_order = od, ...)
	for(an in sapply(ha@anno_list, function(x) x@name)) {
		decorate_annotation(an, {
			grid.text(an, x = unit(-2, "mm"), just = "right")
		})
	}

	par(xpd = FALSE)
	if(is.matrix(x)) {
		den_x = matrix(nrow = 512, ncol = dim(x)[2])
		den_y = matrix(nrow = 512, ncol = dim(x)[2])
		for(i in seq_len(dim(x)[2])) {
			den = density(x[, i], na.rm = TRUE)
			den_x[, i] = den$x
			den_y[, i] = den$y
		}
	} else {
		den_x = matrix(nrow = 512, ncol = length(x))
		den_y = matrix(nrow = 512, ncol = length(x))
		for(i in seq_along(x)) {
			den = density(x[[i]], na.rm = TRUE)
			den_x[, i] = den$x
			den_y[, i] = den$y
		}
	}
	
	matplot(den_x, den_y, type = "l", col = col[anno], xlab = "value", ylab = "density", main = qq("density distribution: @{title}"))

	## MDS plot
	if(is.data.frame(x) || is.matrix(x)) {
		mat = as.matrix(x)
		loc = cmdscale(dist2(t(mat), pairwise_fun = function(x, y) {l = is.na(x) | is.na(y); x = x[!l]; y = y[!l]; sqrt(sum((x-y)^2))}))
		plot(loc[, 1], loc[, 2], pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
		legend("bottomleft", pch = 16, legend = names(col), col = col)

		plot(loc[, 1], loc[, 2], type = "n", pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
		text(loc[, 1], loc[, 2], colnames(x), col = col[anno], cex = 0.8)
	}

	return(od)
}


# == title
# Global methylation distribution
# 
# == param
# -sample_id a vector of sample ids
# -annotation classification information
# -annotation_color color for classifications
# -reorder_column whether reorder the samples
# -ha additional annotation can be specified as a `ComplexHeatmap::HeatmapAnnotation` object
# -chromosome chromosomes
# -by_chr whether make the plot by chromosome
# -max_cov maximum coverage (used to get rid of extremely high coverage which affects visualization of CpG coverage distribution)
# -background background to look into. The value can be a single `GenomicRanges::GRanges` object or a list of `GenomicRanges::GRanges` objects.
# -p probability to randomly sample CpG sites
#
# == details
# It visualize distribution of methylation values and CpG coverages through heatmaps.
#
# == value
# If ``by_chr`` is set to ``FALSE``, it returns a vector of column order.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
global_methylation_distribution = function(sample_id, annotation, 
	annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
	reorder_column = TRUE, ha = NULL, chromosome = paste0("chr", 1:22), by_chr = FALSE, max_cov = 100,
	background = NULL, p = 0.001, meth_range = c(0, 1)) {

	annotation_color = annotation_color[intersect(names(annotation_color), unique(annotation))]
	
	###############################################
	# distribution of global methylation
	if(inherits(background, "list")) {
		meth_list = NULL
		cov_list = NULL

		if(length(background) != length(sample_id)) {
			stop("Since you specified `background` as a list, the length should be same as `sample_id`.")
		}

		for(chr in chromosome) {

			methylation_hooks$set(chr)

			meth_gr = methylation_hooks$GRanges()
			ind_list = lapply(seq_along(sample_id), function(i) {
				mtch = as.matrix(findOverlaps(meth_gr, background[[i]]))
				ind = unique(mtch[, 1])
				nr = length(ind)
				ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]

				message(qq("random sampled @{length(ind)} sites from @{nr} sites on @{chr} in @{sample_id[i]} (with p = @{p})"))
			})

			current_meth_list = lapply(seq_along(ind_list), function(i) {
				m = methylation_hooks$meth(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				cov = methylation_hooks$coverage(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				m[cov == 0] = NA
				m
			})
			current_cov_list = lapply(seq_along(ind_list), function(i) {
				cov = methylation_hooks$coverage(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				cov[cov == 0] = NA
				cov[cov > max_cov] = NA
				log10(cov)
			})

			message(qq("on average there are @{round(mean(sapply(current_meth_list, function(x) sum(is.na(x)))))} CpG with 0 coverage.\n"))

			
			if(by_chr) {
				try(od <- .mat_dist(current_meth_list, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				try(.mat_dist(current_cov_list, reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}

			meth_list = lapply(seq_along(meth_list), function(i) {
				c(meth_list[[i]], current_meth_list[[i]])
			})
			cov_list = lapply(seq_along(cov_list), function(i) {
				c(cov_list[[i]], current_cov_list[[i]])
			})
		}

		if(!by_chr) {
			meth_list = lapply(meth_list, function(meth) {
				n = length(meth)
				if(n > 100000) {
					meth[sample(n, 100000)]
				}
			})
			cov_list = lapply(cov_list, function(cov) {
				n = length(cov)
				if(n > 100000) {
					cov[sample(n, 100000)]
				}
			})
			
			od = .mat_dist(meth_list, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			.mat_dist(cov_list, reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			
			return(invisible(od))
		}
	} else {
		meth_mat = NULL
		cov_mat = NULL
		for(chr in chromosome) {

			methylation_hooks$set(chr)

			meth_gr = methylation_hooks$GRanges()
			if(!is.null(background)) {
				mtch = as.matrix(findOverlaps(meth_gr, background))
				ind = unique(mtch[, 1])
			} else {
				ind = seq_len(length(meth_gr))
			}
			
			nr = length(ind)
			ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]
			
			message(qq("random sampled @{length(ind)} sites from @{nr} sites on @{chr} (with p = @{p})"))
			mm = methylation_hooks$meth(row_index = ind, col_index = sample_id)
			cm = methylation_hooks$coverage(row_index = ind, col_index = sample_id)
			mm[cm == 0] = NA
			cm[cm == 0] = NA
			cm[cm > max_cov] = NA
			
			message(qq("on average there are @{round(mean(apply(mm, 2, function(x) sum(is.na(x)))))} CpG with 0 coverage.\n"))

			meth_mat = rbind(meth_mat, mm)
			cov_mat = rbind(cov_mat, cm)

			if(by_chr) {
				
				try(od <- .mat_dist(mm, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				try(.mat_dist(log10(cm), reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}
		}

		if(!by_chr) {
			nr = nrow(meth_mat)
			if(nr > 100000) {
				meth_mat = meth_mat[sample(nr, 100000), ]
				cov_mat = cov_mat[sample(nr, 100000), ]
			}
			
			od = .mat_dist(meth_mat, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			.mat_dist(log10(cov_mat), reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			
			return(invisible(od))
		}
	}
}
