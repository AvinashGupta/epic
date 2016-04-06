make_transition_matrix_from_chromHMM = function(gr_list_1, gr_list_2, window = NULL, 
	all_states = NULL, min_1 = round(length(gr_list_1)/2), min_2 = round(length(gr_list_2)/3)) {

	if(inherits(gr_list_1, "GRanges")) {
		gr_list_1 = list(gr_list_1)
		min_1 = 0
	}
	if(inherits(gr_list_2, "GRanges")) {
		gr_list_2 = list(gr_list_2)
		min_2 = 0
	}

	if(is.null(window)) {
		window = mGCD(c(sapply(gr_list_1, function(gr) mGCD(unique(width(gr)))),
		                sapply(gr_list_2, function(gr) mGCD(unique(width(gr))))))
		message(paste0("window is set to ", window))
	}

	if(is.null(all_states)) {
		all_states = unique(unlist(lapply(gr_list_1, function(gr) {unique(as.character(mcols(gr)[, 1]))})),
			                unlist(lapply(gr_list_2, function(gr) {unique(as.character(mcols(gr)[, 1]))})))
		all_states = sort(all_states)
		message(paste0(length(all_states), " states in total"))
	}

	message("extract states")
	m1 = as.data.frame(lapply(gr_list_1, function(gr) {
		k = round(width(gr)/window)
		s = as.character(mcols(gr)[, 1])
		as.integer(factor(rep(s, times = k), levels = all_states))
	}))
	m2 = as.data.frame(lapply(gr_list_2, function(gr) {
		k = round(width(gr)/window)
		s = as.character(mcols(gr)[, 1])
		as.integer(factor(rep(s, times = k), levels = all_states))
	}))

	m1 = as.matrix(m1)
	m2 = as.matrix(m2)

	message("count for each state")
	t1 = rowTabulates(m1)
	t2 = rowTabulates(m2)
	l = rowMaxs(t1) >= min_1 & rowMaxs(t2) >= min_2
	t1 = t1[l, , drop = FALSE]
	t2 = t2[l, , drop = FALSE]

	message("determine states")
	states1 = rowWhichMax(t1)
	states2 = rowWhichMax(t2)

	message("generate transition matrix")
	mat = as.matrix(table(states1, states2))
	rownames(mat) = all_states[as.numeric(rownames(mat))]
	colnames(mat) = all_states[as.numeric(colnames(mat))]

	return(mat)
}



# == title
# chord diagram for chromatin states transistion
#
# == param
# -mat the transition matrix
# -max_mat if there are several matrix, set it to the matrix with maximum sum
# -cate1 name of row states
# -cate2 name fo column states
# -grid.col color for grids
# -... pass to `circlize::chordDiagram`
#
chromatin_states_transition_chord_diagram = function(mat, max_mat = mat, cate1, cate2, grid.col, ...) {

	par(xpd = NA)

	tb = mat
	n = nrow(mat)
	rownames(tb) = 1:n
	colnames(tb) = 1:n
	rownames(tb) = paste0("R", 1:n)
	colnames(tb) = paste0("C", 1:n)
	colmat = rep(grid.col, n); dim(colmat) = dim(tb); colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
	qati = quantile(tb, 0.7)
	colmat[tb > qati] = paste0(colmat[tb > qati], "A0")
	colmat[tb <= qati] = paste0(colmat[tb <= qati], "20")
	dim(colmat) = dim(tb);

	de = 360 - (360 - 20 - 30) * sum(mat)/sum(max_mat) - 30
	circos.par(start.degree = -de/4, gap.degree = c(rep(1, n-1), de/2, rep(1, n-1), de/2))
	gcd = rep(grid.col, 2); names(gcd) = c(rownames(tb), colnames(tb))
	chordDiagram(tb, col = colmat, grid.col = gcd,
		directional = TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), ...)

	text(-1, 1, qq("states\nin @{cate1}"), adj = c(0, 1), cex = 1.5)
	text(1, -1, qq("states\nin @{cate2}"), adj = c(1, 0), cex = 1.5)

	for(sn in get.all.sector.index()) {
		if(abs(get.cell.meta.data("cell.start.degree", sector.index = sn) - get.cell.meta.data("cell.end.degree", sector.index = sn)) > 3) {
			circos.text(get.cell.meta.data("xcenter", sector.index = sn, track.index = 2), get.cell.meta.data("ycenter", sector.index = sn, track.index = 2), 
				gsub("C|R", "", sn), col = "white", font = 2, sector.index = sn, track.index = 2, adj = c(0.5, 0.5), niceFacing = TRUE)
			circos.axis(sector.index = sn, track.index = 2, major.tick.percentage = 0.2, labels.away.percentage = 0.2, labels.cex = 0.5)
		}
	}

	circos.clear()
}
