# == title
# Generate transition matrix from chromHMM results
#
# == param
# -gr_list_1 a list of `GenomicRanges::GRanges` objects which contain chromatin states in group 1.
#            The first column in meta columns should be the states.
#            please note the start position when converting bed format to ``GRanges`` format (0-based and 1-based).
# -gr_list_2 a list of `GenomicRanges::GRanges` objects which contains chromatin states in group 2.
# -window window size which was used to do chromHMM states prediction. If it is not specified, the greatest common divisor
#         of all region width is used.
# -min_1 minimal recurrency in ``gr_list_1``
# -min_2 minimal recurrency in ``gr_list_2``
#
# == value
# A transition matrix in which values represent width of regions that transite from one state to the other.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
make_transition_matrix_from_chromHMM = function(gr_list_1, gr_list_2, window = NULL, 
	min_1 = floor(length(gr_list_1)/2), min_2 = floor(length(gr_list_2)/2)) {

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
		if(window == 1) {
			message("when converting bed files to GRanges objects, be careful with the 0-based and 1-based coordinates.")
		}
	}

	# if(is.null(all_states)) {
		all_states = unique(c(unlist(lapply(gr_list_1, function(gr) {unique(as.character(mcols(gr)[, 1]))})),
			                unlist(lapply(gr_list_2, function(gr) {unique(as.character(mcols(gr)[, 1]))}))))
		all_states = sort(all_states)
		message(paste0(length(all_states), " states in total"))
	# }

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
	mat = mat * as.numeric(window)
	class(mat) = "matrix"
	return(mat)
}



# == title
# Chord diagram for chromatin states transistion
#
# == param
# -mat the transition matrix
# -max_mat if there are several matrix to be compared, set it to the matrix with maximum sum
# -remove_unchanged_transition whether to remove regions that states are not changed
# -state_col color for grids. It should be a vector of which names correspond to states
# -legend_position positinos of legends. Possible values are "bottomleft", "bottomright", "topright" and "topleft".
#             If the value is specified as vector with length two, the legend will be split into two parts.
#              Set the value to ``NULL`` to suppress legend.
# -... pass to `circlize::chordDiagram`
#
# == details
# Rows of ``mat`` always locates at the bottom.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
chromatin_states_transition_chord_diagram = function(mat, max_mat = mat, remove_unchanged_transition = TRUE,
	state_col = NULL, legend_position = NULL, ...) {

	op = par(no.readonly = TRUE)
	on.exit(par(op))
	par(xpd = NA)

	mat = as.matrix(mat)

	if(length(intersect(rownames(mat), colnames(mat))) == 0) {
		rownames(mat) = paste0("E", seq_len(nrow(mat)))
		colnames(mat) = paste0("E", seq_len(ncol(mat)))
	}

	all_states = union(rownames(mat), colnames(mat))
	n_states = length(all_states)
	mat2 = matrix(0, nrow = n_states, ncol = n_states)
	rownames(mat2) = all_states
	colnames(mat2) = all_states
	mat2[rownames(mat), colnames(mat)] = mat  

	mat = mat2

	if(is.null(state_col)) {
		col_fun = colorRamp2(0:10, brewer.pal(11, "Spectral"))
		x = (seq_along(all_states)-1)*10/(length(all_states)-1)
		state_col = col_fun(x)
		names(state_col) = all_states
	} else {
		if(is.null(names(state_col))) {
			names(state_col) = all_states
		}
	}

	if(remove_unchanged_transition) {
		for(i in all_states) {
			mat[i, i] = 0
		}
	}

	rownames(mat) = paste0("R_", seq_len(n_states))
	colnames(mat) = paste0("C_", seq_len(n_states))

	state_col2 = c(state_col, state_col)
	names(state_col2) = c(rownames(mat), colnames(mat))

	colmat = rep(state_col2[rownames(mat)], n_states)
	colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
	
	# make thin links very light
	qati = quantile(mat, 0.7)
	colmat[mat > qati] = paste0(colmat[mat > qati], "A0")
	colmat[mat <= qati] = paste0(colmat[mat <= qati], "20")
	dim(colmat) = dim(mat)

	de = 360 - (360 - 20 - 30) * sum(mat)/sum(max_mat) - 30
	circos.par(start.degree = -de/4, gap.degree = c(rep(1, n_states-1), de/2, rep(1, n_states-1), de/2))

	chordDiagram(mat, col = colmat, grid.col = state_col2,
		directional = TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), ...)

	for(sn in get.all.sector.index()) {
		if(abs(get.cell.meta.data("cell.start.degree", sector.index = sn) - get.cell.meta.data("cell.end.degree", sector.index = sn)) > 3) {
			xcenter = get.cell.meta.data("xcenter", sector.index = sn, track.index = 2)
			ycenter = get.cell.meta.data("ycenter", sector.index = sn, track.index = 2)
			i_state = as.numeric(gsub("(C|R)_", "", sn))
			circos.text(xcenter, ycenter, i_state, col = "white", font = 2, cex = 0.7, 
				sector.index = sn, track.index = 2, adj = c(0.5, 0.5), niceFacing = TRUE)
			circos.axis(sector.index = sn, track.index = 2, major.tick.percentage = 0.2, labels.away.percentage = 0.2, labels.cex = 0.5)
		}
	}

	if(length(legend_position) == 1) { 
		add_chord_diagram_legend(legend_position, 1:n_states, all_states, state_col[all_states])
	} else if(length(legend_position)== 2) {
		ib = ceiling(n_states/2)
		ind = 1:ib
		add_chord_diagram_legend(legend_position[1], ind, all_states[ind], state_col[all_states][ind])
		ind = (ib+1):n_states
		add_chord_diagram_legend(legend_position[2], ind, all_states[ind], state_col[all_states][ind])
	}

	circos.clear()
}


add_chord_diagram_legend = function(position = c("bottomleft", "bottomright", "topleft", "topright"), 
	index = seq_along(labels), labels, col) {
	
	position = match.arg(position)[1]
	if(length(index) == 0) {
		return(NULL)
	}

	coor = par("usr")
	n = length(labels)
	text_height = strheight("a")
	labels_max_width = max(strwidth(labels))
	legend_width = text_height*(1+0.5) + labels_max_width
	if(position == "bottomleft") {
		x1 = rep(coor[1], n)
		x2 = x1 + text_height
		y1 = coor[3] + (rev(seq_len(n))-1)*1.5*text_height
		y2 = y1 + text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "bottomright") {
		x1 = rep(coor[2] - labels_max_width, n)
		x2 = x1 + text_height
		y1 = coor[3] + (rev(seq_len(n))-1)*1.5*text_height
		y2 = y1 + text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "topleft") {
		x1 = rep(coor[1], n)
		x2 = x1 + text_height
		y1 = coor[4] - (seq_len(n)-1)*1.5*text_height
		y2 = y1 - text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "topright") {
		x1 = rep(coor[2] - labels_max_width, n)
		x2 = x1 + text_height
		y1 = coor[4] - (seq_len(n)-1)*1.5*text_height
		y2 = y1 - text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	}
}

