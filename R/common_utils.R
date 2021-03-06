######################################################################################
# this file contains functions that wraps functions with same name in other pacakges
######################################################################################


return2 = function(expr, invisible = FALSE) {
	env = parent.frame()
	
	if(identical(env, .GlobalEnv)) {
		base::return(NULL)
	}
	
	# check on.exit
	if(exists(".on.exit.expression", envir = env)) {
		.on.exit.expression = get(".on.exit.expression", envir = env)
		for(i in seq_along(.on.exit.expression)) {
			eval(.on.exit.expression[[i]], envir = env)
		}
	}

	value = eval(substitute(expr), envir = env)
	
	obj = ls(envir = env)
	
	all_obj = ls(envir = env, all.names = TRUE)
	rm(list = all_obj, envir = env)
	gc(verbose = FALSE)
	
	if(invisible) {
		base::return(invisible(value))
	} else {
		base::return(value)
	}
}

# == title
# Set a counter which represents percent finished in a loop
#
# == param
# -n number of iteration in a loop.
# -fmt format of the message, pass to `base::sprintf`. "\%" will be replaced with the percent string generated by the counter.
# 
# == details
# The counter function should be executed in every iteration in the loop.
#
# == value
# A counter function.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# counter = set_counter(1000)
# for(i in 1:1000) {counter()}
# counter = set_counter(1000, fmt = "processing \%s")
# for(i in 1:1000) {counter()}
set_counter = function(n, fmt = "%s") {

	n = as.integer(n)
	i = 1

	f = function() {
		if(interactive()) {
			pct = round(i/n*100, 1)
			str = paste0(i, "/", n, " (", pct, "%)")
			str = sprintf(fmt, str)

			cat(paste(rep("\b", 200), collapse=""))
			cat(str)
			if(i == n) cat("\n")

			i = i + 1
			assign("i", i, envir = parent.env())
			return(invisible(i))
		}
	}
}

diameter = function(x) {
	r = range(x)
	r[2] - r[1]
}

is.file = function(path) {
	if(length(path) == 1) {
		is.atomic(path) && is.character(path) && file.exists(path)
	} else {
		return(FALSE)
	}
}

sleep = function(time) {
	pb = txtProgressBar(style = 3)
	for(i in seq_len(time)/time) {Sys.sleep(1); setTxtProgressBar(pb, i)}
	close(pb)
}

check_system_command = function(cmd) {
	if(Sys.which(cmd) == "") {
		warning(paste0("Cannot find system command ", cmd, ". Please install it or add the path to PATH.\n"))
	}
}


sort_chr = function(x) {
	y = gsub("^chr(\\d)$", "chr0\\1", x)
	y = gsub("^chr(\\d)_", "chr0\\1_", y)
	x[order(y)]
}

order_chr = function(x) {
	y = gsub("^chr(\\d)$", "chr0\\1", x)
	y = gsub("^chr(\\d)_", "chr0\\1_", y)
	order(y)
}


subset_txdb = function(txdb, chromosome = "chr1") {

	txdump = as.list(txdb)
	txdump$transcripts = txdump$transcripts[txdump$transcripts$tx_chrom %in% chromosome, , drop = FALSE]
	txdump$splicings = txdump$splicings[txdump$splicings$tx_id %in% txdump$transcripts$tx_id, , drop = FALSE]
	txdump$genes = txdump$genes[txdump$genes$tx_id %in% txdump$transcripts$tx_id, , drop = FALSE]
	txdump$chrominfo = txdump$chrominfo[txdump$chrominfo$chrom %in% chromosome, , drop = FALSE]
	txdb2 = do.call(makeTranscriptDb, txdump)

	return(txdb2)
}


# == title
# Find neighbour regions
#
# == param
# -query a `GenomicRanges::GRanges` object
# -reference a `GenomicRanges::GRanges` object 
# -upstream upstream that ``query`` is extended
# -downstream downstream that ``query`` is extended
#
# == details
# With a certain extension of ``query``, this funciton looks for ``reference`` which intersects the extended regions.
#
# == value
# A `GenomicRanges::GRanges` object which contains regions in ``query`` for which the extended regisons are overlapped with ``reference``.
# There are three meta columns added:
#
# -distance: distance from the ``query`` to corresponding ``reference``
# -query_index: index of regions in ``query``
# -reference_index: index of regions in ``reference`` 
#
# Note one ``reference`` can correspond to multiple ``query`` regions.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
# gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
# find_neighbours(gr1, gr2, upstream = 3, downstream = 3)
# find_neighbours(gr1, gr2, upstream = 10, downstream = 10)
find_neighbours = function(query, reference, upstream = 1000, downstream = 1000) {

	query_extended = query
	strd = strand(query)
	start(query_extended) = ifelse(strd == "-", start(query) - downstream, start(query) - upstream)
	end(query_extended) = ifelse(strd == "-", end(query) + upstream, end(query) + downstream)

	mtch = findOverlaps(query_extended, reference)

	mtch = as.matrix(mtch)
	neighbours = query[mtch[, 1]]
	neighbours$distance = distance(query[mtch[,1]], reference[mtch[,2]])

	neighbours$query_index = mtch[, 1]
	neighbours$reference_index = mtch[, 2]
	return(neighbours)
}



top_largest_objects = function() {
	envir = parent.frame()
	s = sapply(ls(envir = envir),function(x) eval(parse(text = qq("object.size(@{x})")), envir = envir))
	print(head(sort(s, decreasing = TRUE)))
}


set_proper_seqlengths = function(gr, species) {
	
	chr_len_df = getChromInfoFromUCSC(species)
	
	chr = as.character(chr_len_df[[1]])
	chr_len = chr_len_df[[2]]
	names(chr_len) = chr
	
	slev = seqlevels(gr)
	slev = chr[chr %in% slev]
	seqlevels(gr) = slev
	slen = chr_len[slev]
	seqlengths(gr) = slen
	return(gr)
}

# == title
# Number of columns which are highly correlated to other columns
#
# == param
# -x a matrix, correlation is calculated by columns
# -abs_cutoff cutoff of absolute correlation
# -size size of blocks
# -mc number of cores
# -... pass to `stats::cor`
#
# == details
# For each column, it looks for number of other columns which correlate with absolute correlation coefficient larger tham ``abs_cutoff``.
# The calculation involves pair-wise correlation of all columns in the matrix.
# When number of columns is huge in a matrix, it is out of ability of R to store such a long vector. This function
# solves this problem by splitting the columns into k blocks and looks at each block sequentially or parallel.
#
# The code is partially adapted from https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/
#
# == value
# A vector that represents how many other columns correlate to current column under the correlation cutoff.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# mat = matrix(rnorm(100000 * 10), ncol = 100000, nrow = 20)
# cor_columns(mat)
# }
# NULL
cor_columns = function (x, abs_cutoff = 0.5, size = 1000, mc = 1, ...) {
    
    split_by_block = function(n, size) {
    	size = min(c(n, size))
    	REST <- n%%size
	    LARGE <- n - REST
	    NBLOCKS <- n%/%size
	    GROUP <- rep(1:NBLOCKS, each = size)
	    if (REST > 0)
	        GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
	    split(1:n, GROUP)
    }

    NCOL <- ncol(x)
    SPLIT = split_by_block(NCOL, size)
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)

    nr = nrow(COMBS)

    count_list = mclapply(split_by_block(nr, floor(nr/mc)), function(ind) {

	    count = matrix(0, nrow = NCOL, ncol = length(abs_cutoff))
	    for (i in ind) {
	        COMB <- COMBS[i, ]

	        qqcat("block @{COMB[1]}(row @{(COMB[1]-1)*size+1}~@{COMB[1]*size})/@{max(COMBS[,1])} and @{COMB[2]}(row @{(COMB[2]-1)*size+1}~@{COMB[2]*size})/@{max(COMBS[,2])}\n")
	        G1 <- SPLIT[[COMB[1]]]
	        G2 <- SPLIT[[COMB[2]]]
	        
	        RES <- cor(x[, G1], x[, G2], ...)
	        RES[is.na(RES)] = 0
	        for(k in seq_along(abs_cutoff)) {
	        	tmp_mat = RES
	        	tmp_mat[abs(tmp_mat) > abs_cutoff[k]] = 1
		        tmp_mat[abs(tmp_mat) < abs_cutoff[k]] = 0

		        count[G1, k] = count[G1, k] + rowSums(tmp_mat)
		        if(COMB[1] != COMB[2]) {
		        	count[G2, k] = count[G2, k] + colSums(tmp_mat)
		        }
		    }
	    }
	    return(count)
	}, mc.cores = mc)

    count = matrix(0, nrow = NCOL, ncol = length(abs_cutoff))
    for(i in seq_along(count_list)) {
    	count = count + count_list[[i]]
    }
	dim(count) = c(NCOL, length(abs_cutoff))
	rownames(count) = colnames(x)
	colnames(count) = abs_cutoff
	return(count)
}


# genes = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::genes(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }

# transcripts = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::transcripts(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }

# exons = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::exons(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }

# intronsByTranscript = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::intronsByTranscript(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
# fiveUTRsByTranscript = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::fiveUTRsByTranscript(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
# threeUTRsByTranscript = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::threeUTRsByTranscript(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
# disjointExons = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::disjointExons(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
# transcriptsBy = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- GenomicFeatures::transcriptsBy(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
# GeneRegionTrack = function(..., max_try = 100, sleep = 60) {
# 	n_try = 0
# 	while(n_try <= max_try) {
# 		oe = try(res <- Gviz::GeneRegionTrack(...))
# 		if(inherits(oe, "try-error")) {
# 			Sys.sleep(sleep)
# 			n_try = n_try + 1
# 		} else {
# 			return(res)
# 		}
# 	}
# }
