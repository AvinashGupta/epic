

# == title
# Annotate to gene models
#
# == param
# -gr    a `GenomicRanges::GRanges` object
# -txdb a `GenomicFeatures::TxDb` object. At least the object should contain 'gene_id', 'tx_name'.
#       ('gene_id' and 'tx_name' columns are used to identify genes and transcripts)
# -gene_model type of gene model. By transcripts or by genes
# -species We need this information to find out proper intergenic regions
# -promoters_upstream length of upstream promoter from TSS, pass to `GenomicRanges::promoters`
# -promoters_downstream length of downstream promoter from TSS, pass o `GenomicRanges::promoters`
# -annotation_type pass to `annotate_to_genomic_features`
# -annotation_prefix pass to `annotate_to_genomic_features`
#
# == value
# Following columns are attached to ``gr``:
#
# -nearest_tss the nearest tss (depending on ``gene_model``)
# -dist_to_tss distance to the closest tss (depending on ``gene_model``)
# -nearest_gm the closest gene model (depending on ``gene_model``)
# -dist_to_gm distance to teh closest gene model (depending on ``gene_model``)
# -'prefix_to'_exon percent of the region which is covered by exons or number of exons overlapped to the region
# -'prefix_to'_intron percent of the region which is covered by introns or number of introns overlapped to the region
# -'prefix_to'_promoter percent of the region which is covered by promoters or number of promoters overlapped to the region
# -'prefix_to'_intergenic percent of the region which is covered by intergenic regions or number of intergenic regions overlapped to the region
# -'prefix_to'_fiveUTR percent of the region which is covered by 5'UTRs or number of 5'UTRs overlapped to the region
# -'prefix_to'_threeUTR percent of the region which is covered by 3'UTRs or number of 3'UTRs overlapped to the region
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
annotate_to_gene_models = function(gr, txdb, gene_model =c("tx", "gene"), 
	species = "hg19", promoters_upstream = 2000, promoters_downstream = 200,
    annotation_type = c("percent", "number"), 
    annotation_prefix = "overlap_to_") {
	
	qq.options(LOCAL = TRUE)
	qq.options(code.pattern = "@\\{CODE\\}")

	gene_model = match.arg(gene_model)[1]
	message("preparing annotations for different gene regions.")
	
	### need to make sure 'gene' also has 'transcript' attribute !!! remove it

	message("extracting genes")
	gene = genes(txdb)
	
	# intergenic is gaps between 'gene'
	# we re-define `gene` here to set all strand to '*' (including its levels)
	gene2 = GRanges(seqnames = seqnames(gene),
		            ranges = ranges(gene))
	gene2 = set_proper_seqlengths(gene2, species)

	message("extracting intergenic regions")
	intergenic = gaps(sort(gene2))
	intergenic = intergenic[ strand(intergenic) == "*"]
	
	gene_id = names(gene)

	message("extracting transcripts")
	# note tx also contain 'gene', because for some case, gene contains transcript_id information.
    # need to remove them
	tx = transcripts(txdb)

	# for gene, gene_id == transcript_id
	l = ! mcols(tx)$tx_name %in% gene_id  # real transcript
	# just in case all transcript are all genes
	if(sum(l)) {
		tx = tx[l]
	}
	names(tx) = mcols(tx)$tx_name

	# select a gene model
	if(gene_model == "gene") {
		gm = gene
	} else {
		gm = tx
	}

	message("extracting tss, promoter, exon, intron, 5UTR and 3UTR")
	# tss and promoter are defined according to gene model (gene or transcript)
	tss = promoters(gm, upstream = 0, downstream = 1) # `promoters` is a GRanges method
	promoter = promoters(gm, upstream = promoters_upstream, downstream = promoters_downstream)
	
	# exon, intron, ... are defined from transcript
    # so this would be a little unconsistent here
	#if(gene_model == "tx") {
		exon = exons(txdb)
		intronList = intronsByTranscript(txdb, use.name = TRUE)
		intron = unlist(intronList[ names(intronList) %in% names(tx)])
		fiveUTRList = fiveUTRsByTranscript(txdb, use.name = TRUE)
		fiveUTR = unlist(fiveUTRList[ names(fiveUTRList) %in% names(tx)])
		threeUTRList = threeUTRsByTranscript(txdb, use.name = TRUE)
		threeUTR = unlist(threeUTRList[ names(threeUTRList) %in% names(tx)])
	#} else {
	#	exon = disjointExons(txdb)
	#}
	
    n_gr = length(gr)

	# closest tss
	message("annotating to closest tss")
	nst = as.matrix(nearest(gr, tss, select = "all"))
	m1 = gr[ nst[, 1] ]
	m2 = tss[ nst[, 2] ]
	gi = tapply(names(m2), nst[, 1], paste, collapse = ",")
	suppressWarnings(dst <- tapply(distance(m1, m2), nst[, 1], unique))
    gi2 = rep(NA, n_gr)
	gi2[as.numeric(names(gi))] = gi
    dst2 = rep(NA, n_gr)
	dst2[as.numeric(names(dst))] = dst
	mcols(gr)$nearest_tss = gi2
	mcols(gr)$dist_to_tss = dst2  
	
	# closest transcript
	message(qq("annotating to closest gene model (@{gene_model})"))
	nst = as.matrix(nearest(gr, gm, select = "all"))
	m1 = gr[ nst[, 1] ]
	m2 = gm[ nst[, 2] ]
	gi = tapply(names(m2), nst[, 1], paste, collapse = ",")
	suppressWarnings(dst <- tapply(distance(m1, m2), nst[, 1], unique))
	gi2 = rep(NA, n_gr)
    gi2[as.numeric(names(gi))] = gi
    dst2 = rep(NA, n_gr)
    dst2[as.numeric(names(dst))] = dst
	mcols(gr)$nearest_gm = gi2
	mcols(gr)$dist_to_gm = dst2  
	
	# overlap with exon, intron, intergenic
	gr = annotate_to_genomic_features(gr, gene, "gene", type = annotation_type, prefix = annotation_prefix)
    gr = annotate_to_genomic_features(gr, tx, "tx", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, exon, "exon", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, intron, "intron", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, promoter, "promoter", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, intergenic, "intergenic", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, fiveUTR, "fiveUTR", type = annotation_type, prefix = annotation_prefix)
	gr = annotate_to_genomic_features(gr, threeUTR, "threeUTR", type = annotation_type, prefix = annotation_prefix)
	
	qq.options(LOCAL = FALSE)
	return2(gr)
}

# == title
# Annotate to genomic features
#
# == param
# -gr           a `GenomicRanges::GRanges` object
# -genomic_features a single `GenomicRanges::GRanges` object or a list of `GenomicRanges::GRanges` objects
# -name         names for the genomic features if there is no name in ``genomic_features``
# -type  How to calculate the values for the annotation.
#        ``number`` means numbers of genomic features that each region in ``gr`` overlap; ``percent`` means the 
#        percent of each region in ``gr`` that is overlapped by genomic features
# -prefix prefix for names of the annotation columns
# -... pass to `GenomicRanges::countOverlaps` or `percentOverlaps`
#
# == details
# It adds new columns in ``gr`` which tell you how ``gr`` is overlaped by ``genomic_features``.
#
# Note for the annotation, strand information is ignored.
#
# == value
# A `GenomicRanges::GRanges` with additional columns of annotations.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# df1 = generateRandomBed(nr = 1000)
# df2 = generateRandomBed(nr = 1000)
# df3 = generateRandomBed(nr = 1000)
# gr1 = GRanges(seqnames = df1[[1]], ranges = IRanges(df1[[2]], df1[[3]]))
# gr2 = GRanges(seqnames = df2[[1]], ranges = IRanges(df2[[2]], df2[[3]]))
# gr3 = GRanges(seqnames = df3[[1]], ranges = IRanges(df3[[2]], df3[[3]]))
# annotate_to_genomic_features(gr1, list(gr2 = gr2, gr3 = gr3))
# annotate_to_genomic_features(gr1, list(gr2 = gr2, gr3 = gr3), type = "number", prefix = "#")
annotate_to_genomic_features = function(gr, genomic_features, 
	name = NULL, type = c("percent", "number"), prefix = "overlap_to_", ...) {
	
	op = qq.options()
	qq.options(code.pattern = "@\\{CODE\\}")
	on.exit(qq.options(op))

	# check names
	if(is.null(name)) {
		name = deparse(substitute(genomic_features))
		if(inherits(genomic_features, "list")) {
			if(is.null(names(genomic_features))) {
				name = paste0(name, "_", seq_along(genomic_features))
			} else {
				name = names(genomic_features)
			}
		}
	}

    type = match.arg(type)[1]
	
    if(inherits(genomic_features, "GRanges")) {
        
        message(qq("annotating to @{name}"))

        ostrand = strand(gr)
        strand(gr) = "*"
        strand(genomic_features) = "*"
        if(type == "percent") {
            s2 = percentOverlaps(gr, genomic_features, ...)
        } else if(type == "number") {
            s2 = countOverlaps(gr, genomic_features, ...)
        }

        mcols(gr)[ paste0(prefix, name) ] = s2
        strand(gr) = ostrand

        return(gr)

    }

    # if genomic_features is a list, first annotate to the first one and send the rest recursively
	if(inherits(genomic_features, "list")) {
		if(length(genomic_features) == 0) {
			return(gr)
		}

		if(length(genomic_features) != length(name)) {
			stop("Length of `genomic_features` should be equal to the length of `name`.\n")
		}

		gr = annotate_to_genomic_features(gr, genomic_features[[1]], name = name[1], type = type, prefix = prefix, ...)
		if(length(genomic_features) == 1) {
			return(gr)
		} else {
			genomic_features = genomic_features[-1]
			name = name[-1]
			gr = annotate_to_genomic_features(gr, genomic_features, name = name, type = type, prefix = prefix, ...)
            return(gr)
		}
	}

}

# == title
# Find overlapping genomic regions
#
# == param
# -query a `GenomicRanges::GRanges` object
# -subject a `GenomicRanges::GRanges` object
# -... pass to `GenomicRanges::findOverlaps`
#
# == details
# For every interval in ``query``, it calculates the percent that is covered by ``subject``.
#
# Be careful with ``strand`` in your `GenomicRanges::GRanges` object!!
#
# == value
# A numeric vector which has the same length as ``query``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
# gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
# percentOverlaps(gr1, gr2)
# percentOverlaps(gr2, gr1)
percentOverlaps = function(query, subject, ...) {
	
    subject = reduce(sort(subject))
	mtch = as.matrix(findOverlaps(query, subject, ...))
	
	# intervals in query that overlap with subject
	m1 = query[ mtch[, 1] ]
	# intervals in subject that overlap with query
	m2 = subject[ mtch[, 2] ]

	# overlap between m1 and m2
	m3 = pintersect(m1, m2)

	# tapply on index of `query`
	# it is same as tapply(width(m3), names(m1), sum)
	# but some gr may not have name
	ta = tapply(width(m3), mtch[, 1], sum)
	# width of intervals in query
	w1 = width(query)
	pct = numeric(length(query))
	pct[ as.numeric(names(ta)) ] = ta
	pct = pct / w1 

	# if there is name attribute
	names(pct) = names(query)

	return(pct)
}
