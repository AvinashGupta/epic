
## test data from Roadmap project

############################################################
#
methylation_hooks$set = function(chr) {

    if(!is.null(methylation_hooks$obj)) {
        if(attr(methylation_hooks$obj, "chr") == chr) {
            qqcat("[@{chr}] @{chr} is already set.\n")
            return(invisible(NULL))
        }
    }

    qqcat("[@{chr}] loading /icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/methylation/@{chr}_roadmap_merged.rds\n")

    obj = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/methylation/@{chr}_roadmap_merged.rds"))
    attr(obj, "chr") = chr
    methylation_hooks$obj = obj

    return(invisible(NULL))
}

methylation_hooks$meth = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        obj$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$meth[row_index, , drop = FALSE]
    } else {
        obj$meth[row_index, col_index, drop = FALSE]
    }

}

methylation_hooks$raw = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        obj$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$meth[row_index, , drop = FALSE]
    } else {
        obj$meth[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$site = function(obj = methylation_hooks$obj, index = NULL) {
    if(is.null(index))
        start(obj$gr)
    else start(obj$gr[index])
}

methylation_hooks$GRanges = function(obj = methylation_hooks$obj) {
    obj$gr
}

methylation_hooks$coverage = function(obj = methylation_hooks$obj,
    row_index = NULL, col_index = NULL) {

    if(is.null(row_index) && is.null(col_index)) {
        obj$cov[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$cov[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$cov[row_index, , drop = FALSE]
    } else {
        obj$cov[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$set("chr21")
sample_id = colnames(methylation_hooks$meth(row_index = 1:2))
sample_id = sample_id[!sample_id %in% c("E013", "E054", "E070", "E084", "E085")]
SAMPLE = data.frame(id = sample_id, class = rep("roadmap", length(sample_id)), stringsAsFactors = FALSE)
rownames(SAMPLE) = sample_id
COLOR = list(class = c("roadmap" = "red"))

###########################################################################
gencode_gtf_file = "/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10.long.gtf"

cat("Loading gencode...\n")
load("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10_long_transcript_merged.RData")
map = structure(names(gene_annotation$gtf), names = gsub("\\.\\d+$", "", names(gene_annotation$gtf)))


cat("load expression...\n")
expression = list()
expression$count = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/57epigenomes.N.pc.gz", row.names = 1, header = TRUE))
expression$rpkm = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/57epigenomes.RPKM.pc.gz", row.names = 1, header = TRUE))
rownames(expression$count) = map[rownames(expression$count)]
rownames(expression$rpkm) = map[rownames(expression$rpkm)]

cat("load txdb...\n")
TXDB = loadDb("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10.long.sqlite")

gt = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_type")
gt = gt[gt == "protein_coding"]
EXPR = expression$rpkm[intersect(rownames(expression$rpkm), names(gt)), ]

CHROMOSOME = paste0("chr", c(1:22))

GENOME = "hg19"

OUTPUT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epic_test"

GENOMIC_FEATURE_LIST = c(
    gene                   = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/gene.bed"),
    exon                   = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/exon.bed"),
    intron                 = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/intron.bed"),
    tss_2k                 = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/tss_2k.bed"),
    intergenic             = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/intergenic.bed"),
    cgi                    = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/cpgIslandExt.bed"),
    cgi_shore              = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/cgi_shore_2k.bed"),
    dnase                  = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/dnase.bed"),
    enhancer               = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/enhancer.bed"),
    repeats_LINE           = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/repeats_LINE.bed"),
    repeats_SINE           = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/repeats_SINE.bed"),
    tfbs                   = qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/bed/encode_uniform_tfbs_merged_1kb.bed")        # too large for memory
)
cat("loading genomic features...\n")
GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(x) {
    qqcat("reading @{x} as GRanges object...\n")
    df = read.table(x, sep = "\t", stringsAsFactors = FALSE)
    df = df[df[[1]] %in% CHROMOSOME, , drop = FALSE]
    makeGRangesFromDataFrame(df[1:3],
                seqnames.field = "V1",
                start.field = "V2",
                end.field = "V3")
})


MARKS = c("H3K4me3", "H3K4me1")

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks", pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
    intersect(sample_id, rownames(SAMPLE))
}

chipseq_hooks$peak = function(mark, sid) {
    qqcat("reading peaks: @{sid}, @{mark}")
    df = read.table(qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
    gr[seqnames(gr) %in% CHROMOSOME]
}

