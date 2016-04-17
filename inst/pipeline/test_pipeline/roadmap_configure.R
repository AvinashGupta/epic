
## test data from Roadmap project

############################################################
#
methylation_hooks$get_data = function(chr) {

    qqcat("[@{chr}] loading /icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/methylation/@{chr}_roadmap_merged.rds\n")

    obj = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/methylation/@{chr}_roadmap_merged.rds"))
    return(obj)
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


cat("Loading gencode...\n")
load("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10_long_transcript_merged.RData")
map = structure(names(gene_annotation$gtf), names = gsub("\\.\\d+$", "", names(gene_annotation$gtf)))

cat("load expression...\n")
expression = list()
count = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/57epigenomes.N.pc.gz", row.names = 1, header = TRUE))
rpkm = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/57epigenomes.RPKM.pc.gz", row.names = 1, header = TRUE))
rownames(count) = map[rownames(count)]
rownames(rpkm) = map[rownames(rpkm)]
count = count[, intersect(colnames(count), sample_id), drop = FALSE]
rpkm = rpkm[, intersect(colnames(rpkm), sample_id), drop = FALSE]

l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)
expr = rpkm[l, , drop = FALSE]
EXPR = log2(expr + 1)


cat("load txdb...\n")
TXDB = loadDb("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10.long.sqlite")
GTF_FILE = "/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gen10.long.gtf"

GENE_TYPE = "protein_coding"


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

con = pipe("ls /icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/*.narrowPeak.gz")
MARKS = scan(con, what = "character"); close(con)
MARKS = gsub("^.*E\\d+-(.*?)\\.narrowPeak\\.gz$", "\\1", MARKS)
MARKS = sort(unique(MARKS))

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks", pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
}

chipseq_hooks$peak = function(mark, sid) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
}

