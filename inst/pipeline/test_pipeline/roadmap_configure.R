
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

df = read.table(textConnection(
"eid    group   color
E017    IMR90   #E41A1C
E002    ESC #924965
E008    ESC #924965
E001    ESC #924965
E015    ESC #924965
E014    ESC #924965
E016    ESC #924965
E003    ESC #924965
E024    ESC #924965
E020    iPSC    #69608A
E019    iPSC    #69608A
E018    iPSC    #69608A
E021    iPSC    #69608A
E022    iPSC    #69608A
E007    ES-deriv    #4178AE
E009    ES-deriv    #4178AE
E010    ES-deriv    #4178AE
E013    ES-deriv    #4178AE
E012    ES-deriv    #4178AE
E011    ES-deriv    #4178AE
E004    ES-deriv    #4178AE
E005    ES-deriv    #4178AE
E006    ES-deriv    #4178AE
E062    Blood_&_T-cell  #55A354
E034    Blood_&_T-cell  #55A354
E045    Blood_&_T-cell  #55A354
E033    Blood_&_T-cell  #55A354
E044    Blood_&_T-cell  #55A354
E043    Blood_&_T-cell  #55A354
E039    Blood_&_T-cell  #55A354
E041    Blood_&_T-cell  #55A354
E042    Blood_&_T-cell  #55A354
E040    Blood_&_T-cell  #55A354
E037    Blood_&_T-cell  #55A354
E048    Blood_&_T-cell  #55A354
E038    Blood_&_T-cell  #55A354
E047    Blood_&_T-cell  #55A354
E029    HSC_&_B-cell    #678C69
E031    HSC_&_B-cell    #678C69
E035    HSC_&_B-cell    #678C69
E051    HSC_&_B-cell    #678C69
E050    HSC_&_B-cell    #678C69
E036    HSC_&_B-cell    #678C69
E032    HSC_&_B-cell    #678C69
E046    HSC_&_B-cell    #678C69
E030    HSC_&_B-cell    #678C69
E026    Mesench #B65C73
E049    Mesench #B65C73
E025    Mesench #B65C73
E023    Mesench #B65C73
E052    Myosat  #E67326
E055    Epithelial  #FF9D0C
E056    Epithelial  #FF9D0C
E059    Epithelial  #FF9D0C
E061    Epithelial  #FF9D0C
E057    Epithelial  #FF9D0C
E058    Epithelial  #FF9D0C
E028    Epithelial  #FF9D0C
E027    Epithelial  #FF9D0C
E054    Neurosph    #FFD924
E053    Neurosph    #FFD924
E112    Thymus  #DAB92E
E093    Thymus  #DAB92E
E071    Brain   #C5912B
E074    Brain   #C5912B
E068    Brain   #C5912B
E069    Brain   #C5912B
E072    Brain   #C5912B
E067    Brain   #C5912B
E073    Brain   #C5912B
E070    Brain   #C5912B
E082    Brain   #C5912B
E081    Brain   #C5912B
E063    Adipose #AF5B39
E100    Muscle  #C2655D
E108    Muscle  #C2655D
E107    Muscle  #C2655D
E089    Muscle  #C2655D
E090    Muscle  #C2655D
E083    Heart   #D56F80
E104    Heart   #D56F80
E095    Heart   #D56F80
E105    Heart   #D56F80
E065    Heart   #D56F80
E078    Sm._Muscle  #F182BC
E076    Sm._Muscle  #F182BC
E103    Sm._Muscle  #F182BC
E111    Sm._Muscle  #F182BC
E092    Digestive   #C58DAA
E085    Digestive   #C58DAA
E084    Digestive   #C58DAA
E109    Digestive   #C58DAA
E106    Digestive   #C58DAA
E075    Digestive   #C58DAA
E101    Digestive   #C58DAA
E102    Digestive   #C58DAA
E110    Digestive   #C58DAA
E077    Digestive   #C58DAA
E079    Digestive   #C58DAA
E094    Digestive   #C58DAA
E099    Other   #999999
E086    Other   #999999
E088    Other   #999999
E097    Other   #999999
E087    Other   #999999
E080    Other   #999999
E091    Other   #999999
E066    Other   #999999
E098    Other   #999999
E096    Other   #999999
E113    Other   #999999
"), header = TRUE, row.names = 1, stringsAsFactors = FALSE, comment.char = "")
SAMPLE = data.frame(id = sample_id, class = df[sample_id, "group"], stringsAsFactors = FALSE)
rownames(SAMPLE) = sample_id
COLOR = list(class = structure(unique(df[, "color"]), names = unique(df[, "group"])))


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
MARKS = c("DNase.macs2", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks", pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
	intersect(sample_id, rownames(SAMPLE))
}

chipseq_hooks$peak = function(mark, sid) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
}

