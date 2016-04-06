context("test genomic annotations")

gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))

test_that("test percentOverlaps", {
	expect_that(percentOverlaps(gr1, gr2), equals(c(0, 4/7)))
	expect_that(percentOverlaps(gr2, gr1), equals(c(0, 4/8)))
})

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_list_volker.RData")
gr = makeGRangesFromDataFrameWithFirstThreeColumns(PMD_list[[1]])
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/gr_list_1.RData")
genomic_features = lapply(gr_list_1[1:2], makeGRangesFromDataFrameWithFirstThreeColumns)

###
annotate_to_genomic_features(gr, genomic_features[[1]])
annotate_to_genomic_features(gr, genomic_features[[1]], name = "gene")
annotate_to_genomic_features(gr, genomic_features[[1]], name = "gene", type = "number")

annotate_to_genomic_features(gr, genomic_features)
annotate_to_genomic_features(gr, genomic_features, name = c("G", "E"))
annotate_to_genomic_features(gr, genomic_features, prefix = "")

### build a transcriptDb object
library(GenomicFeatures)
txdb = loadDb("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode.v17.sqlite")

annotate_to_gene_models(gr, txdb, gene_model = "tx")
annotate_to_gene_models(gr, txdb, gene_model = "gene")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_type = "number")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_prefix = "")

}