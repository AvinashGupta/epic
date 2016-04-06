context("test genomic_regions_correlation")

gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
background = GRanges(seqnames = "chr1", ranges = IRanges(start = 7, end = 15))

test_that("test genomic correlation functions", {

	expect_that(genomic_corr_jaccard(gr1, gr2), equals(4/16))
	expect_that(genomic_corr_jaccard(gr1, gr2, background = background), equals(3/8))

	expect_that(genomic_corr_absdist(gr1, gr2), equals(2.5))

	expect_that(genomic_corr_nintersect(gr1, gr2), equals(1))
	expect_that(genomic_corr_pintersect(gr1, gr2), equals(0.4))
	expect_that(genomic_corr_sintersect(gr1, gr2), equals(4))
})

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

### load test data
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/LMR_real_volker.RData")
gr_list = lapply(LMR_list[1:2], makeGRangesFromDataFrameWithFirstThreeColumns)
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/gr_list_1.RData")
genomic_features = lapply(gr_list_1[1:2], makeGRangesFromDataFrameWithFirstThreeColumns)

# general test

gr1 = gr_list[[1]]
gr2 = genomic_features[[1]]
genomicCorr.reldist(gr1, gr2)
genomicCorr.jaccard(gr1, gr2)
genomicCorr.absdist(gr1, gr2)
genomicCorr.absdist(gr1, gr2, method = median)
genomicCorr.absdist(gr1, gr2, trim = 0.1)
genomicCorr.nintersect(gr1, gr2)
genomicCorr.pintersect(gr1, gr2)
genomicCorr.pintersect(gr1, gr2, method = median)
genomicCorr.sintersect(gr1, gr2)

# test if some chromosomes donot exist in the gr2
gr1 = gr1[seqnames(gr1) %in% c("chr1", "chr2")]
gr2 = gr2[seqnames(gr2) %in% c("chr2", "chr3")]
genomicCorr.reldist(gr1, gr2)
genomicCorr.jaccard(gr1, gr2)
genomicCorr.absdist(gr1, gr2)
genomicCorr.nintersect(gr1, gr2)
genomicCorr.pintersect(gr1, gr2)
genomicCorr.sintersect(gr1, gr2)

# test with restrict set
gr1 = gr1[seqnames(gr1) == "chr2"]
gr2 = gr2[seqnames(gr2) == "chr2"]
gr3 = bs.fit@gr
genomicCorr.jaccard(gr1, gr2, restrict = gr3)
genomicCorr.sintersect(gr1, gr2, restrict = gr3)

### test correlation main function
genomic_regions_correlation(gr_list, genomic_features, nperm = 10)
genomic_regions_correlation(gr_list, genomic_features, nperm = 10, chromosome = c("chr1", "chr2"))
genomic_regions_correlation(gr_list, genomic_features, nperm = 10, mc.cores = 2)
genomic_regions_correlation(gr_list, genomic_features, nperm = 10, stat_fun = genomicCorr.reldist)
genomic_regions_correlation(gr_list, genomic_features, nperm = 10, stat_fun = genomicCorr.absdist, trim = 0.1)

### with background
background = systemdf("bedtools random -l 1000 -n 10000 -g /icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/bed/hg19.len | sort -V -k1,1 -k2,2 | bedtools merge -i stdin")
background = makeGRangesFromDataFrameWithFirstThreeColumns(background)
genomic_regions_correlation(gr_list, genomic_features, nperm = 10, background = background)


}