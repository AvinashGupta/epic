context("test basic_genomic_regions_stat")

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_list_volker.RData")
gr_list = lapply(PMD_list, makeGRangesFromDataFrameWithFirstThreeColumns)

basic_genomic_regions_stat(gr_list, type = "proportion")
basic_genomic_regions_stat(gr_list, type = "number")
basic_genomic_regions_stat(gr_list, type = "median_width")
basic_genomic_regions_stat(gr_list, type = "proportion", by_chr = TRUE)
basic_genomic_regions_stat(gr_list, type = "number", by_chr = TRUE)
basic_genomic_regions_stat(gr_list, type = "median_width", by_chr = TRUE)


}