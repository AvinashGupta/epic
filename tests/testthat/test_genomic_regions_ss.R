context("test gr_ss")

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_list_volker.RData")
PMD_list = lapply(PMD_list[1:36], makeGRangesFromDataFrameWithFirstThreeColumns)

cr = common_regions(PMD_list)
cr = common_regions(PMD_list, min_width = 1000)
cr = common_regions(PMD_list, min_width = 10000, min_coverage = 4)

factors = rep(c("T1", "T2", "T3"), each = 12)
res = subgroup_specific_genomic_regions(cr, factors = factors)
res = subgroup_specific_genomic_regions(cr, factors = factors, type = "001")
res = subgroup_specific_genomic_regions(cr, factors = factors, present = 0.8, absent = 0.2)
res = subgroup_specific_genomic_regions(cr, factors = factors, present = function(x) sum(x > 0)/length(x) > 0.5, absent = function(x) sum(x > 0)/length(x) < 0.3)

plot_subgroup_specific_heatmap(res)

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/gr_list_1.RData")
genomic_features = lapply(gr_list_1[1:2], makeGRangesFromDataFrameWithFirstThreeColumns)

plot_subgroup_specific_heatmap(res, genomic_features = genomic_features)

}
