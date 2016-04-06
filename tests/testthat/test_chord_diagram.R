context("test chromatin transition diagram")

if(Sys.getenv("IS_PBS") != "") {

files = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/", pattern = "bed$")
gr_list_1 = lapply(files[1:5], function(f) {
	df = read.table(paste0("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/", f),stringsAsFactors = FALSE)
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), states = df[[4]])
})
gr_list_2 = lapply(files[6:10], function(f) {
	df = read.table(paste0("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/", f),stringsAsFactors = FALSE)
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), states = df[[4]])
})

mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)
chromatin_states_transition_chord_diagram(mat)
chromatin_states_transition_chord_diagram(mat, remove_unchanged_transition = TRUE)
}
