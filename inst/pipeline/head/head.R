
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

    qqcat("[@{chr}] loading /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/rds/@{chr}_roadmap_merged.rds\n")

    bs.fit = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/rds/@{chr}_roadmap_merged.rds"))

    class(bs.fit) = "bs_fit"
    attr(bs.fit, "chr") = chr

    methylation_hooks$obj = bs.fit

    return(invisible(NULL))
}

methylation_hooks$meth = function(bs_fit = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        bs_fit$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$meth[row_index, , drop = FALSE]
    } else {
        bs_fit$meth[row_index, col_index, drop = FALSE]
    }

}

methylation_hooks$raw = function(bs_fit = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        bs_fit$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$meth[row_index, , drop = FALSE]
    } else {
        bs_fit$meth[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$site = function(bs_fit = methylation_hooks$obj, index = NULL) {
    if(is.null(index))
        start(bs_fit$gr)
    else start(bs_fit$gr[index])
}

methylation_hooks$GRanges = function(bs_fit = methylation_hooks$obj) {
    bs_fit$gr
}

methylation_hooks$coverage = function(bs_fit = methylation_hooks$obj,
    row_index = NULL, col_index = NULL) {

    if(is.null(row_index) && is.null(col_index)) {
        bs_fit$cov[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$cov[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$cov[row_index, , drop = FALSE]
    } else {
        bs_fit$cov[row_index, col_index, drop = FALSE]
    }
}

sample_id = dir("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/WGBS_coverage/", pattern = "\\.bed$")
sample_id = gsub("^(E\\d+).*$", "\\1", sample_id)
sample_id = sample_id[!sample_id %in% c("E013", "E054", "E070", "E084", "E085")]
SAMPLE = data.frame(id = sample_id, type = rep("roadmap", length(sample_id)), stringsAsFactors = FALSE)
rownames(SAMPLE) = sample_id
SAMPLE_COLOR = c("roadmap" = "red")

###########################################################################

# wgbs_qcplot("AK100")
# plot_coverage_and_methylation_on_genome("AK100", chromosome = c("chr21", "chr22"))

# global_methylation_distribution(SAMPLE$id)
# global_methylation_distribution(SAMPLE$id, by_chr = TRUE)

gencode_gtf_file = "/icgc/dkfzlsdf/analysis/B080/guz/gencode/gen10.long.gtf"

##############
if(!exists("gene_annotation", envir = .GlobalEnv, inherits = FALSE)) {
    cat("Loading gencode...\n")
    load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gen10_long__transcript_merged.RData")
    map = structure(names(gene_annotation$gtf), names = gsub("\\.\\d+$", "", names(gene_annotation$gtf)))
}
if(!exists("expression", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load expression...\n")
    expression = list()
    expression$count = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/expression/57epigenomes.N.pc.gz", row.names = 1, header = TRUE))
    expression$rpkm = as.matrix(read.table("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/expression/57epigenomes.RPKM.pc.gz", row.names = 1, header = TRUE))
    rownames(expression$count) = map[rownames(expression$count)]
    rownames(expression$rpkm) = map[rownames(expression$rpkm)]
}
if(!exists("txdb", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load txdb...\n")
    txdb = loadDb("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gen10.long.sqlite")
}

gt = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_type")
gt = gt[gt == "protein_coding"]
expr = expression$rpkm[intersect(rownames(expression$rpkm), names(gt)), ]

genes = genes(txdb)
chr = as.vector(seqnames(genes))
names(chr) = names(genes)

chromosome = paste0("chr", c(1:22))
expr = expr[chr[rownames(expr)] %in% chromosome, ]

genome = "hg19"

output_dir = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/image"
RDS_FOLDER = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/rds"

dir.create(output_dir, showWarnings = FALSE)
dir.create(paste0(output_dir, "/gviz"), showWarnings = FALSE)

common_sample_id = intersect(SAMPLE$id, colnames(expr))
SAMPLE = SAMPLE[common_sample_id, ]
expr = expr[, common_sample_id]

if(!exists("GENOMIC_FEATURE_LIST", envir = .GlobalEnv)) {
    # annotate to other regions
    GENOMIC_FEATURE_LIST = c(
        gene                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/gene.bed"),
        exon                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/exon.bed"),
        intron                 = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/intron.bed"),
        tss_2k                 = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/tss_2k.bed"),
        intergenic             = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/intergenic.bed"),
        cgi                    = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/cpgIslandExt.bed"),
        cgi_shore              = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/cgi_shore_2k.bed"),
        dnase                  = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/dnase.bed"),
        enhancer               = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/enhancer.bed"),
        repeats_LINE           = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/repeats_LINE.bed"),
        repeats_SINE           = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/repeats_SINE.bed"),
        tfbs                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/encode_uniform_tfbs_merged_1kb.bed")        # too large for memory
    )
    cat("loading genomic features...\n")
    GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(x) {
        qqcat("reading @{x} as GRanges object...\n")
        df = read.table(x, sep = "\t", stringsAsFactors = FALSE)
        df = df[df[[1]] %in% chromosome, , drop = FALSE]
        makeGRangesFromDataFrame(df[1:3],
                    seqnames.field = "V1",
                    start.field = "V2",
                    end.field = "V3")
    })
}

get_hm_sample = function(hm) {
    hm_sample = dir("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/narrow_peaks", pattern = qq("E\\d+-@{hm}.narrowPeak.gz"))
    hm_sample = gsub(qq("-@{hm}.narrowPeak.gz"), "", hm_sample)
    intersect(SAMPLE$id, hm_sample)
}

get_hm = function(sid, histone_mark) {
    df = read.table(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/roadmap_processed/narrow_peaks/@{sid}-@{histone_mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    df = df[df[[1]] %in% chromosome, ]
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
}

get_hm_list = function(histone_mark) {
    sample = get_hm_sample(histone_mark)

    hm_list = list()
    for(i in seq_along(sample)) {
        cat(sample[i], "\n")
        gr = get_hm(sample[i], histone_mark)
        hm_list[[i]] = gr
    }
    names(hm_list) = sample
    hm_list
}