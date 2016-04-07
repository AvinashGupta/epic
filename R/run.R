
# == title
# run pipeline
#
# == details
#
#  Usage: Rscript -e "epic::run()" cmd [options]
#
#  Available cmd:
#
#      chromatin_states_transitions
#      correlated_enriched
#      correlated_regions
#      correlated_regions_downstream
#      correlated_regions_filter
#      correlated_regions_gviz
#      correlated_regions_reduce
#      differential_methylation_in_cgi_and_shore
#      differential_methylation_in_genomic_features
#      general_methylation_distribution
#      methylation_subtype_classification_in_cgi_and_shore
#
# For each cmd, use `Rscript -e "epic::run()" cmd --help` to get help
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
run = function() {

	all_cmds = dir(system.file("pipeline", package = "epic"), pattern = "\\.R$")
	all_cmds = gsub("\\.R$", "", all_cmds)
	
	msg = 'Usage: Rscript -e "epic::run()" cmd [options]\n\nAvailable cmd:\n\n'
	msg = paste0(msg, qq("  @{all_cmds}\n"))

	x = commandArgs(trailingOnly = TRUE)
	R_binary = file.path(R.home("bin"), "R")

	if(length(x) == 0) {
		cat(msg)
		if(interactive()) {
			return(invisible(NULL))
		} else {
			q(save = "no")
		}
	}

	if(!x[1] %in% all_cmds) {
		cat(x, "is not supported.\n\n")
		cat(msg)
		if(interactive()) {
			return(invisible(NULL))
		} else {
			q(save = "no")
		}
	}
	# x = ifelse(grepl(" ", x), paste0("\"", x, "\""), x)

	# cmd = qq("\"@{R_binary}\" --vanilla --slave --args @{paste(x[-1], collapse=' ')} < \"@{system.file('pipeline', package = 'epic')}/@{x[1]}.R\"")
	# cat(cmd, "\n")

	GetoptLong:::source(qq("@{system.file('pipeline', package = 'epic')}/@{x[1]}.R"), argv = paste(x[-1], collapse=' '))

}
