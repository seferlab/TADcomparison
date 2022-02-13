library(CaTCH)
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA,
              help = "normalized HiC matrix"),
  make_option(c("-o", "--output"), type = "character", default = NA,
              help = "output domains")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( is.na(opt$input) | is.na(opt$output)) {
  stop("Missing required parameters. See usage (--help)")
}

mat_in <- opt$input
domain_out <- opt$output

if ( !file.exists(mat_in) ) {
  stop("Cannot open inputfile.")
}

#fileinput=system.file(mat_in,package="CaTCH")
res=domain.call(mat_in)
result=res$clusters
result=result[which(result$RI > 0.649),]
result=result[which(result$RI < 0.651),]
write.table(result[,c('start','end')], file=domain_out, append = FALSE, row.names = FALSE, col.names = FALSE)
