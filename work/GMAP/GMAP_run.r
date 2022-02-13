library(rGMAP)
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA, help = "normalized HiC matrix"),
  make_option(c("-o", "--output"), type = "character", default = NA, help = "output domains"),
  make_option(c("-r", "--resolution"), type = "numeric", default = NA, help = "resolution")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( is.na(opt$input) | is.na(opt$output)) {
  stop("Missing required parameters. See usage (--help)")
}

mat_in <- opt$input
domain_out <- opt$output
reso <- as.numeric(opt$resolution)

if ( !file.exists(mat_in) ) {
  stop("Cannot open inputfile.")
}

datax<-read.table(mat_in)
#res = rGMAP(datax, index_obj = NULL, resl = 50000, dom_order = 2)
res = rGMAP(datax, resl = reso)
#print(res$hierTads)
#print(res$tads)
result=res$hierTads
result=result[,c('start','end')]
write.table(result, file = domain_out,sep = '\t',row.names=FALSE,col.names=FALSE)


