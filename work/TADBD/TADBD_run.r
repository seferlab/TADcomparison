library(TADBD)
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA, help = "normalized HiC matrix"),
  make_option(c("-o", "--output"), type = "character", default = NA, help = "output domains")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( is.na(opt$input) | is.na(opt$output)) {
  stop("Missing required parameters. See usage (--help)")
}
mat_in=opt$input
domain_out=opt$output
#species <- "hg19"
#chr <- "chr18"
#resolution <- 50000
#options(scipen = 999)
#data(hicdata)
#hicmat <- DataLoad(hicdata, bsparse = FALSE, species, chr, resolution)
hicmat=read.table(mat_in)
hicmat=as.matrix(hicmat)
df_result=TADBD(hicmat)
#df_result$domains
write.table(df_result$domains, file = domain_out,sep = '\t',row.names=FALSE,col.names=FALSE)
