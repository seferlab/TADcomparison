library(HiCseg)
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
  stop("Cannot open normalized HiC matices.")
}

print(mat_in)
  
txtdata<-read.table(mat_in)
n_row=nrow(txtdata)
colnames(txtdata)<-c(1:n_row)
mat<-as.matrix(txtdata)
  
Kmax=as.integer(n_row/3)
res=HiCseg_linkC_R(n_row,Kmax,"G",mat,"D")
write.table(res$t_hat,file=domain_out,row.names=F,col.names=F,sep="\t")
