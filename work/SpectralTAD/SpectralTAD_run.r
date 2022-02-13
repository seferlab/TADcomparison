library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA, help = "normalized HiC matrix"),
  make_option(c("-o", "--output"), type = "character", default = NA, help = "output domains"),
  make_option(c("-r", "--resolution"), type = "numeric", default = NA, help = "resolution"),
  make_option(c("-c", "--chromosome"), type = "character", default = NA, help = "chromosome name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( is.na(opt$input) | is.na(opt$output)) {
  stop("Missing required parameters. See usage (--help)")
}

mat_in <- opt$input
domain_out <- opt$output

reso <- as.numeric(opt$resolution)
chrom <- opt$chromosome

if ( !file.exists(mat_in) ) {
  stop("Cannot open inputfile.")
}

library(SpectralTAD)
mat<-read.table(mat_in)
mat<-data.matrix(mat)
maxn<-ncol(mat)-1
colnames(mat)=rownames(mat)=c(0:maxn)*reso
tads = SpectralTAD(mat, chr = chrom, levels = 2, qual_filter = FALSE, resolution=reso)
write.table(tads$Level_1, file=domain_out, append = FALSE, row.names = FALSE, col.names = FALSE)
write.table(tads$Level_2, file=domain_out, append = TRUE, row.names = FALSE, col.names = FALSE)
