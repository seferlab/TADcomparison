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

ndomains=res$ncluster[,'ndomains']
diff_len=length(ndomains)-1
ndomains_diff=array(0,dim=c(diff_len))
for(i in 1:diff_len){
ndomains_diff[i]<-(ndomains[i]-ndomains[i+1])
}
ndomains_diff_2=ndomains_diff[20:(diff_len-20)]
max_n=which(ndomains_diff_2==max(ndomains_diff_2),arr.ind=TRUE)
max_n=max_n+20-1
final_RI=res$ncluster[max_n,'RI']

result=res$clusters
result=result[which(result$RI >= final_RI ),]
result=result[which(result$RI <= final_RI),]
write.table(result[,c('start','end')], file=domain_out, append = FALSE, row.names = FALSE, col.names = FALSE)

