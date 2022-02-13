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

library(adjclust)

mat=read.table(mat_in,header=F)
mat=as.matrix(mat)
colnames(mat)<-c(1:nrow(mat))
rownames(mat)<-c(1:nrow(mat))
mat[is.na(mat)] <- 0

h=floor(nrow(mat)*0.1)
res <- hicClust(mat,h,log = TRUE)
selected.capushe <- select(res)
#TADs=numeric(length(selected.capushe))
TADs<-matrix(0,nrow=length(unique(selected.capushe)),ncol=2)
last=-1

end_i<-length(selected.capushe)
for (i in 1:end_i){
	if (i==1){
		TADs[selected.capushe[i],1]<-i;
	}
	else if(i==end_i){
		TADs[selected.capushe[i],2]<-i;
	}
	else if (last != selected.capushe[i]){
		
		TADs[selected.capushe[i-1],2]<-(i-1);
		TADs[selected.capushe[i],1]<-i;
    }
	last=selected.capushe[i];
}
write.table(TADs, file = domain_out,sep="\t",row.names = FALSE, col.names = FALSE);

