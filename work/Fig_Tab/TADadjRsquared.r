library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA,
              help = "Hi-C matrix path"),
  make_option(c("-d", "--dataname"), type = "character", default = NA,
              help = "data name"),
  make_option(c("-c", "--cntchr"), type = "numeric", default = NA,
              help = "number of chromosomes"),
  make_option(c("-o", "--output"), type = "character", default = NA,
              help = "output name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( is.na(opt$input) | is.na(opt$dataname) | is.na(opt$cntchr) | is.na(opt$output)) {
  stop("Missing required parameters. See usage (--help)")
}


source('ari.R')
res= 50000
dis=200
getAUC<-function(ari)
{	AUC=c()
	for (i in 1:dim(ari$pr)[1]){
		o = 0
		x <- c(0,as.numeric(ari$rc[i,]),1)
		y <- c(1,as.numeric(ari$pr[i,]),0)
		for(i in 1:(length(x)-1)){o = o+ (x[i+1]-x[i])*(y[i+1]+y[i])/2}
		AUC <- c(AUC,o)
                #id <- order(x)
                #AUC <- c(AUC,sum(diff(x[id])*rollmean(y[id],2)))
	}
	AUC[is.na(AUC)] <- 0
	return(as.numeric(AUC))
}

ArmatusAUC=rep(0,dis+1)
ArrowheadAUC=rep(0,dis+1)
CaTCHAUC=rep(0,dis+1)
CHACAUC=rep(0,dis+1)
CHDFAUC=rep(0,dis+1)
ClusterTADAUC=rep(0,dis+1)
deDocAUC=rep(0,dis+1)
DIAUC=rep(0,dis+1)
EASTAUC=rep(0,dis+1)
GMAPAUC=rep(0,dis+1)
HiCExplorerAUC=rep(0,dis+1)
HiCsegAUC=rep(0,dis+1)
IC_FinderAUC=rep(0,dis+1)
InsulationScoreAUC=rep(0,dis+1)
MatryoshkaAUC=rep(0,dis+1)
MrTADFinderAUC=rep(0,dis+1)
MSTDAUC=rep(0,dis+1)
OnTADAUC=rep(0,dis+1)
SpectralAUC=rep(0,dis+1)
SpectralTADAUC=rep(0,dis+1)
TADBDAUC=rep(0,dis+1)
TADbitAUC=rep(0,dis+1)
TADtreeAUC=rep(0,dis+1)
TopDomAUC=rep(0,dis+1)





dname=as.character(opt$dataname)
current_method=''

n=0
n_arrowhead=0
n_TADbit=0
chrs=c(1:opt$cntchr)
chrs=as.character(chrs)
chrs=c(chrs,'X')

for (chrnum in chrs){
	n = n + 1
	c_file=gsub("chr1",paste('chr',chrnum,sep=''),opt$input)
	hmat <- read.table(c_file)
	hmat <- data.matrix(hmat)
	max_bound<-nrow(hmat)
	
	current_method='Armatus'
	ArmatusTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	ArmatusTAD[ArmatusTAD>max_bound]<-max_bound
	Armatusari = getadjr2(hmat,ArmatusTAD,dis)
	ArmatusAUC = ArmatusAUC+ Armatusari$adjr2[,1]
	
	current_method='Arrowhead'
	mat_in=paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep='')
	if ( file.exists(mat_in) ) {
	n_arrowhead =n_arrowhead+1
	ArrowheadTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	ArrowheadTAD[ArrowheadTAD>max_bound]<-max_bound
	Arrowheadari = getadjr2(hmat,ArrowheadTAD,dis)
	ArrowheadAUC = ArrowheadAUC+ Arrowheadari$adjr2[,1]
	}
	
	current_method='CaTCH'
	CaTCHTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	CaTCHTAD[CaTCHTAD>max_bound]<-max_bound
	CaTCHari = getadjr2(hmat,CaTCHTAD,dis)
	CaTCHAUC = CaTCHAUC+ CaTCHari$adjr2[,1]
	
	current_method='CHAC'
	CHACTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	CHACTAD[CHACTAD>max_bound]<-max_bound
	CHACari = getadjr2(hmat,CHACTAD,dis)
	CHACAUC = CHACAUC+ CHACari$adjr2[,1]
	
	current_method='CHDF'
	CHDFTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	CHDFTAD[CHDFTAD>max_bound]<-max_bound
	CHDFari = getadjr2(hmat,CHDFTAD,dis)
	CHDFAUC = CHDFAUC+ CHDFari$adjr2[,1]
	
	current_method='ClusterTAD'
	ClusterTADTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	ClusterTADTAD[ClusterTADTAD>max_bound]<-max_bound
	ClusterTADari = getadjr2(hmat,ClusterTADTAD,dis)
	ClusterTADAUC = ClusterTADAUC+ ClusterTADari$adjr2[,1]
	
	current_method='deDoc'
	deDocTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	deDocTAD[deDocTAD>max_bound]<-max_bound
	deDocari = getadjr2(hmat,deDocTAD,dis)
	deDocAUC = deDocAUC+ deDocari$adjr2[,1]
	
	current_method='DI'
	DITAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	DITAD[DITAD>max_bound]<-max_bound
	DIari = getadjr2(hmat,DITAD,dis)
	DIAUC = DIAUC+ DIari$adjr2[,1]
	
	current_method='EAST'
	EASTTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	EASTTAD[EASTTAD>max_bound]<-max_bound
	EASTari = getadjr2(hmat,EASTTAD,dis)
	EASTAUC = EASTAUC+ EASTari$adjr2[,1]
	
	current_method='GMAP'
	GMAPTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	GMAPTAD[GMAPTAD>max_bound]<-max_bound
	GMAPari = getadjr2(hmat,GMAPTAD,dis)
	GMAPAUC = GMAPAUC+ GMAPari$adjr2[,1]
	
	current_method='HiCExplorer'
	HiCExplorerTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	HiCExplorerTAD[HiCExplorerTAD>max_bound]<-max_bound
	HiCExplorerari = getadjr2(hmat,HiCExplorerTAD,dis)
	HiCExplorerAUC = HiCExplorerAUC+ HiCExplorerari$adjr2[,1]
	
	current_method='HiCseg'
	HiCsegTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	HiCsegTAD[HiCsegTAD>max_bound]<-max_bound
	HiCsegari = getadjr2(hmat,HiCsegTAD,dis)
	HiCsegAUC = HiCsegAUC+ HiCsegari$adjr2[,1]
	
	current_method='IC-Finder'
	IC_FinderTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	IC_FinderTAD[IC_FinderTAD>max_bound]<-max_bound
	IC_Finderari = getadjr2(hmat,IC_FinderTAD,dis)
	IC_FinderAUC = IC_FinderAUC+ IC_Finderari$adjr2[,1]
	
	current_method='InsulationScore'
	InsulationScoreTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	InsulationScoreTAD[InsulationScoreTAD>max_bound]<-max_bound
	InsulationScoreari = getadjr2(hmat,InsulationScoreTAD,dis)
	InsulationScoreAUC = InsulationScoreAUC+ InsulationScoreari$adjr2[,1]
	
	current_method='Matryoshka'
	MatryoshkaTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	MatryoshkaTAD[MatryoshkaTAD>max_bound]<-max_bound
	Matryoshkaari = getadjr2(hmat,MatryoshkaTAD,dis)
	MatryoshkaAUC = MatryoshkaAUC+ Matryoshkaari$adjr2[,1]
	
	current_method='MrTADFinder'
	MrTADFinderTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	MrTADFinderTAD[MrTADFinderTAD>max_bound]<-max_bound
	MrTADFinderari = getadjr2(hmat,MrTADFinderTAD,dis)
	MrTADFinderAUC = MrTADFinderAUC+ MrTADFinderari$adjr2[,1]
	
	current_method='MSTD'
	MSTDTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	MSTDTAD[MSTDTAD>max_bound]<-max_bound
	MSTDari = getadjr2(hmat,MSTDTAD,dis)
	MSTDAUC = MSTDAUC+ MSTDari$adjr2[,1]
	
	current_method='OnTAD'
	OnTADTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	OnTADTAD[OnTADTAD>max_bound]<-max_bound
	OnTADari = getadjr2(hmat,OnTADTAD,dis)
	OnTADAUC = OnTADAUC+ OnTADari$adjr2[,1]
	
	current_method='Spectral'
	Spectral_TAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	Spectral_TAD[Spectral_TAD>max_bound]<-max_bound
	Spectralari = getadjr2(hmat,Spectral_TAD,dis)
	SpectralAUC = SpectralAUC+ Spectralari$adjr2[,1]
	
	current_method='SpectralTAD'
	SpectralTADTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	SpectralTADTAD[SpectralTADTAD>max_bound]<-max_bound
	SpectralTADari = getadjr2(hmat,SpectralTADTAD,dis)
	SpectralTADAUC = SpectralTADAUC+ SpectralTADari$adjr2[,1]
	
	current_method='TADBD'
	TADBDTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	TADBDTAD[TADBDTAD>max_bound]<-max_bound
	TADBDari = getadjr2(hmat,TADBDTAD,dis)
	TADBDAUC = TADBDAUC+ TADBDari$adjr2[,1]
	
	current_method='TADbit'
	#mat_in=paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep='')
	#if ( file.exists(mat_in) ) {
	#n_TADbit =n_TADbit+1
	TADbitTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	TADbitTAD[TADbitTAD>max_bound]<-max_bound
	TADbitari = getadjr2(hmat,TADbitTAD,dis)
	TADbitAUC = TADbitAUC+ TADbitari$adjr2[,1]
	#}
	
	current_method='TADtree'
	TADtreeTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	TADtreeTAD[TADtreeTAD>max_bound]<-max_bound
	TADtreeari = getadjr2(hmat,TADtreeTAD,dis)
	TADtreeAUC = TADtreeAUC+ TADtreeari$adjr2[,1]
	
	current_method='TopDom'
	TopDomTAD=read.table(paste('../all_TADs/bin/',current_method,'/',dname,'_',current_method,'.chr',chrnum,sep=''),sep='\t')
	TopDomTAD[TopDomTAD>max_bound]<-max_bound
	TopDomari = getadjr2(hmat,TopDomTAD,dis)
	TopDomAUC = TopDomAUC+ TopDomari$adjr2[,1]
	
	

}

ArrowheadAUC = ArrowheadAUC*n/n_arrowhead
#if (n_TADbit!=0){
#TADbitAUC = TADbitAUC*n/n_TADbit
#}

out = rbind(ArmatusAUC,ArrowheadAUC,CaTCHAUC,CHACAUC,CHDFAUC,ClusterTADAUC,deDocAUC,DIAUC,HiCsegAUC,IC_FinderAUC,InsulationScoreAUC,EASTAUC,GMAPAUC,HiCExplorerAUC,MatryoshkaAUC,MrTADFinderAUC,MSTDAUC,OnTADAUC,SpectralAUC,SpectralTADAUC,TADBDAUC,TADbitAUC,TADtreeAUC,TopDomAUC)/n
write.table(out,file=opt$output,row.names = F, col.names = F, quote = F,sep = '\t')

