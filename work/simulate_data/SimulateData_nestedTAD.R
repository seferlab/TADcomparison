

###################################################
###################################################
#
# This script generates simulated Hi-C data matrices as per the procedure used and described in
# Forcato et al. Nature Methods 2017 - Comparison of computational methods for Hi-C data analysis.
# PMID: 28604721 DOI: 10.1038/nmeth.4325 
# 
# CREDIT - AKCNOWLEDGEMENTS - DISCLAIMER
#
# The procedure adopted here is an extension of the original Hi-C data simulation procedure 
# proposed by diffHiC developer (Aaron Lun) and presented in 
# Lun AT, Smyth GK. diffHic: a Bioconductor package to detect differential
# genomic interactions in Hi-C data. BMC Bioinformatics. 2015 Aug 19;16:258. doi:
# 10.1186/s12859-015-0683-0. PubMed PMID: 26283514
# (see in particular the procedure described in suplementary information page 5 of this article)
#
# We thank Aaron Lun (University of Cambridge) for sharing the code used to simulate Hi-C data in the diffHic article. 
# Some portions of Aaron Lun's script original code are included in this script and marked with comments where
# they are more similar to the original code.
# More in general the structure of the simulation procedure 
# The portions of this script that most retain the original code are marked with comments for proper credit.
# However, this does not imply that Dr Lun is in any way responsible for the code 
#
###################################################
###################################################


## Version 7
# adding a control parameter to force a capping of nestedTadSizes

## Version 6
# adding the code and parameters to simulate multiple levels of nested TADs 
 

### NOTE FOR COMPATIBILITY WITH OTHER PACKAGES
# In this version we use R 3.2.0 version as diffHiC version 1.0.0 has been used to analyze the data 
# This detail is important as output files are also saved in diffHic format,
# and the data structures are different in older vs newer diffHiC versions
# as such we must use the same diffHic version to save the h5 data format
# R version = R 3.2.0
# diffHiC version = should be 1.0.0
#
# However, right now  diffHiC version 1.0.0 is no more available.
# Instead the version 1.0.1 is available in bioconductor (BioC version 3.1)
# There's also a problem of conflicting versions in bioconductor because R.3.2 is associated to both BioC3.1 and BioC3.2. 
# When installing packages with biocLite() there is apparently no way to force the installation of BioC3.1. instead than Bioc3.2 when any R 3.2.x version is present
# Thus I had to install manually diffHic version 1.0.1 (corresponding to Bioc3.1)
# 
# Nevertheless, the diffHiC version installed on the system is only relevant for the command (savePairs)
# used to format output object in the diffHiC format. In case the output is only saved in text format, this is not relevant.
# In case a different (newer) diffHiC version is used to save the output data, then you just have to make sure 
# that the function savePairs is called with the right (newer) set of parameters,
# that will be largely the same, but may have slightly different names
#
###

 



################################################################################
################################################################################
###
### Other notes on parameters with examples
###
   
# 
# ## label for simulated chromosome
# simChromName<-"ChrS"

# ## genomics bin size for the simulate matrix
# BINSIZE<-40000

# # max size of TADs (in number of bins)
# max.tad <- 50  
# # minimum size of TADs (in bins
# min.tad<-3

# # number of  contiguous, non-overlapping TADs.
# #T ARGET.chrsize<-180.92*1e6 # same as chr5 in hg19
# #  to simulate a chromosome as big as human chr5 the number of tads (ranging from 3 to 50 bins) should be around 170
# # > mean(3:50)
# # [1] 26.5
# # > 180.92*1e6
# # [1] 180920000
# # > 180.92*1e6/40000/26.5
# # [1] 170.6792
	
# # signal decay in TAD
# base.tad.signal<-28
# interaction.decay.rate<-0.69

## random noise: here we define the random noise as a percentage of data points to which noise is addedd. 
## the percentage is (number of data points to which noise is added)/(number of data points in upper triangular matrix)
## note: the actual number of data points to which noise is added will actually be smaller,
## as some data points may be sampled twice or more to add them noise twice or more, depending on the sampling strategy
# in the original procedure (By Aaron Lun) settings, given upper triangle matrix size = ((3112^2)-3112)/2 = 4840716
# and the ratio > 200000/4840716 = [1] 0.0413162 --> 4% of sampling with noise
# percent_of_sampling_points_nonSpecificsignal<-0.04

# # estimated_binsperchromsome<-round((TARGET.chrsize/BINSIZE))
# # estimated_matrix_size_upperTriangle<-(((estimated_binsperchromsome)^2-(estimated_binsperchromsome))/2)
# # 
# # SampledPointNonSpecificSignal<-round(percent_of_sampling_points_nonSpecificsignal * estimated_matrix_size_upperTriangle) ## %%% changed to make this proportional to changing parameters
## we changed this compared to the original procedure to make this noise signal level proportional to changing parameters of base TAD signal
# # base.nonspecific.signal<-round((5/80)*base.tad.signal)


# ## interactions
# ## original was 100 --> over matrix size 100/4840716 =  2.06581e-05 (i.e. 2 every 100k data points)
## we changed to make this proportional to changing parameters
# n.interactions <-  round(2e-5*estimated_matrix_size_upperTriangle)
# # > round(2e-5*estimated_matrix_size_upperTriangle)
# # [1] 205

# here I'm using just half of the interactions used in the simulation for diffHiC (they were 200)
# because in those simulations they wanted to test differential interaction calls and for that
# half of the 200 interactions were simulated in one group of samples
# whereas the other half was simulated in the other group of samples


## max distance from diagonal for the interactions
# (no more than 10 bins beyond TAD max size as there are few interactions very far from diagonal in real data)
## this was changed to make this proportional to changing parameters 
### not going beyond diagonal #40... we may consider also changing this parameter
# max.range <- max.tad + round(max.tad/3) 
# # signal at interaction
## %%% changed to make this proportional to changing parameters
# base.interaction.signal<-base.tad.signal*2 

# ## read counts simulation
# we kept the original parameter for this
# dispersion <- 0.01 




#####################################################
#####################################################
## accessory functions

##
## Function from original "comparator.R" script by Aaron Lun
##
altmatch <- function(x1, x2, y1, y2) {
	# Quick dirty matcher.
	pair.x <- paste(x1, x2, sep=".")
	pair.y <- paste(y1, y2, sep=".")
	matched <- match(pair.x, pair.y)
	keep <- !is.na(matched)
	return(list(unique1=which(!keep),
				unique2=(1:length(pair.y))[-matched[keep]],
				match1=which(keep),
				match2=matched[keep]))
}


#####################################################
#####################################################
##
## the following functions were provided by Aaron Lun "xscss" R package
## --> from "simulateHiC.R" file inside xscss
## but were modified to implement additional options, such as the definition of nested TADs

# ## see for comparison the original spawnTADs function here
# spawnTADs <- function(ntads, min.tad=10, max.tad=30) 
# # This function spawns counts for TADs. It returns a list of matrices
# # indicating the pairs of restriction site indices in each TAD.
# {
# 	ntads <- as.integer(ntads)
# 	mysize <- as.integer(runif(ntads, min.tad, max.tad))
# 	tads <- list()
# 	last.launch <- 0L
# 	for (i in 1:ntads) { 
# 		current <- last.launch + mysize[i]:1   ## %%% catch max min of this to return the TADs start end positions
# 		tads[[i]] <- t(combn(current, 2L))
# 		colnames(tads[[i]]) <- c("anchor", "target")
# 		last.launch <- last.launch + mysize[i]
# 	}
# 	tads
# }


## this "dual output" is to return both TAD borders and TAD intervals 
# This function spawns counts for TADs. It returns a list of matrices
# indicating the pairs of restriction site indices in each TAD.
# with nesting is to have also nested TADs
spawnTADs_withNesting_dualOutput <- function(ntads, min.tad=10, max.tad=30, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL) 
{

    #check parameters
    if (!is.null(nestedSizeCap)) {
        if (nestedSizeCap<=max.tad) {
        stop("nestedSizeCap should be larger than max.tad")
        }
    }

	ntads <- as.integer(ntads)
	mysize <- as.integer(runif(ntads, min.tad, max.tad))
	tads <- list()
	tads_intervals<-NULL

	# copy object for subsequent manipulations
	mysize_ref<-mysize
	
	for (nestingLevel in 1:nestedLevels) {
	    if (nestingLevel==1) {
	        borderSkippingProbs=0
	    } else {
    	    borderSkippingProbs=nestedProb
	    }


        # sampling TADs to be removed with probability nestedProb
        # the selected TADs will be merged to the following one, as such we sample up to (mysize length -1)
        current_BordersRemoval<-sort(sample(1:(length(mysize_ref)-1), size=round(borderSkippingProbs*(length(mysize_ref)-1))))

        current_mysize<-NULL
        ## TADs selected for removal will just be merged to the subsequent one
        for (imerge in 1:length(mysize_ref)) {
            if (imerge %in% current_BordersRemoval) {
                # cap the maximum size of nested TADs
                if ( !is.null(nestedSizeCap)) {
                    if ((mysize_ref[imerge]+mysize_ref[(imerge+1)]) > nestedSizeCap) {
                    # skip the merging it the merged TAD would be larger than the cap size
                    current_mysize<-c(current_mysize, mysize_ref[imerge])
                    } else {
                    mysize_ref[(imerge+1)]<-(mysize_ref[imerge]+mysize_ref[(imerge+1)])
                    }
                } else {
                mysize_ref[(imerge+1)]<-(mysize_ref[imerge]+mysize_ref[(imerge+1)])
                }
            } else {
                current_mysize<-c(current_mysize, mysize_ref[imerge])
            }
        }
        mysize_ref<-current_mysize

	
# > current_BordersRemoval
#  [1]  49  44  51  33 148 134  13 111 169  99  35  39  72 108  65  90   2  78 121
# [20]  31   4 125 124 150  16 145 113  26  54  95  83 105  58 117 137 163 155  56
# [39]  25  38  85 154  43 106  73 118 153  70  74  14 133
# > sort(current_BordersRemoval)
#  [1]   2   4  13  14  16  25  26  31  33  35  38  39  43  44  49  51  54  56  58
# [20]  65  70  72  73  74  78  83  85  90  95  99 105 106 108 111 113 117 118 121
# [39] 124 125 133 134 137 145 148 150 153 154 155 163 169
# 
	
    	last.launch <- 0L
   	    pre_existing_tads<-length(tads)
        for (i in 1:length(mysize_ref)) { 
            current <- last.launch + mysize_ref[i]:1   ## %%% catch max min of this to return the TADs start end positions
            tads_intervals<-rbind(tads_intervals, range(current))
            tads[[i+pre_existing_tads]] <- t(combn(current, 2L))
            colnames(tads[[i+pre_existing_tads]]) <- c("anchor", "target")
            last.launch <- last.launch + mysize_ref[i]
        }
    }	
	
	
	colnames(tads_intervals)<-c("startbin", "endbin")
	return(list(tads_list=tads, tads_intervals=tads_intervals))
}

## this "dual output" is to return both TAD borders and TAD intervals 
# This function spawns counts for TADs. It returns a list of matrices
# indicating the pairs of restriction site indices in each TAD.
spawnTADs_dualOutput <- function(ntads, min.tad=10, max.tad=30) 
{
	ntads <- as.integer(ntads)
	mysize <- as.integer(runif(ntads, min.tad, max.tad))
	tads <- list()
	last.launch <- 0L
	tads_intervals<-NULL
	for (i in 1:ntads) { 
		current <- last.launch + mysize[i]:1   ## %%% catch max min of this to return the TADs start end positions
		tads_intervals<-rbind(tads_intervals, range(current))
		tads[[i]] <- t(combn(current, 2L))
		colnames(tads[[i]]) <- c("anchor", "target")
		last.launch <- last.launch + mysize[i]
	}
	colnames(tads_intervals)<-c("startbin", "endbin")
	return(list(tads_list=tads, tads_intervals=tads_intervals))
}


### This is the same as original Aaron Lun function
getMeans <- function(mat, base=10, decay=0.5, prior=1, rate=1) 
# This makes mean values for the requested matrices that are coming in,
# based on the log-log relationship observed between mean and distance.
{
	distance <- mat[,1]-mat[,2] + prior
	true.mean <- base*distance^(-decay)
	if (is.na(rate)) { return(true.mean) } else {
	return(rgamma(nrow(mat), shape=true.mean*rate, rate=rate)) } ## rate is set to NA below thus it is returning the treu mean
}


### This is the same as original Aaron Lun function
spawnDiagon <- function(len) 
# Adds coordinates for the diagonals. 
{
	len <- as.integer(len)
	cbind(anchor=1:len, target=1:len)
}

### This is the same as original Aaron Lun function
# but with more comments added
spawnLigations2 <- function(nspec, len, diag=FALSE, prior=1, decay=0.5) 
# This adds ligation events where the weighting of the selected ligation events is
# skewed towards the diagonal. This avoids putting too much work into the rest of the
# contact matrix that is unlikely to be interesting.
{ 
	len <- as.integer(len) ## %%% number of rows (diagonals) in the matrix
	all.dists <- (len-1):ifelse(diag, 0, 1) ## %%%  define whether to include or not the diagonal x=y (no in the script below)
	n.entries <- len - all.dists ## %%%  define how many data points are available at each considered diagonal (distance)
	chosen.dists <- sample(all.dists, nspec, p=n.entries*(all.dists+prior)^(-decay), replace=TRUE)  ## %%% sample the data points with replacement and probability proportional to how many data points are available at each diagonal and decreasing with power law (x-y+p)^c

	chosen <- split(1:nspec, chosen.dists) ## %%% group the selected data points by diagonals and format as list
	anchor <- target <- integer(nspec) ## %%% define empty vector
	for (dist in names(chosen)) { ## %%% for each diagonal
		relevants <- chosen[[dist]] ## %%% how many data points needs to be samples (result from sample above)
		allowable.anchors <- (as.integer(dist)+1L):len ## %%%  how many possible data points in that diagonal?
		cur.anchors <- sample(allowable.anchors, length(relevants), replace=TRUE) ## %%% sample with replacement among possible data points of that diagonals, sample length(relevants) elements which is the result of the previous sample function call
		cur.targets <- cur.anchors - as.integer(dist)  ## %%% define the y coordinates by subtracting the diagonals
		anchor[relevants] <- cur.anchors ## %%%  save value sin vector
		target[relevants] <- cur.targets ## %%%  save value sin vector 
	}

	o <- order(anchor, target) ## %%% sort and return
	cbind(anchor=anchor[o], target=target[o])
}


### This is the same as original Aaron Lun function
spawnLigations <- function(nspec, max.range, len, diag=TRUE) 
# This adds ligation events within a certain range of each other. This avoids
# simulating the far-out interactions which are mostly irrelevant.
{ 
	max.range <- as.integer(max.range)
	len <- as.integer(len)
	if (max.range > len) { stop("max.range must be less than the total length of the chromosome") }

	per.height <- pmin(max.range + 1L, len-1:len+1L)
	if (!diag) { per.height <- per.height[-len] - 1L }
	total.possibles <- sum(per.height)
	if (total.possibles <= nspec) { 
		warning("reporting everything as requested") 
		chosen <- 1:total.possibles
	} else {
		chosen <- sample(total.possibles, nspec) 
	}
	
	ratchet <- c(1L, cumsum(per.height)+1L)
	target <- findInterval(chosen, ratchet) 
	anchor <- chosen - ratchet[target] + target
	if (!diag) { anchor <- anchor + 1L }
	cbind(anchor=anchor, target=target)
}

### This is the same as original Aaron Lun function
dummyFragments <- function(len, name="chrA") {
	names(len) <- name
	GRanges(name, IRanges(1:len, 1:len), seqlengths=len)
}


### This is the same as original Aaron Lun function
aggregateAll <- function(...) 
# This takes a bunch of (coordinate matrix, count vector) lists, and
# merges them to compute the aggregate mean.
{
	everything <- list(...)
	all.anchors <- all.targets <- all.counts <- list()
	counter <- 1L
	for (x in everything) { 
		all.anchors[[counter]] <- x[[1]][,1]
		all.targets[[counter]] <- x[[1]][,2]
		all.counts[[counter]] <- x[[2]]
		counter <- counter + 1L
	}

	# Reordering.
	all.anchors <- unlist(all.anchors)
	all.targets <- unlist(all.targets)
	all.counts <- unlist(all.counts)
	o <- order(all.anchors, all.targets)
	all.anchors <- all.anchors[o]
	all.targets <- all.targets[o]
	all.counts <- all.counts[o]
	
	# Finding the sum.
	keep <- which(diff(all.anchors)!=0L | diff(all.targets)!=0L)
	first <- c(1L, keep + 1L)
	last <- c(keep, length(o))
	cum.count <- cumsum(all.counts)
	collected <- cum.count[last] - cum.count[first] + all.counts[first]

	# Returning the aggregated output.
	list(cbind(anchor=all.anchors[first], target=all.targets[first]), collected)
}



#####################################################
#####
##### FROM NOW ON THE CODE IS MORE DISTAND FROM THE ORIGINAL CODE BY Aaron Lun
##### still the structure of the procedure (in term of steps and sequence of calls 
##### to the accessory functions defined above) is conserved as the underlying rationale
##### of the procedure is used here as well
#####



#####################################################
#####################################################
#####
##### ADDING VANILLA COVERAGE FUNCTION FOR NORMALIZATION
#####


# x is the squared symmetric contact matrix
VanillaCoverage<-function(OBS) {

# just make sure it is a matrix
if (!is.matrix(OBS)) {
stop("OBS must be a matrix")
}

# grand total
T<-sum(OBS)

# rows and column sums
Cs<-colSums(OBS)
Rs<-rowSums(OBS)

# expected counts given the total by column/row and the grandtotal of read counts
EXP<-((Cs) %*% t(Rs))/T

# returning the normalized matrix
return(OBS/EXP)

}



#####################################################
#####################################################


#### SAMPLE PARAMETERS FOR DEBUGGING
# simChromName<-"ChrS"
# BINSIZE<-40000
# max.tad <- 20  
# min.tad<-3
# TARGET.chrsize<-180.92*1e6 # same as chr5 in hg19
# base.tad.signal<-80
# interaction.decay.rate<-0.8
# percent_of_sampling_points_nonSpecificsignal<-0.04
# n.interactions=NULL
# base.diagonal.signal=NULL
# max.range=NULL
# base.interaction.signal=NULL
# dispersion <- 0.01 
# PNGplot=FALSE
# saveVanilla=FALSE
# saveDiffHiC=FALSE
# saveHOMER=FALSE
# addIntConstant=FALSE
# base.nonspecific.signal=NULL
# nestedTADs=TRUE
# nestedLevels=3
# nestedProb=0.25
# nestedSizeCap=75


runHiCSimulation<-function(BINSIZE, simChromName="ChrS", max.tad=50, min.tad=3,
 TARGET.chrsize=180.92*1e6, base.tad.signal=80, interaction.decay.rate=0.8,
 percent_of_sampling_points_nonSpecificsignal=0.04, n.interactions=NULL,
 base.diagonal.signal=NULL, max.range=NULL, base.interaction.signal=NULL, base.nonspecific.signal=NULL,
 dispersion=0.01, PNGplot=FALSE, saveVanilla=TRUE, saveDiffHiC=FALSE, saveHOMER=FALSE, addIntConstant=FALSE,
 nestedTADs=FALSE, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL) {



require(diffHic)
require(edgeR)


TIME<-(Sys.time()) # can be used to set seed
# in this version it is just used to define the output filename


## make this linked to the target chromosome size
n.tads <- round((TARGET.chrsize/BINSIZE)/mean(min.tad:max.tad))


#  signal in diagonal
if (is.null(base.diagonal.signal)) {
# proportional to original Aron Lun choices (28*100/80)
base.diagonal.signal<-(base.tad.signal*100/80)
}

estimated_binsperchromsome<-round((TARGET.chrsize/BINSIZE))
estimated_matrix_size_upperTriangle<-(((estimated_binsperchromsome)^2-(estimated_binsperchromsome))/2)

## %%% changed to make this proportional to changing parameters
SampledPointNonSpecificSignal<-round(percent_of_sampling_points_nonSpecificsignal * estimated_matrix_size_upperTriangle)

if (is.null(base.nonspecific.signal)) {
## %%% changed to make this proportional to changing parameters
base.nonspecific.signal<-round((5/80)*base.tad.signal) 
}

# ## original was 100 --> over matrix size 100/4840716 =  2.06581e-05 (i.e. 2 every 100k)
if (is.null(n.interactions)) {
## %%% changed to make this proportional to changing parameters
    n.interactions <-  round(2e-5*estimated_matrix_size_upperTriangle)  
}



if (is.null(max.range)) {
# max distance from diagonal for the interactions (no more than 10 bins beyond TAD max size as there are few interactions very far from diagonal in real data)
## %%% changed to make this proportional to changing parameters 
### %%% not going beyond diagonal #40... we may consider also changing this parameter
max.range <- max.tad + round(max.tad/3) 
}

if (is.null(base.interaction.signal)) {
# signal at interaction
## %%% changed to make this proportional to changing parameters
base.interaction.signal<-base.tad.signal*2
}


print("Simulate data...")

## generate random TAD partitions
## please note that in this simulation scheme the total chromosome size is not defined (thus might vary)
## instead we are setting the total number of TADs	

if (nestedTADs) {
    TadsSimulated <- spawnTADs_withNesting_dualOutput(ntads=n.tads, min.tad=min.tad, max.tad=max.tad, nestedLevels=nestedLevels, nestedProb=nestedProb, nestedSizeCap=nestedSizeCap) 
    tad.co<-TadsSimulated$tads_list
    tad.co_intervals<-TadsSimulated$tads_intervals

} else {
    TadsSimulated <- spawnTADs_dualOutput(ntads=n.tads, min.tad=min.tad, max.tad=max.tad) 
    tad.co<-TadsSimulated$tads_list
    tad.co_intervals<-TadsSimulated$tads_intervals
}


# max without capping
# > max(t(diff(t(tad.co_intervals))))
# [1] 91



### %%%  this tad.co object contains a list of TADs
### each element of the list (each TAD) contain all the pairwise coordinate of pixels belonging to the TAD
### 

# > length(tad.co)
# [1] 160
# > class(tad.co)
# [1] "list"
# > str(tad.co[[1]])
#  int [1:351, 1:2] 27 27 27 27 27 27 27 27 27 27 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : NULL
#   ..$ : chr [1:2] "anchor" "target"
# > str(tad.co[[2]])
#  int [1:105, 1:2] 42 42 42 42 42 42 42 42 42 42 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : NULL
#   ..$ : chr [1:2] "anchor" "target"
# > 
# > 
# > head(tad.co[[2]])
#      anchor target
# [1,]     42     41
# [2,]     42     40
# [3,]     42     39
# [4,]     42     38
# [5,]     42     37
# [6,]     42     36
# > head(tad.co[[1]])
#      anchor target
# [1,]     27     26
# [2,]     27     25
# [3,]     27     24
# [4,]     27     23
# [5,]     27     22
# [6,]     27     21
# > tail(tad.co[[1]])
#        anchor target
# [346,]      4      3
# [347,]      4      2
# [348,]      4      1
# [349,]      3      2
# [350,]      3      1
# [351,]      2      1
# > tail(tad.co[[2]])
#        anchor target
# [100,]     31     30
# [101,]     31     29
# [102,]     31     28
# [103,]     30     29
# [104,]     30     28
# [105,]     29     28
# > 



tad.co <- do.call(rbind, tad.co)
## %%% this matrix will contain all the pairs of coordinates corresponding to points within TADs 
# 	> head(tad.co)
#      anchor target
# [1,]     27     26
# [2,]     27     25
# [3,]     27     24
# [4,]     27     23
# [5,]     27     22
# [6,]     27     21

	
	
## this is  generating the background signal with interaction decay rate -0.8
## and constant scaling factor kt = 80, as described in the supplementary page 5
## of the original article
## here we use our own parameters

## this is generating the means of inside TAD signals as for the second equation in the first term of simulation (diffHiC suppl page 5)
tad.mu <- getMeans(tad.co, base=base.tad.signal, decay=interaction.decay.rate, rate=NA)  ## %%% also this function is not present in diffHiC
	
# Adding the diagonals proper.
nlen <- max(tad.co)
self.co <- spawnDiagon(nlen) ## 

# 
# ## NOVEL fix: adding a boost of signal to fix increased signal due to nested TADs
# if (nestedTADs) {
# 	self.mu <- rep(base.diagonal.signal*nestedLevels, nlen)   ### diagonal signal
# } else {
	self.mu <- rep(base.diagonal.signal, nlen)   ### diagonal signal
# }

	# Adding some random noise. 
	nnspec <- SampledPointNonSpecificSignal
	random.co <- spawnLigations2(nnspec, nlen, decay=interaction.decay.rate) ## %%% not defined  
	
	## The function above is just selecting 200000 data points across the entire matrix and 
	## the step below is defining the value 5 that will be added randomly to the mean of the 200k data points selected above
	## the same bloc of x,y coordinates in principle could be sampled twice
	## in that case the mu value 5 will be added multiple times
	## the only restriction is that x!=y
	##moreover the probability of being sampled is proportional to (x-y+p)^c
# 	 > table(duplicated(paste(random.co[,2], random.co[,1], sep="_")))
# 
#  FALSE   TRUE 
# 176734  23266 



	random.mu <- rep(base.nonspecific.signal, nnspec)


### %%% just checking what is the maximum number of added noise background

randomNoiseTotals<-tapply(random.mu, INDEX=paste("anch", random.co[,1], "tar", random.co[,2], "bin", sep="_"), FUN=sum)

## some of them are as big as 30% of the signal of true interactions (160 in the default settings)
# > head(sort(randomNoiseTotals))
#   anch_10_tar_1_bin   anch_10_tar_3_bin   anch_10_tar_8_bin anch_100_tar_15_bin 
#                   5                   5                   5                   5 
# anch_100_tar_42_bin anch_100_tar_52_bin 
#                   5                   5 
# > tail(sort(randomNoiseTotals))
#   anch_801_tar_800_bin   anch_813_tar_812_bin anch_1035_tar_1034_bin 
#                     40                     40                     45 
# anch_2836_tar_2835_bin   anch_301_tar_300_bin   anch_896_tar_894_bin 
#                     45                     50                     50 

	# Adding some differential interactions 
	## %%% in the original simulation procedure this was used to define points of differential interactions between two conditions (with replicates
	## but we can use it also to just define points of interactions and to a peak call. 
	## later we may need to adjust the parameters because the statistics to detect differential interactions might be more sensitive than the ones to call peaks in individual samples
	## as such we may need to increase the signal for each interaction point
	
	
 	inter.co <- spawnLigations(n.interactions, max.range, nlen, diag=FALSE)  ## %%% not defined
	inter.mu <- getMeans(inter.co, base=base.interaction.signal, decay=interaction.decay.rate, rate=NA) ## %%% not defined

## add an extra constant to the interaction signal if required
## same intensity as the nonspecific signal
## to get better separation from nonspecific signal
if (addIntConstant) {
    inter.mu <- (inter.mu + base.nonspecific.signal)
}


## not using this as it was part of two groups simulation for differential interactions calls
# 	chosen <- logical(n.interactions)
# 	chosen[sample(n.interactions, n.interactions/2)] <- TRUE  ## %%% picking just half of the sampled interactions beacause below it is adding them only to one of the two  groups alteranted
# 	## %% need to adjust to have only one set of TRUE interactions


	# Running through for each width of interaction.
	blocks <- dummyFragments(nlen, name=simChromName) ## %%% not defined
	xparam <- pairParam(blocks) 


## %% need to adjust to have only one set of TRUE interactions
# 		for (grouping in 1:2) {
cur.mus <- inter.mu
# 		if (grouping==1L) {
cur.inters <- inter.co # keeping everything as we simulate just one group
# 			cur.mus[chosen] <- 0
# 		} else {
# 			cur.inters <- inter.co[chosen,]
# 			cur.mus[!chosen] <- 0
# 		}

## summing up all the contributions to mean values from the various equations
collated <- aggregateAll( list(tad.co, tad.mu), list(self.co, self.mu), list(random.co, random.mu), list(cur.inters, cur.mus) )


# > table(duplicated(paste(collated[[1]][,1], collated[[1]][,2], sep="_")))
# 
#  FALSE 
# 193221 


## generate read counts
## it looks like the size parameter here is defined as 1/theta (theta=dispersion)
## there seems to be some confusion in the field as
## QUOTE
## in NB distribution the variance can then be written as (m + (m^2)/r)
## some authors set a= 1/r thus the variance is (m +am^2)
## In this context, and depending on the author, either the parameter r or its reciprocal a is referred to as the “dispersion parameter”, “shape parameter”)
## end QUOTE
## in the rnbinom function the size parameter should be the "r" shape parameter
## Aron default dispersion 0.01 instead is most likely referring to 1/r as he is taking the reciprocal to get the shape parameter.
## Indeed this is the case as the shape parameter should be an integer (although not necessarily according to rnbinom manual) in this case 1/0.01 = 100



counts <- as.integer(rnbinom(length(collated[[2]]), mu=collated[[2]], size=1/dispersion))  ## %%% here it is generating the 
## read counts based on negative binomial distributions are generated by sampling from neg binom distribution with means defined above
			

# > length(counts)
# [1] 193221
# > nrow(collated[[1]])
# [1] 193221



####################################
####################################
###
### FORMATTING THE OUTPUT AS FULL (SQUARED) INTERACTION MATRIX
###

COUNT.MATRIX<-matrix(0, nrow=nlen, ncol=nlen)


COUNT.MATRIX[(collated[[1]])]<-counts

## we want a symmetric matrix
COUNT.MATRIX[((collated[[1]])[,c(2,1)])]<-counts


## we want a symmetric matrix

# > head(collated[[1]])
#      anchor target
# [1,]      1      1
# [2,]      2      1
# [3,]      2      2
# [4,]      3      1
# [5,]      3      2
# [6,]      3      3
# > head(((collated[[1]])[,c(2,1)]))
#      target anchor
# [1,]      1      1
# [2,]      1      2
# [3,]      2      2
# [4,]      1      3
# [5,]      2      3
# [6,]      3      3

# > COUNT.MATRIX[(collated[[1]])]<-counts
# > COUNT.MATRIX[1:10,1:10]
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#  [1,]  103    0    0    0    0    0    0    0    0     0
#  [2,]   57  107    0    0    0    0    0    0    0     0
#  [3,]   37   65  101    0    0    0    0    0    0     0
#  [4,]   33   39   50   97    0    0    0    0    0     0
#  [5,]   30   21   41   63  104    0    0    0    0     0
#  [6,]   31   25   24   36   56  105    0    0    0     0
#  [7,]   23   30   22   40   51   50  114    0    0     0
#  [8,]   25   23   18   36   30   31   52  105    0     0
#  [9,]   27   18   26   27   28   20   31   62   85     0
# [10,]   18   15   20   27   12   16   25   48   73    92

# > ## we want a symmetric matrix
# > COUNT.MATRIX[((collated[[1]])[,c(2,1)])]<-counts
# > COUNT.MATRIX[1:10,1:10]
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#  [1,]  103   57   37   33   30   31   23   25   27    18
#  [2,]   57  107   65   39   21   25   30   23   18    15
#  [3,]   37   65  101   50   41   24   22   18   26    20
#  [4,]   33   39   50   97   63   36   40   36   27    27
#  [5,]   30   21   41   63  104   56   51   30   28    12
#  [6,]   31   25   24   36   56  105   50   31   20    16
#  [7,]   23   30   22   40   51   50  114   52   31    25
#  [8,]   25   23   18   36   30   31   52  105   62    48
#  [9,]   27   18   26   27   28   20   31   62   85    73
# [10,]   18   15   20   27   12   16   25   48   73    92


### pairs of x and y indexes for a matrix in R are handled with cbind (or a matrix with two columns) see example below 
# > 
# > COUNT.MATRIX
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#  [1,]    1   11   21   31   41   51   61   71   81    91
#  [2,]    2   12   22   32   42   52   62   72   82    92
#  [3,]    3   13   23   33   43   53   63   73   83    93
#  [4,]    4   14   24   34   44   54   64   74   84    94
#  [5,]    5   15   25   35   45   55   65   75   85    95
#  [6,]    6   16   26   36   46   56   66   76   86    96
#  [7,]    7   17   27   37   47   57   67   77   87    97
#  [8,]    8   18   28   38   48   58   68   78   88    98
#  [9,]    9   19   29   39   49   59   69   79   89    99
# [10,]   10   20   30   40   50   60   70   80   90   100
# > 
# > 
# > COUNT.MATRIX[c(1,3,5), c(5,3,1)]
#      [,1] [,2] [,3]
# [1,]   41   21    1
# [2,]   43   23    3
# [3,]   45   25    5
# > COUNT.MATRIX[cbind(c(1,3,5), c(5,3,1))]
# [1] 41 23  5


print("Save counts...")

Bin_labels<-paste(simChromName, 1:nlen, sep="_")

colnames(COUNT.MATRIX)<-Bin_labels


OUT.COUNT.MATRIX<-data.frame("Bin"=Bin_labels, COUNT.MATRIX, stringsAsFactors=FALSE)
# Bin	ChrX_1	ChrX_2	ChrX_3	ChrX_4	ChrX_5	ChrX_6	ChrX_7	ChrX_8	ChrX_9	ChrX_10	


OUTFILES_suffix<-paste(simChromName, "Bin", BINSIZE, "params", paste(max.tad, min.tad, n.tads, base.tad.signal,
 interaction.decay.rate, base.diagonal.signal, SampledPointNonSpecificSignal,
  base.nonspecific.signal, n.interactions, max.range, base.interaction.signal,
  dispersion, ifelse(addIntConstant, yes="addIntConstant", no=""), 
  ifelse(nestedTADs, yes=paste("nestedTADs", nestedLevels, nestedProb, sep="_"), no=""),
  ifelse((!is.null(nestedSizeCap)), yes=paste("cap", nestedSizeCap, sep=""), no=""),
  sep="_"), "TIME", format(TIME, "%Y-%m-%d_%Hh%Mm%Ss"), sep=".")


## write counts
outfilename_counts<-paste("simHiC_countMatrix", OUTFILES_suffix, "txt", sep=".")
write.table(OUT.COUNT.MATRIX, file=outfilename_counts, sep="\t", row.names=FALSE, quote=FALSE)


if (saveVanilla) {

    print("Save vanilla...")

    ## write Vanilla Coverage normalized matrix

    VANILLA.MATRIX<-VanillaCoverage(COUNT.MATRIX)
    colnames(VANILLA.MATRIX)<-Bin_labels

    OUT.VANILLA.MATRIX<-data.frame("Bin"=Bin_labels, VANILLA.MATRIX, stringsAsFactors=FALSE)
    outfilename_vanilla<-paste("simHiC_vanillaCoverageMatrix", OUTFILES_suffix, "txt", sep=".")
    write.table(OUT.VANILLA.MATRIX, file=outfilename_vanilla, sep="\t", row.names=FALSE, quote=FALSE)

}



if (saveDiffHiC | saveHOMER) {
   
print("Prepare Data for HOMER or diffHiC")
   
    # data to be saved in output (more than zero counts)
    diffHic_record <- counts > 0L
    diffHic_anchor.id <- collated[[1]][diffHic_record,1]
    diffHic_target.id <- collated[[1]][diffHic_record,2]
    diffHic_counts <- counts[diffHic_record]
    diffHic_anchor.id_Final <- rep(diffHic_anchor.id, diffHic_counts)
    diffHic_target.id_Final <- rep(diffHic_target.id, diffHic_counts)
    
    # generate precise reads location by sampling positions from 25% to 75% of bin positions
    diffHic_anchor.loc <- diffHic_anchor.id_Final*BINSIZE - as.integer(runif(sum(diffHic_counts), 0.25*BINSIZE, 0.75*BINSIZE))
    diffHic_target.loc <- diffHic_target.id_Final*BINSIZE - as.integer(runif(sum(diffHic_counts), 0.25*BINSIZE, 0.75*BINSIZE))
    diffHic_anchor.str <- sample(c("+", "-"), sum(diffHic_counts), replace=TRUE)
    diffHic_target.str <- sample(c("+", "-"), sum(diffHic_counts), replace=TRUE)
}


### save data in diffHiC format 
if (saveDiffHiC) {

    print("Save diffHiC Format...")
 
    outfilename_diffHic<-paste("simHiC_diffHiCformat", OUTFILES_suffix, "h5", sep=".")

    # Saving for use in diffHic.
    savePairs(x=data.frame(anchor.id=diffHic_anchor.id_Final, target.id=diffHic_target.id_Final, 
    anchor.pos=diffHic_anchor.loc*ifelse(diffHic_anchor.str=="-", -1L, 1L), anchor.len=10,
    target.pos=diffHic_target.loc*ifelse(diffHic_target.str=="-", -1L, 1L), target.len=10), 
    file=outfilename_diffHic, param=xparam)

    outfilename_diffHic_xparam<-paste("simHiC_diffHiC_xparam", OUTFILES_suffix, "RData", sep=".")

    # Saving for use in diffHic.
    save(xparam, file=outfilename_diffHic_xparam)

}



### save data in HOMER format
if (saveHOMER) {

    print("Save HOMER Format...")

  outfilename_HOMER<-paste("simHiC_HOMERformat", OUTFILES_suffix, "txt", sep=".")

    # Saving into an alternative form for HOMER (saving as characters to avoid abbreviation, e.g., to 5e5).
outfilename_HOMER_bzipped<-paste(outfilename_HOMER, "bz2", sep=".")


t0<-proc.time()[3]
    HOMERdata<-data.frame("", simChromName, as.character(diffHic_anchor.loc), diffHic_anchor.str,
      simChromName, as.character(diffHic_target.loc), diffHic_target.str)
t1<-proc.time()[3]
print(paste("elapsed", (t1-t0)))
     
    outfilename_HOMER_bzipped_handler<-bzfile(outfilename_HOMER_bzipped, open = "w")
    write.table(HOMERdata, file=outfilename_HOMER_bzipped_handler, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(outfilename_HOMER_bzipped_handler)
rm(HOMERdata)
gc()

# HOMER FORMATTED DATA MUST BE post-processed and formatted with this command
# system(sprintf("makeTagDirectory %s -format HiCsummary %s", curfname, curfile.homer))
# makeTagDirectory [yourChoicOfDirectoryName] -format HiCsummary [nomefile simHiC_HOMERformat.txt]



}





### MUST WRITE OUT ALSO THE TAD BORDERS AND THE INTERACTION POSITIONS
print("Save TADs...")

## write TADS
outfilename_tads<-paste("simHiC_TADintervals", OUTFILES_suffix, "txt", sep=".")
write.table(tad.co_intervals, file=outfilename_tads, sep="\t", row.names=FALSE, quote=FALSE)


print("Save interactions...")

## write interactions
outfilename_interactions<-paste("simHiC_interactions", OUTFILES_suffix, "txt", sep=".")
write.table(inter.co, file=outfilename_interactions, sep="\t", row.names=FALSE, quote=FALSE)
# already sorted
# > table(inter.co[,1] >= inter.co[,2])
# 
# TRUE 
#  100 


## %%% better to chek how much 5x noise is added to individual points?
## write logs
MaxRandomNoisePerBin<-max(randomNoiseTotals)
outfilename_logs<-paste("simHiC_logs", OUTFILES_suffix, "txt", sep=".")
write(file=outfilename_logs, paste("MaxRandomNoisePerBin", MaxRandomNoisePerBin))



##########################################################################################
##########################################################################################
######## 
######## PLOT SIMULATED DATA
######## 
######## 

capping_percentile<-1
Binsize<-BINSIZE
PEAKscolor<-"red"

# log transform the data
data.forPLOT<-log10(COUNT.MATRIX+1)

### Fetch the top 97th quantile
TopNthPercentile <- quantile(data.forPLOT,capping_percentile)


print("Get colors")
BreaksRange <- c(seq(from=min(data.forPLOT),to=TopNthPercentile,length.out=100), max(data.forPLOT))

	# color palette
	Colors <- c("#ffffff", "#353535")
	InterpolateColours <- colorRampPalette(Colors)(length(BreaksRange)-1)


   # plot the heatmap
    x_coordinates_bins<-1:nrow(data.forPLOT)
    y_coordinates_bins<-1:ncol(data.forPLOT)

    # bin number*Binsize is the end of the bin interval coordinate thus need to add the first one and subtract one
    x_coordinates<-c((x_coordinates_bins[1]-1), x_coordinates_bins)*Binsize
    y_coordinates<-c((y_coordinates_bins[1]-1), y_coordinates_bins)*Binsize


  # step for axes labels
    AxisLabel.step.Genome<-10000000

# define the coordinates for peaks plotting
filteredPeaks<-data.frame(plotting1=((inter.co[,1]*Binsize)-(Binsize/2)), plotting2=((inter.co[,2]*Binsize)-(Binsize/2)))

if (PNGplot) {

print("Save PNG...")

    ##############################
    # draw the plot FOR THE WHOLE CHROMOSOME
#   pdf(file=paste("simHiC_PLOT_interactions", OUTFILES_suffix, "pdf", sep="."), height=20, width=22)
     png(file=paste("simHiC_PLOT_interactions", OUTFILES_suffix, "png", sep="."), height=1500, width=1650)
     layout(matrix(1:2, nrow=1), widths=c(10,3))

    #making fonts smaller
    par(cex=1.5, cex.axis=1.5, mar=c(5,4,2,0.5))

    ## plot the heatmap
    image(z=data.forPLOT, x=x_coordinates, y=y_coordinates, col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))
#     image(z=data.forPLOT[1:20,1:20], x=x_coordinates[1:21], y=y_coordinates[1:21], col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))
# LIMIT_PLOT<-401
#     image(z=data.forPLOT[1:LIMIT_PLOT,1:LIMIT_PLOT], x=x_coordinates[1:(LIMIT_PLOT+1)], y=y_coordinates[1:(LIMIT_PLOT+1)], col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))

   # place x-axis labels
    x_labels_positions<-x_coordinates[which((x_coordinates %% AxisLabel.step.Genome)==0)]
    x_labels<-paste(x_labels_positions/1000000, "Mb", sep="")
    axis(side=1, at=x_labels_positions, labels=x_labels, las=2)

   # place y-axis labels
    y_labels_positions<-y_coordinates[which((y_coordinates %% AxisLabel.step.Genome)==0)]
    y_labels<-paste(y_labels_positions/1000000, "Mb", sep="")
    axis(side=2, at=y_labels_positions, labels=x_labels, las=2)

    points(filteredPeaks$plotting2,filteredPeaks$plotting1,pch=16,cex=0.5,col=PEAKscolor)


        par(mar=c(15,2,15,10))
          image(y=1:length(BreaksRange), z=matrix(1:(length(BreaksRange)-1), nrow=1), col=InterpolateColours, xlab="", ylab="", yaxt="n", xaxt="n")
          box()

    	  legend_ticks_selected<-seq(1,(length(BreaksRange)-1), length.out=4)
    	  legend_ticks_labels<-round(BreaksRange[legend_ticks_selected], digits=3)
    	  legend_ticks_labels[length(legend_ticks_labels)]<-paste(">=", legend_ticks_labels[length(legend_ticks_labels)], sep="")
          axis(side=4, at=legend_ticks_selected, labels=legend_ticks_labels, las=2)


    dev.off()

}

    ##############################
    ##############################
    ##############################
    # draw the plot FOR a ZOOM IN REGIONS
print("Save PDF...")

   pdf(file=paste("simHiC_PLOT_zoomIN", OUTFILES_suffix, "pdf", sep="."), height=20, width=22)
#     png(file=paste("simHiC_PLOT_interactions", OUTFILES_suffix, "png", sep="."), height=1500, width=1650)
     layout(matrix(1:2, nrow=1), widths=c(10,3))

    #making fonts smaller
    par(cex=1.5, cex.axis=1.5, mar=c(5,4,2,0.5))

    ## plot the heatmap
#     image(z=data.forPLOT, x=x_coordinates, y=y_coordinates, col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))
#     image(z=data.forPLOT[1:20,1:20], x=x_coordinates[1:21], y=y_coordinates[1:21], col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))
 LIMIT_PLOT<-400
     image(z=data.forPLOT[1:LIMIT_PLOT,1:LIMIT_PLOT], x=x_coordinates[1:(LIMIT_PLOT+1)], y=y_coordinates[1:(LIMIT_PLOT+1)], col=InterpolateColours, breaks=BreaksRange, yaxt="n", xaxt="n", xlab="", ylab="", main=paste(simChromName, n.interactions, "interactions"))

   # place x-axis labels
    x_labels_positions<-x_coordinates[which((x_coordinates %% AxisLabel.step.Genome)==0)]
    x_labels<-paste(x_labels_positions/1000000, "Mb", sep="")
    axis(side=1, at=x_labels_positions, labels=x_labels, las=2)

   # place y-axis labels
    y_labels_positions<-y_coordinates[which((y_coordinates %% AxisLabel.step.Genome)==0)]
    y_labels<-paste(y_labels_positions/1000000, "Mb", sep="")
    axis(side=2, at=y_labels_positions, labels=x_labels, las=2)

    points(filteredPeaks$plotting2,filteredPeaks$plotting1,pch=16,cex=0.5,col=PEAKscolor)


        par(mar=c(15,2,15,10))
          image(y=1:length(BreaksRange), z=matrix(1:(length(BreaksRange)-1), nrow=1), col=InterpolateColours, xlab="", ylab="", yaxt="n", xaxt="n")
          box()

    	  legend_ticks_selected<-seq(1,(length(BreaksRange)-1), length.out=4)
    	  legend_ticks_labels<-round(BreaksRange[legend_ticks_selected], digits=3)
    	  legend_ticks_labels[length(legend_ticks_labels)]<-paste(">=", legend_ticks_labels[length(legend_ticks_labels)], sep="")
          axis(side=4, at=legend_ticks_selected, labels=legend_ticks_labels, las=2)


    dev.off()



return(OUTFILES_suffix)

}



###########################################################################
###########################################################################
###########################################################################
###########################################################################
#####
#####
##### RUN SIMULATION LOOP
#####



KtBaselineTADs<-28
nestedLevels<-3
 
for (iteration in (1:5)[1]) {
  for (noisePercentage in c(0.04, 0.08, 0.12, 0.16, 0.20)) {

    print(paste("iteration:", iteration))
    print(paste("noisePercentage:", noisePercentage))


    check<-runHiCSimulation(BINSIZE=40000, simChromName="chr5", max.tad=50, min.tad=3,
     TARGET.chrsize=180.92*1e6, base.tad.signal=round(KtBaselineTADs/nestedLevels), interaction.decay.rate=0.69,
     percent_of_sampling_points_nonSpecificsignal=noisePercentage, n.interactions=205, 
     base.diagonal.signal=(KtBaselineTADs*100/80), max.range=(50+round(50*1/3)), base.interaction.signal=(KtBaselineTADs*2), 
     base.nonspecific.signal=round((5/80)*KtBaselineTADs), dispersion=0.01, PNGplot=FALSE, saveVanilla=FALSE,
     saveDiffHiC=FALSE, saveHOMER=FALSE, addIntConstant=FALSE,
      nestedTADs=TRUE,nestedLevels=nestedLevels, nestedProb=0.25, nestedSizeCap=75)

    print(paste("completed", check))
  }
}






