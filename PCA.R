#all code performed on R

#may be beneficial to add more comments throughout the code

#set working directory
setwd("Path/To/File")

#for warning message generation
options(warn=1)

#install all relevant packages 
install.package("adegenet")
install.package("adegraphics")
install.package("vcfR")
install.package("pegas")
install.package("stAMPP")
install.package("ade4")
install.package("MASS")

#load all relevant packages
library(adegenet)
library(adegraphics)
library(vcfR)
library(pegas)
library(stAMPP)
library(ade4)
library(MASS)

#################################################################################
# start of modified function used to convert tetraploids into genlight objects  #
#################################################################################

vcfR2genlight.tetra <- function(x, n.cores = 1)

{
	bi <- is.biallelic(x)

	if (sum(!bi > 0)) {
		msg <- paste("Found", sum!(bi), "Loci with more than two alleles.")
		msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
		msg <- c(msg, "/n", paste(sum(!bi), "Loci will be omittd from the genlight object."))
		warning(msg)
		x <- x[bi, ]
	}

	x <- addID(x)

	CHROM <- x@fix[, "CHROM"]
	POS <- x@fix[, "POS"]
	ID <- x@fix[, "ID"]

	x <- extract.gt(x)

	x[x == "0|0"] <- 0
  	x[x == "0|1"] <- 1
  	x[x == "1|0"] <- 1
  	x[x == "1|1"] <- 2
  	x[x == "0/0"] <- 0
  	x[x == "0/1"] <- 1
  	x[x == "1/0"] <- 1
  	x[x == "1/1"] <- 2
  	x[x == "1/1/1/1"] <- 4
  	x[x == "0/1/1/1"] <- 3
  	x[x == "0/0/1/1"] <- 2
  	x[x == "0/0/0/1"] <- 1
  	x[x == "0/0/0/0"] <- 0
  	x[x == "0/0/0/0/0/0"] <- 0
  	x[x == "0/0/0/0/0/1"] <- 1
  	x[x == "0/0/0/0/1/1"] <- 2
  	x[x == "0/0/0/1/1/1"] <- 3
  	x[x == "0/0/1/1/1/1"] <- 4
  	x[x == "0/1/1/1/1/1"] <- 5
  	x[x == "1/1/1/1/1/1"] <- 6

  	if (requireNamespace("adgenet")) {
  		x <- new("genlight", t(x), n.cores = n.cores)
  	}

  	else {
  		warning("adgenet not installed")
  	}

  	adgenet::chromosome(x) <- CHROM
  	adgenet::position(x) <- POS
  	adgenet::locNames(x) <- ID
  	return(x)
}

##############################
#  end of modified function  #
##############################

########################################################################
#  start of function to calculate PCA on genlight objects much faster  #
########################################################################

glPcaFast < function(x, center = TRUE, scale = FALSE, nf = NULL, loadings = TRUE, alleleAsUnit = FALSE, returnDotProd = FALSE){
	if(!inherits(x, "genlight")) stop("x is not a genlight object.")
	if (center) {
		vecMeans <- glMean(x, alleleAsUnit = alleleAsUnit)
		if(any(is.na(vecMeans))) stop("NAs detected in the vector of means.")
	}
	if(scale){
		vecVar <- glVar(x, alleleAsUnit = alleleAsUnit)
	}
	if(any(is.na(vecVar))) stop("NAs detected in the vector of variances.")
}

#dividing full data by the ploidy data will keep the non applicable (NA) values
mx <- t(sapply(x$gen, as.integer)) / ploidy(x)

NAidx <- which(is.na(mx), arr.ind = T)
if(center) {
	mx[NAidx] <- vecMeans[NAidx[,2]]
} else {
	mx[NAidx] <- 0
}

#center and scale the data
mx <- scale(mx, center = if(center) vecMeans else F, scale = if(scale) vecVar else F)

#conduct eigen analysis
allProd <- tcrossprod(mx) / nInd(x)

eigRes <- eigen(allProd, symmetric = TRUE, only.values = FALSE)

rank <- sum(eigRes$values > 1e-12)

eigRes$values <- eigRes$vectors[, 1:rank, drop = FALSE]

if(is.null(nf)){
	barplot(eigRes$values, main = "Eigen Values", col = heat.colors(rank))
	cat("Select the number of axes: ")
	nf <- as.integer(readLines(n = 1))
}

#rescale the PCs
res <- list()
res$eig <- eigRes$values
nf <- min(nf, sum(res$eig > 1e-10))

#use for debugging
res$matprod <- allProd 
li = XQU = V/Lambda^(1/2)

eigRes$vectors <- eigRes$vectors * sqrt(nInd(x))

res$scores <- sweep(eigRes$vectors[, 1:nf, drop = FALSE], 2, sqrt(eigRes$valies[1:nf]), FUN = "*")

if(loadings){
	if(scale){
		vecSd <- sqrt(vecVar)
	} 
	res$loadings <- matrix(0, nrow = nloc(x), ncol = nf) #make an empty matrix
	}
	myPloidy <- ploidy(x)
	for(k in 1:nInd(x)){
		temp <- as.integer(X@gen[[k]]) / myPloidy[k]
		if(center) {
			temp[is.na(temp)] <- vecMeans[is.na(temp)]
			temp <- temp - vecMeans
		} else {
			temp[is.na(temp)] <- 0
			}
		if(scale) {
			temp <- temp/vecSd
		}
		res$loadings <- res$loadings + matrix(temp) %>% eigRes$vectors[k, 1:nf, drop = FALSE]
		}
		res$loadings <- res$loadings / Nind(x)
		res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN = "/")
	}

#format the output
colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)

#########################
#  end of PCA function  #
#########################

#import SNP data from VCF file and read into an object
vcf <- read.vcfR("file_name.vcf.gz")

#convert VCF to genlight with modified function
aa.genlight <- vcfR2genlight.tetra(vcf) 

#add SNP names
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")

#add population names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3) 

#check
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)

#use the modified PCA function for increased efficiency
pca.1 < glPcaFast(aa.genlight, nf = 300)

#plot the PCA scatter graph with individual labels
scatter(pca.1, posi = "bottomright")
loadigplot(pca.1)

#explained variance proportion
pca.1$eig[1]/sum(pca.1$eig) #1st axis
pca.1$eig[2]/sum(pca.1$eig) #2nd axis
pca.1$eig[3]/sum(pca.1$eig) #3rd axis

#colour the populations in the graph
col <- funky(10)
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F)

#save the figures in a PDF format
pdf ("PCA.pdf", width=14, height=7)

g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)

g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)

ADEgS(c(g1, g2), layout = c(1, 2))

dev.off()
