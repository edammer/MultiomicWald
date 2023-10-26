## The Wald test used by DESeq 2 is converting a Z statistic (log2FC/SE(Log2FC)) to p value, and multiplying by 2 for both tails of the distribution.
# e.g., this is the DESeq Z (stat) for a decreasing gene in a 2-group contrast with p=0.05
# pnorm(-1.96,lower.tail=TRUE)*2

## Goals:
## Demonstrate what is the appropriateness of different SE or log2FC calculations contributing to the above Z statistic used for Wald p value.
## Test the interoperability of the negative binomial framework for modeling data that can be used to produce the ratioed values for Wald statistic, using either integer or continuous input.

## Literature Resources:
# DESeq2 methods (2014):  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8  #* This is the article describing the math in DESeq2.
# Ref #58 from DESeq2 methods:  https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372


## Set up working folder
rootdir="E:/OneDrive/statTest_dammer/"
setwd(rootdir)

## Load Traits
#traits<-read.csv("Exp. 9-2 & 10-1 combined traits.csv",header=TRUE,row.names=1,check.names=FALSE)
repsPerGroup=3
traits<-data.frame(group.simple=c( rep("Group1",repsPerGroup),
                                   rep("Group2",repsPerGroup)) )
traits$group.simple<-factor(traits$group.simple)

####################################
## Read integer counts matrix
#counts <- read.csv("../Exp. 9-2 & 10-1 mRNA-seq analysis/Exp. 9-2 & 10-1 mRNA-seq analysis/Exp. 9-2 & 10-1 combined counts_filtered_05222023.csv", header = TRUE, row.names = 1,check.names=FALSE)
#
# # Drop rows with NA
#check <- counts
#cleanDat.RNA.integer <- check[rowSums(is.na(check)) == 0, ]       
#
#dim(cleanDat.RNA.integer)
# #[1] 17038 29

## Simulate integer counts matrix similar to above -- using PROPER package by Hao Wu @ emory.edu
#seqDepth=16000000
simOptions<-PROPER::RNAseq.SimOptions.2grp(ngenes=46000,seqDepth=1,lBaselineExpr=function(x) rgamma(x/2,shape=4.85,rate=0.525), lOD=function(x) rgamma(x,shape=0.85,rate=6), p.DE=0.05, lfc=function(x) rnorm(x, mean=0, sd=1.75), sim.seed=1)
simCounts<- PROPER::simRNAseq(simOptions, n1=repsPerGroup,n2=repsPerGroup)
counts<-simCounts$counts
dim(counts)
#[1] 46000     6
top40pct.idx<-which(rowSums(counts)>quantile(rowSums(counts),probs=0.60))
bottom50pct.idx<-which(rowSums(counts)<quantile(rowSums(counts),probs=0.50))
counts<-counts[-top40pct.idx,]
counts<-counts[-sample(bottom50pct.idx,length(bottom50pct.idx)*0.42857*1.75),]
dim(counts)
#[1] 17299     6
counts<-round(counts/2,0)
colSums(counts)  # read depth simulated, per sample
#[1] 21819205 20838332 21462486 21381287 21337288 21193346

#y=rowMeans(counts)
#y[y>2000]<-NA
#hist(y,breaks=100)
#hist(z,breaks=100,add=TRUE,col="#FF888855")  # baseMean distribution of actual counts (17038 row means), matches histogram after modeling and reshaping

cleanDat.RNA.integer<-apply(counts,2,as.integer)
rownames(traits)<-colnames(cleanDat.RNA.integer)<-paste0("sample",1:ncol(cleanDat.RNA.integer))
rownames(cleanDat.RNA.integer)<-paste0("gene",1:nrow(cleanDat.RNA.integer))

####################################

library(DESeq2)

## Run DESeq2 for Wald priors, and the native DESeq2 calculations of Wald statistics
dds2 <- DESeqDataSetFromMatrix(countData = cleanDat.RNA.integer,    #[,which(traits$group.simple=="Group1" | traits$group.simple=="Group2")],
                              colData = traits,                     #[which(traits$group.simple=="Group1" | traits$group.simple=="Group2"),],
                              design = ~group.simple) 

dds2 <- DESeq(dds2,test="Wald",fitType="mean")

as.data.frame(mcols(mcols(dds2), use.names=TRUE))
#                                                     type                                            description
#baseMean                                     intermediate              mean of normalized counts for all samples
#baseVar                                      intermediate          variance of normalized counts for all samples
#allZero                                      intermediate                         all counts for a gene are zero
#dispGeneEst                                  intermediate                      gene-wise estimates of dispersion
#dispGeneIter                                 intermediate                     number of iterations for gene-wise
#dispFit                                      intermediate                            fitted values of dispersion
#dispersion                                   intermediate                           final estimate of dispersion
#dispIter                                     intermediate                                   number of iterations
#dispOutlier                                  intermediate                          dispersion flagged as outlier
#dispMAP                                      intermediate                          maximum a posteriori estimate
#Intercept                                         results                      log2 fold change (MLE): Intercept
#group.simple_Group2_vs_Group1                    results log2 fold change (MLE): group.simple Group2 vs Group1
#SE_Intercept                                      results                              standard error: Intercept
#SE_group.simple_Group2_vs_Group1                 results         standard error: group.simple Group2 vs Group1
#WaldStatistic_Intercept                           results                              Wald statistic: Intercept
#WaldStatistic_group.simple_Group2_vs_Group1      results         Wald statistic: group.simple Group2 vs Group1
#WaldPvalue_Intercept                              results                           Wald test p-value: Intercept
#WaldPvalue_group.simple_Group2_vs_Group1         results      Wald test p-value: group.simple Group2 vs Group1
#betaConv                                          results                                   convergence of betas
#betaIter                                          results                                   iterations for betas
#deviance                                          results                          deviance for the fitted model
#maxCooks                                          results                        maximum Cook's distance for row

data.frame(mcols(dds2,use.names=TRUE))[1:4,4:10]
#DataFrame with 4 rows and 7 columns
#      dispGeneEst dispGeneIter dispFit dispersion dispIter dispOutlier  dispMAP
#gene1    1.854090            8 1.15171   1.341609       12       FALSE 1.341609
#gene2    1.672975            8 1.15171   1.371343        9       FALSE 1.371343
#gene3    1.037232           11 1.15171   1.103074        8       FALSE 1.103074
#gene4    1.435150           12 1.15171   1.270084        7       FALSE 1.270084

## Pull out dispersion statistics for use in our own glm.nb calculations.
dispMat=data.frame(mcols(dds2,use.names=TRUE))[,4:10]


## Set up parallel backend
parallelThreads=31

require(doParallel, quietly=TRUE)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)


## Attempt 1 of 2  to reproduce Wald Z stat calculations of DESeq2 using integer counts
statOut <-foreach(this.gene=as.character(rownames(cleanDat.RNA.integer)), .combine='rbind') %dopar% {
	#this.gene="Gene1"
	
	thisGene.exp<-unlist(as.vector(cleanDat.RNA.integer[this.gene,which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]))
	Grouping=traits$group.simple[which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]

	lmdat.int=data.frame(exp=thisGene.exp,Grouping=Grouping)

	## Try negative binomial glm with theta from final DESeq2 gene specific dispersion; that failing, fallback to glm.nb with automatic theta estimation; or that failing, use lm
	model <- tryCatch( glm(exp~Grouping,data=lmdat.int,family=MASS::negative.binomial(dispMat[this.gene,"dispersion"])), error= function(e) tryCatch( MASS::glm.nb(exp ~ Grouping, data = lmdat.int), error= function(e2) lm(exp~Grouping,data = lmdat.int) ))

	SE=sqrt(diag(vcov(model))[2])
#	LFC=coef(model)[2]  # this is used in statOut2
	#ALT LFC  -- mean(log2(counts)), with 0 counts replaced by 0.5
	group1.counts.0to0.5<-lmdat.int$exp[which(Grouping==names(table(Grouping))[1])]
	group2.counts.0to0.5<-lmdat.int$exp[which(Grouping==names(table(Grouping))[2])]
	group1.counts.0to0.5[group1.counts.0to0.5==0]<- 0.5
	group2.counts.0to0.5[group2.counts.0to0.5==0]<- 0.5
	LFC=mean(log2(group1.counts.0to0.5))-mean(log2(group2.counts.0to0.5))
	stat=LFC/SE
	lowerTail=if(LFC<0) { TRUE } else { FALSE }
	pvalue=pnorm(LFC/SE,lower.tail=lowerTail)*2
	if(pvalue>1) pvalue=1

	return(c(Gene=this.gene,baseMean=mean(lmdat.int$exp),LFC=LFC,lfcSE=SE,stat=stat,pvalue=pvalue))
}
statOut.rownames=statOut[,"Gene"]
statOut<-apply(statOut[,2:ncol(statOut)],2,as.numeric)
rownames(statOut)<-statOut.rownames
colnames(statOut)<-as.data.frame(do.call(rbind,strsplit(colnames(statOut),"[.]")))[,1]

statOut<-cbind(statOut, p.adjust(statOut[,"pvalue"],method="BH"))
colnames(statOut)[ncol(statOut)]<-"padj"

#write.csv(statOut,"statOut1.integer.csv")



## Repeat stats table calculation with negative binomial glm used for both log2FC and lfcSE
statOut2 <-foreach(this.gene=as.character(rownames(cleanDat.RNA.integer)), .combine='rbind') %dopar% {
	#this.gene="Gene1"
	
	thisGene.exp<-unlist(as.vector(cleanDat.RNA.integer[this.gene,which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]))
	Grouping=traits$group.simple[which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]

	lmdat.int=data.frame(exp=thisGene.exp,Grouping=Grouping)

	## Try negative binomial glm with theta from final DESeq2 gene specific dispersion; that failing, fallback to glm.nb with automatic theta estimation; or that failing, use lm
	model <- tryCatch( glm(exp~Grouping,data=lmdat.int,family=MASS::negative.binomial(dispMat[this.gene,"dispersion"])), error= function(e) tryCatch( MASS::glm.nb(exp ~ Grouping, data = lmdat.int), error= function(e2) lm(exp~Grouping,data = lmdat.int) ))

	SE=sqrt(diag(vcov(model))[2])
	LFC=coef(model)[2]  #(statOut2 LFC)
	#ALT LFC (statOut LFC)
#	group1.counts.0to0.5<-lmdat.int$exp[which(Grouping==names(table(Grouping))[1])]
#	group2.counts.0to0.5<-lmdat.int$exp[which(Grouping==names(table(Grouping))[2])]
#	group1.counts.0to0.5[group1.counts.0to0.5==0]<- 0.5
#	group2.counts.0to0.5[group2.counts.0to0.5==0]<- 0.5
#	LFC=mean(log2(group1.counts.0to0.5))-mean(log2(group2.counts.0to0.5))
	stat=LFC/SE
	lowerTail=if(LFC<0) { TRUE } else { FALSE }
	pvalue=pnorm(LFC/SE,lower.tail=lowerTail)*2
	if(pvalue>1) pvalue=1

	return(c(Gene=this.gene,baseMean=mean(lmdat.int$exp),LFC=LFC,lfcSE=SE,stat=stat,pvalue=pvalue))
}
statOut2.rownames=statOut2[,"Gene"]
statOut2<-apply(statOut2[,2:ncol(statOut2)],2,as.numeric)
rownames(statOut2)<-statOut2.rownames
colnames(statOut2)<-as.data.frame(do.call(rbind,strsplit(colnames(statOut2),"[.]")))[,1]

statOut2<-cbind(statOut2, p.adjust(statOut2[,"pvalue"],method="BH"))
colnames(statOut2)[ncol(statOut2)]<-"padj"

#write.csv(statOut2,"statOut2.integer.csv")


## Cap SE at 3, larger than any value output in DESeq2  -- affects a handful of values
statOut<-as.data.frame(statOut)
statOut2<-as.data.frame(statOut2)
statOut$lfcSE[statOut$lfcSE>3]<-3
statOut2$lfcSE[statOut2$lfcSE>3]<-3



## precalculated DESeq2 Stats Data:
#nativeStat<-read.csv("C11_Astrocyte-CIBOP pulldown RNA vs Group2-CIBOP pulldown RNA_rowmeans_17038_Diffex.csv",header=TRUE,row.names=1)
# in-session calculations used:

## or use in-line calculation:
nativeStat<-as.data.frame(results(dds2))



pdf("[simData]Comparison of stats underlying P value calcs by NB GLM in DESeq2.pdf",width=12,height=18)
  par(mfrow=c(3,2))
  par(mar=c(5,6,3.5,2))
  WGCNA::verboseScatterplot(nativeStat$stat,as.data.frame(statOut)$stat*(-1),pch=21,bg="#FFFFFF33",col="darkmagenta",xlab="DESeq2 Native Z stat",ylab="NB glm Z stat(1),\nusing DESeq2 final gene-specific dispersion", main="Z stat calc using nonZero mean log2\n")
  WGCNA::verboseScatterplot(nativeStat$stat,as.data.frame(statOut2)$stat,pch=21,bg="#FFFFFF33",col="darkmagenta",xlab="DESeq2 Native Z stat",ylab="NB glm Z stat(2),\nusing DESeq2 final gene-specific dispersion", main="Z stat calc using coef(model)[Group2vsGroup1]\n")
  
  WGCNA::verboseScatterplot(nativeStat$log2FoldChange,as.data.frame(statOut)$LFC*(-1),pch=21,bg="#FFFFFF33",col="darkred",xlab="DESeq2 Native log2FC",ylab="mean(log2( nonZero (0.5 if 0) )) difference", main="LFC calc using nonZero mean log2\n")
  WGCNA::verboseScatterplot(nativeStat$log2FoldChange,as.data.frame(statOut2)$LFC,pch=21,bg="#FFFFFF33",col="darkred",xlab="DESeq2 Native log2FC",ylab="LFC as coef(model)[Group2vsGroup1]", main="LFC calc using coef(model)[Group2vsGroup1]\n")
  
  frame()
  WGCNA::verboseScatterplot(nativeStat$lfcSE,as.data.frame(statOut2)$lfcSE,pch=21,bg="#FFFFFF33",col="darkslateblue",xlab="DESeq2 Native Log2FC Std Error",ylab="sqrt(vcov(model))[Group2vsGroup1]", main="LFC SE calc using sqrt(vcov(model))[Group2vsGroup1]\n")
dev.off()



####################################
## Repeat stats with normalized continuous counts

# Pre-calculated:
#cleanDat.RNA.continuous<-read.csv("Exp. 9-2 & 10-1 DESeq2_norm counts_all groups_rowmeans_17038_10062023.csv",header=TRUE,row.names=1,check.names=FALSE)
#cleanDat.RNA.continuous<-cleanDat.RNA.continuous[,this.gene,which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]

# Calculated in place (6 columns only)
cleanDat.RNA.continuous<-as.data.frame(counts(dds2, normalized=TRUE))


dim(cleanDat.RNA.continuous)
#[1] 17038    6
cleanDat.RNA.continuous<-as.matrix(log2(cleanDat.RNA.continuous))
cleanDat.RNA.continuous[!is.finite(cleanDat.RNA.continuous)]<-0  #above operation set -Inf to 0 on row 63 and elsewhere
#cleanDat.RNA<-na.omit(cleanDat.RNA.continuous)
dim(cleanDat.RNA.continuous)
#[1] 17038    29




## Repeat stats table calculation with negative binomial glm used for both log2FC and lfcSE  -- USING CONTINUOUS NORMALIZED COUNTS DATA
statOut.continuous <-foreach(this.gene=as.character(rownames(cleanDat.RNA.continuous)), .combine='rbind') %dopar% {
	#this.gene="Gene1"
	
	thisGene.exp<-unlist(as.vector(cleanDat.RNA.continuous )) #[this.gene,which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]))
	Grouping=traits$group.simple[which(traits$group.simple=="Group1" | traits$group.simple=="Group2")]

	lmdat.int=data.frame(exp=thisGene.exp,Grouping=Grouping)

	## Try negative binomial glm with theta from final DESeq2 gene specific dispersion; that failing, fallback to glm.nb with automatic theta estimation; or that failing, use lm
	model <- tryCatch( glm(exp~Grouping,data=lmdat.int,family=MASS::negative.binomial(dispMat[this.gene,"dispersion"])), error= function(e) tryCatch( MASS::glm.nb(exp ~ Grouping, data = lmdat.int), error= function(e2) lm(exp~Grouping,data = lmdat.int) ))

	SE=sqrt(diag(vcov(model))[2])
	LFC=coef(model)[2]  #(statOut2 LFC)

	stat=LFC/SE
	lowerTail=if(LFC<0) { TRUE } else { FALSE }
	pvalue=pnorm(LFC/SE,lower.tail=lowerTail)*2
	if(pvalue>1) pvalue=1

	return(c(Gene=this.gene,baseMean=mean(lmdat.int$exp),LFC=LFC,lfcSE=SE,stat=stat,pvalue=pvalue))
}
statOut.continuous.rownames=statOut.continuous[,"Gene"]
statOut.continuous<-apply(statOut.continuous[,2:ncol(statOut.continuous)],2,as.numeric)
rownames(statOut.continuous)<-statOut.continuous.rownames
colnames(statOut.continuous)<-as.data.frame(do.call(rbind,strsplit(colnames(statOut.continuous),"[.]")))[,1]

statOut.continuous<-cbind(statOut2, p.adjust(statOut.continuous[,"pvalue"],method="BH"))
colnames(statOut.continuous)[ncol(statOut.continuous)]<-"padj"

#write.csv(statOut.continuous,"statOut.continuous.csv")


## Cap SE at 3, larger than any value output in DESeq2
statOut.continuous<-as.data.frame(statOut.continuous)
statOut.continuous$lfcSE[statOut.continuous$lfcSE>3]<-3

pdf("[simData]Comparison of stats underlying P value calcs by NB GLM in DESeq2[compareINTEGERvsCONTINUOUS_input].pdf",width=12,height=18)
  par(mfrow=c(3,2))
  par(mar=c(5,6,6,2))
  WGCNA::verboseScatterplot(nativeStat$stat,as.data.frame(statOut.continuous)$stat,pch=21,bg="#FFFFFF33",col="darkmagenta",xlab="DESeq2 Native Z stat",ylab="CONTINUOUS NORM DATA; NB glm Z stat(1),\nusing DESeq2 final gene-specific dispersion", main="Z stat via DESeq2 INTEGER (x)\nvs. NB glm CONTINUOUS (y)\n")
  WGCNA::verboseScatterplot(as.data.frame(statOut2)$stat,as.data.frame(statOut.continuous)$stat,pch=21,bg="#FFFFFF33",col="darkmagenta",xlab="NB glm INTEGER Z stat Calc",ylab="NB glm Continuous Z stat,\nusing DESeq2 final gene-specific dispersion", main="Z stat via NB glm INTEGER (x)vs. NB glm CONTINUOUS (y)\n")
  
  WGCNA::verboseScatterplot(nativeStat$log2FoldChange,as.data.frame(statOut.continuous)$LFC,pch=21,bg="#FFFFFF33",col="darkred",xlab="DESeq2 Native log2FC",ylab="mean(log2( nonZero (0.5 if 0) )) difference", main="LFC calc via DESeq2 INTEGER (x)\nvs. NB glm CONTINUOUS (y)\n")
  WGCNA::verboseScatterplot(as.data.frame(statOut2)$LFC,as.data.frame(statOut.continuous)$LFC,pch=21,bg="#FFFFFF33",col="darkred",xlab="NB glm INTEGER log2FC",ylab="NB glm Continuous log2FC", main="LFC calc NB glm INTEGER (x)\nvs. NB glm CONTINUOUS (y)\n")
  
  WGCNA::verboseScatterplot(nativeStat$lfcSE,as.data.frame(statOut.continuous)$lfcSE,pch=21,bg="#FFFFFF33",col="darkslateblue",xlab="DESeq2 Native Log2FC Std Error",ylab="NB glm Continuous Log2FC Std Error", main="LFC SE calc via DESeq2 INTEGER (x)\nvs. NB glm CONTINUOUS (y)\n")
  WGCNA::verboseScatterplot(as.data.frame(statOut2)$lfcSE,as.data.frame(statOut.continuous)$lfcSE,pch=21,bg="#FFFFFF33",col="darkslateblue",xlab="NB glm INTEGER Log2FC Std Error",ylab="NB glm Continuous Log2FC Std Error", main="LFC SE calc via NB glm INTEGER (x)\nvs. NB glm CONTINUOUS (y)\n")
dev.off()
