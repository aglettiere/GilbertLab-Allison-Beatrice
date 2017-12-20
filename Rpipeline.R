#===================================================================================
#
#		MICROBIAL COMMUNITY STUDY FROM 16s rRNA
#
#===================================================================================
# Takes biom, map file and tree from QIIME and use the phyloseq package to obtain
# alpha-diversity, beta-diversity, heatmaps of share OTUs and analyze some results
# with DESeq and XgBoost

# By Beatriz Penalver Bernabe 06/14/2016
#===================================================================================
#===================================================================================


#===================================================================================
#===================================================================================

# Upload all the required libraries
library(phyloseq)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(viridis)
library(vegan)
library("DESeq2")
library(xgboost)
library('data.table')
library(Matrix)
library(readr)
library(stringr)
library(caret)
library(car)
library(gplots)
library(RColorBrewer)
library(limma)
library(metagenomeSeq)
library(plyr)
library(statmod)
library(randomForest)
library(biomformat)

# Default preference for ggplots
#===================================================================================
theme_set(theme_bw())

#Import data from QIIME
#===================================================================================
#Import biom file
biomfile = "/Users/beatrizp/Dropbox/SCFA_Obesity/otu_table_openref.json" # It is a json file
biom = import_biom(biomfile)

#Import tree
tree=read_tree_greengenes('/Users/beatrizp/Dropbox/SCFA_Obesity/OTUclustering_openreference_usearch/rep_set.tre')

#Import mapfile
mapfile = import_qiime_sample_data("/Users/beatrizp/Dropbox/SCFA_Obesity/Dugas_Demo_mapping_completev5.txt")
	 
#Read the taxonomy file
tax <- read.table(file='/Users/beatrizp/Dropbox/SCFA_Obesity/OTUclustering_openreference_usearch/uclust_assigned_taxonomy/rep_set_tax_assignments.txt', sep="\t",stringsAsFactors=FALSE)
	
#Split the file so it can serve as input for phyloseq
taxonomy_split<-NULL
for (itax in 1:nrow(tax)){
	ids<-strsplit(tax[itax,2],";",fixed=TRUE)[[1]]
	vec_ids<-c()
	for (ii in 1:length(ids)){
		vec_ids<-c(vec_ids, strsplit(ids[ii],"__",fixed=TRUE)[[1]][2])
		print(vec_ids)
	}
	if (length(vec_ids)<7){
		vec_ids<-c(vec_ids, rep(NA,(7-length(vec_ids))))
	}
	if (itax==1){
		taxonomy_split<-vec_ids
	} else{
		taxonomy_split<-rbind(taxonomy_split,vec_ids)
	}
}
rownames(taxonomy_split)<-tax[,1]
colnames(taxonomy_split)<-c("Kingdom","Phylum","Class", "Order","Family","Genus","Species") 
# tax_split<-tax_split[na.omit(match(rownames(biom),
	# rownames(tax_split))),]
tax_split<-tax_table(taxonomy_split)
#Merge all the files
final = merge_phyloseq(biom, mapfile,tax_split,tree)


#Heatmap shared OTUs
#===================================================================================
# Get otu_table
otu_table_final<-otu_table(final)

# Removed samples that have less reads than the water control
otu_table_final_filtersamples<-otu_table_final[,colSums(otu_table_final)>500]

#Remove sequences per sample that are really close to the sequencing thresolhd
otu_table_final_filtersamples_OTUs<-otu_table_final_filtersamples[
	apply(otu_table_final_filtersamples,	1,function(x){
		rr<-ifelse(all(as.numeric(x)<10), FALSE,TRUE)
		return(rr)
		}),]
		
#Remove those that are present at a very low abundance
otu_table_final_filtersamples_OTUs<-
	otu_table_final_filtersamples[rowSums(otu_table_final_filtersamples)/sum(rowSums(otu_table_final_filtersamples))*100>0.1,]
	
otu_table_final_filtersamples_OTUs<-otu_table(otu_table_final_filtersamples_OTUs)


#Generate the new pyloseq file for normalization 
final_filtered = merge_phyloseq(otu_table_final_filtersamples_OTUs, mapfile,tax_split,tree)

# #Determine the shared OTUs for each sample pair
# shared_otus<-matrix(0,ncol=ncol(otu_table_final),nrow=ncol(otu_table_final))
# for (isample in 1:ncol(otu_table_final)){
	# for (jsample in 1:ncol(otu_table_final)){
		# minOTUS<-length(which(otu_table_final[,isample]*otu_table_final[,jsample]>0))
		# shared_otus[isample,jsample]<-minOTUS
	# }
# }
# colnames(shared_otus)<-rownames(shared_otus)<-
	# mapfile$N1_gh_lean_N2_gh_obese_N3_us_lean_N4_us_obese[match(colnames(otu_table_final),
		# rownames(mapfile))]
	
# # Calculate the distance between samples based on their shared OTUs
# distance_otus<-dist(shared_otus,method="euclidean") #distance matrix

# #Clasify them using hierchical clustering
# fit <-hclust(distance_otus,method="ward.D")
# plot(fit)

# #Identified the number of groups
# groups<-cutree(fit,k=4)

# #Re-arrange the files
# shared_otus_order<-shared_otus[order(groups),order(groups)]
# rownames(shared_otus_order)<-rownames(shared_otus)[order(groups)]
# colnames(shared_otus_order)<-rownames(shared_otus)[order(groups)]

# #Plot heat maps using gplots
# colheat=c(brewer.pal (8,'Blues'))
# pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/heatmap_sharedOTUS.pdf')
# heatmap.2(shared_otus_order, Rowv=FALSE, Colv=FALSE, 
	# dendrogram="none", col=colheat,key=T, keysize=1.5, density.info="none", 
    # trace="none", labRow=TRUE, breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4))
# dev.off()

#Abundance plots as a function of the different meta-data variables
#===================================================================================

#Create a directory to save the abundance plots
dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref")

#Order data

sample_data_order<-data.frame(sample_data(final_filtered),stringsAsFactors=FALSE)

sample_data_order<-sample_data_order[order(sample_data_order[,"Country"],sample_data_order[,"BMI"]),]
sample_data_order_phy<-sample_data(sample_data_order)
final_order<-merge_phyloseq(otu_table(final_filtered), sample_data_order_phy,tax_split,tree)
bb<-psmelt(final_order)
cc=merge(bb,sample_data_order[,c("X.SampleID","BMI","Country","obese","AGE","Percentbodyfatmass","formicacid",
	"aceticacid","propionicacid","butyricacid","isovalericacid","GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")],all.x=TRUE,by.x="Sample",by.y="X.SampleID")

cc[,"BMI"]=cc[,"BMI.x"]
cc[,"Country"]=cc[,"Country.x"]
cc[,"obese"]=cc[,"obese.x"]
cc[,"AGE"]=cc[,"AGE.x"]
cc[,"Percentbodyfatmass"]=cc[,"Percentbodyfatmass.x"]
cc[,"formicacid"]=cc[,"formicacid.x"]
cc[,"aceticacid"]=cc[,"aceticacid.x"]
cc[,"propionicacid"]=cc[,"propionicacid.x"]
cc[,"butyricacid"]=cc[,"butyricacid.x"]
cc[,"isovalericacid"]=cc[,"isovalericacid.x"]
cc[,"adiponectin"]=cc[,"adiponectin.x"]
cc[,"GlucoseResultGLUCRSLT"]=cc[,"GlucoseResultGLUCRSLT.x"]
cc[,"leptin"]=cc[,"leptin.x"]
cc[,"insulin"]=cc[,"insulin.x"]
cc[,"glucosemmol"]=cc[,"glucosemmol.x"]
cc[,"insulinmmol"]=cc[,"insulinmmol.x"]
cc[,"homa_ir"]=cc[,"homa_ir.x"]
cc[,"ratio"]=cc[,"ratio.x"]
dd<-cc[order(cc[,"Country"],cc[,"BMI"]),]

# Phylum level
#--------------------------------------------------------------------------------------------
nphylum<-sort(tapply(dd$Abundance,dd$Phylum,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
phyum_plot_data<-NULL
for (isample in nsample){
	for (iphylum in names(nphylum)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_phylum<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_iphylum<-isample_melt[which(isample_melt$Phylum==iphylum),]
		count_iphylum<-sum(as.numeric(isample_iphylum[,"Abundance"]))
		abundance_iphylum<-sum(as.numeric(isample_iphylum[,"Abundance"]))/total_phylum*100
		if (is.null(phyum_plot_data)){
			phyum_plot_data<-c(isample,iphylum,count_iphylum,abundance_iphylum,
				isample_melt[which(isample_melt$Phylum==iphylum),"BMI"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"Country"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"obese"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"AGE"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"formicacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"aceticacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"propionicacid"][1],							
				isample_melt[which(isample_melt$Phylum==iphylum),"butyricacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"isovalericacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"adiponectin"][1],			
				isample_melt[which(isample_melt$Phylum==iphylum),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"leptin"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"insulin"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"glucosemmol"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"insulinmmol"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"homa_ir"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"ratio"][1])
		} else {
			phyum_plot_data<-rbind(phyum_plot_data,
				c(isample,iphylum,count_iphylum,abundance_iphylum,
				isample_melt[which(isample_melt$Phylum==iphylum),"BMI"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"Country"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"obese"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"AGE"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"formicacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"aceticacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"propionicacid"][1],							
				isample_melt[which(isample_melt$Phylum==iphylum),"butyricacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"isovalericacid"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"adiponectin"][1],			
				isample_melt[which(isample_melt$Phylum==iphylum),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"leptin"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"insulin"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"glucosemmol"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"insulinmmol"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"homa_ir"][1],
				isample_melt[which(isample_melt$Phylum==iphylum),"ratio"][1]))
		}
	}
}
colnames(phyum_plot_data)<- c("Sample","Phylum","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(phyum_plot_data[,"BMI"])-min(as.numeric(phyum_plot_data[,"BMI"])))/
	(max(as.numeric(phyum_plot_data[,"BMI"]))-min(as.numeric(phyum_plot_data[,"BMI"])))*100
phyum_plot_data<-data.frame(cbind(phyum_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(phyum_plot_data)<- c("Sample","Phylum","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Phylum=phyum_plot_data$Phylum[1:length(unique(phyum_plot_data$Phylum))],
	Phylum_dummy=phyum_plot_data$Phylum[1:length(unique(phyum_plot_data$Phylum))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(phyum_plot_data$Counts),nrow=length(unique(phyum_plot_data$Phylum)),byrow=FALSE)
rownames(matrixData)<-phyum_plot_data$Phylum[1:length(unique(phyum_plot_data$Phylum))]
colnames(matrixData)<-phyum_plot_data$Sample[seq(1,nrow(phyum_plot_data),
	by=length(unique(phyum_plot_data$Phylum)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
phylum_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(sample_data(metadata)), 
	featureData = AnnotatedDataFrame(featureData))
phylum_metagenomeSeq_filter = filterData(phylum_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(phylum_metagenomeSeq_filter )
phylum_metagenomeSeq_filter <- cumNorm(phylum_metagenomeSeq_filter, p =p) 
phylum_metagenomeSeq_filter_nor = MRcounts(phylum_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(phylum_metagenomeSeq_filter))
pd <- pData(phylum_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(phylum_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(phylum_metagenomeSeq_filter)$obese
AGE = pData(phylum_metagenomeSeq_filter)$AGE 
normFactor = normFactors(phylum_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = phylum_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(phylum_metagenomeSeq_filter)$counts))
  
obese = pData(phylum_metagenomeSeq_filter)$obese
AGE = pData(phylum_metagenomeSeq_filter)$AGE 
Country= pData(phylum_metagenomeSeq_filter)$Country
normFactor = normFactors(phylum_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = phylum_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(phylum_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(phylum_metagenomeSeq_filter)$obese,
	pData(phylum_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age")
# fitting the ZIG model 
res = fitZig(obj = phylum_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	(obesity1_USA + obesity1_Ghana+obesity0_USA)/3-obesity0_Ghana,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_phyum<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
phylum_mG<-cbind("TaxonomicLevel"=rep("Phylum",nrow(pvalues_phyum)),
	 fit2$coefficients, pvalues_phyum)
rownames(phylum_mG)<-rownames(fit2$coefficients)

# Fitting linear models with limma
dummy_matrix<-matrix(0,ncol=length(unique(phyum_plot_data[,"Phylum"])), 
	nrow=length(unique(phyum_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(phyum_plot_data[counter,"Phylum"]==unique(phyum_plot_data[,"Phylum"])[icol],
				phyum_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(phyum_plot_data[,"Phylum"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(phyum_plot_data[,"Phylum"])
rownames(dummy_matrix)<-unique(phyum_plot_data[,"Sample"])
  
design_data<-phyum_plot_data[!duplicated(phyum_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
phyum_matrix<-t(dummy_matrix)
fit_phyum<-lmFit(phyum_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_phyum<-contrasts.fit(fit_phyum, cont.matrix)
fit2_phyum<-eBayes(fit2_phyum)
pvalues_phyum<-apply(fit2_phyum$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
phylum_limma<-cbind("TaxonomicLevel"=rep("Phylum",nrow(pvalues_phyum)),
	 fit2_phyum$coefficients,pvalues_phyum)
rownames(phylum_limma)<-rownames(fit2_phyum$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
data_deseq2<-DESeqDataSetFromMatrix(matrixData, 
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="mean")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
phylum_DESeq2<-cbind("TaxonomicLevel"=rep("Phylum", nrow(res_final)), res_final)
rownames(phylum_DESeq2)<-rownames(matrixData)

phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

phyum_plot_data[,"Abundance"]<-as.numeric(phyum_plot_data[,"Abundance"])
phyum_plot_data[,"BMInorm"]<-as.numeric(phyum_plot_data[,"BMInorm"])
phyum_plot_data$Sample <- factor(phyum_plot_data$Sample, levels = phyum_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/PhylumvsSample.pdf')              
p = ggplot(phyum_plot_data, aes(x = Sample, y = Abundance, fill = Phylum, show.legend=FALSE))
p = p + geom_bar(aes(color=Phylum, fill=Phylum),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") +scale_fill_manual (values=phy)+ 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
	geom_point(data=phyum_plot_data, aes(x = Sample, y = as.numeric(BMInorm),show.legend=FALSE),colour="black",
	shape = 2,  fill = "white", size = 2) 
print(p)
dev.off()

phyum_plot_data[,"Abundance"]<-as.numeric(phyum_plot_data[,"Abundance"])
phyum_plot_data$Sample <- factor(phyum_plot_data$Sample, levels = phyum_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/PhylumvsGroup_facetwrap.pdf')              
p = ggplot(phyum_plot_data, aes(x = Sample, y = Abundance, fill = Phylum, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) +
	facet_wrap(~Phylum, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/SamplevsBMI.pdf')          
 minBMI=min(as.numeric(phyum_plot_data$BMI),na.rm=TRUE)
 maxBMI=max(as.numeric(phyum_plot_data$BMI),na.rm=TRUE)
limits_BMI=seq(minBMI,maxBMI,by=5)
p1=ggplot(data=phyum_plot_data, aes(x = Sample, y = as.numeric(BMI),colour=Country,show.legend=FALSE))+geom_point()+ 	
	theme_bw() + ylab("BMI") +scale_fill_manual (values=c("red","blue"))+ 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +theme(text = element_text(size = 12)) + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
print(p1)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/PhylumvsObese_facetwrap.pdf')              
p = ggplot(phyum_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Phylum, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()

# Class level
#--------------------------------------------------------------------------------------------------------------------------------------------
nclass<-sort(tapply(dd$Abundance,dd$Class,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
class_plot_data<-NULL
for (isample in nsample){
	for (iclass in names(nclass)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_class<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_iclass<-isample_melt[which(isample_melt$Class==iclass),]
		count_iclass<-sum(as.numeric(isample_iclass[,"Abundance"]))
		abundance_iclass<-sum(as.numeric(isample_iclass[,"Abundance"]))/total_class*100
		if (is.null(class_plot_data)){
			class_plot_data<-c(isample,iclass,count_iclass,abundance_iclass,
				isample_melt[which(isample_melt$Class==iclass),"BMI"][1],
				isample_melt[which(isample_melt$Class==iclass),"Country"][1],
				isample_melt[which(isample_melt$Class==iclass),"obese"][1],
				isample_melt[which(isample_melt$Class==iclass),"AGE"][1],
				isample_melt[which(isample_melt$Class==iclass),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Class==iclass),"formicacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"aceticacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"propionicacid"][1],							
				isample_melt[which(isample_melt$Class==iclass),"butyricacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"isovalericacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"adiponectin"][1],			
				isample_melt[which(isample_melt$Class==iclass),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Class==iclass),"leptin"][1],
				isample_melt[which(isample_melt$Class==iclass),"insulin"][1],
				isample_melt[which(isample_melt$Class==iclass),"glucosemmol"][1],
				isample_melt[which(isample_melt$Class==iclass),"insulinmmol"][1],
				isample_melt[which(isample_melt$Class==iclass),"homa_ir"][1],
				isample_melt[which(isample_melt$Class==iclass),"ratio"][1])
		} else {
			class_plot_data<-rbind(class_plot_data,
				c(isample,iclass,count_iclass,abundance_iclass,
				isample_melt[which(isample_melt$Class==iclass),"BMI"][1],
				isample_melt[which(isample_melt$Class==iclass),"Country"][1],
				isample_melt[which(isample_melt$Class==iclass),"obese"][1],
				isample_melt[which(isample_melt$Class==iclass),"AGE"][1],
				isample_melt[which(isample_melt$Class==iclass),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Class==iclass),"formicacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"aceticacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"propionicacid"][1],							
				isample_melt[which(isample_melt$Class==iclass),"butyricacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"isovalericacid"][1],
				isample_melt[which(isample_melt$Class==iclass),"adiponectin"][1],			
				isample_melt[which(isample_melt$Class==iclass),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Class==iclass),"leptin"][1],
				isample_melt[which(isample_melt$Class==iclass),"insulin"][1],
				isample_melt[which(isample_melt$Class==iclass),"glucosemmol"][1],
				isample_melt[which(isample_melt$Class==iclass),"insulinmmol"][1],
				isample_melt[which(isample_melt$Class==iclass),"homa_ir"][1],
				isample_melt[which(isample_melt$Class==iclass),"ratio"][1]))
		}
	}
}
colnames(class_plot_data)<- c("Sample","Class","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(class_plot_data[,"BMI"])-min(as.numeric(class_plot_data[,"BMI"])))/
	(max(as.numeric(class_plot_data[,"BMI"]))-min(as.numeric(class_plot_data[,"BMI"])))*100
class_plot_data<-data.frame(cbind(class_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(class_plot_data)<- c("Sample","Class","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Class=class_plot_data$Class[1:length(unique(class_plot_data$Class))],
	Class_dummy=class_plot_data$Class[1:length(unique(class_plot_data$Class))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(class_plot_data$Counts),nrow=length(unique(class_plot_data$Class)),byrow=FALSE)
rownames(matrixData)<-class_plot_data$Class[1:length(unique(class_plot_data$Class))]
colnames(matrixData)<-class_plot_data$Sample[seq(1,nrow(class_plot_data),by=length(unique(class_plot_data$Class)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
class_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
class_metagenomeSeq_filter = filterData(class_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(class_metagenomeSeq_filter )
class_metagenomeSeq_filter <- cumNorm(class_metagenomeSeq_filter, p =p) 
class_metagenomeSeq_filter_nor = MRcounts(class_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(class_metagenomeSeq_filter))
pd <- pData(class_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(class_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(class_metagenomeSeq_filter)$obese
AGE = pData(class_metagenomeSeq_filter)$AGE 
normFactor = normFactors(class_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = class_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(class_metagenomeSeq_filter)$counts))
  
obese = pData(class_metagenomeSeq_filter)$obese
AGE = pData(class_metagenomeSeq_filter)$AGE 
Country= pData(class_metagenomeSeq_filter)$Country
normFactor = normFactors(class_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = class_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(class_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(class_metagenomeSeq_filter)$obese,
	pData(class_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age")
# fitting the ZIG model 
res = fitZig(obj = class_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	(obesity1_USA + obesity1_Ghana+obesity0_USA)/3-obesity0_Ghana,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_class<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
class_mG<-cbind("TaxonomicLevel"=rep("Class",nrow(pvalues_class)),
	 fit2$coefficients,pvalues_class)
rownames(class_mG)<-rownames(fit2$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
matrixData<-matrix(as.numeric(class_plot_data$Counts),nrow=length(unique(class_plot_data$Class)),byrow=FALSE)
rownames(matrixData)<-class_plot_data$Class[1:length(unique(class_plot_data$Class))]
colnames(matrixData)<-class_plot_data$Sample[seq(1,nrow(class_plot_data),by=length(unique(class_plot_data$Class)))]
data_deseq2<-DESeqDataSetFromMatrix(matrixData,
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="mean")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
class_DESeq2<-cbind("TaxonomicLevel"=rep("Class", nrow(res_final)), res_final)
 
# Fitting models with limma
dummy_matrix<-matrix(0,ncol=length(unique(class_plot_data[,"Class"])), 
	nrow=length(unique(class_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(class_plot_data[counter,"Class"]==unique(class_plot_data[,"Class"])[icol],
				class_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(class_plot_data[,"Class"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(class_plot_data[,"Class"])
rownames(dummy_matrix)<-unique(class_plot_data[,"Sample"])
  
design_data<-class_plot_data[!duplicated(class_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
class_matrix<-t(dummy_matrix)
fit_class<-lmFit(class_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_class<-contrasts.fit(fit_class, cont.matrix)
fit2_class<-eBayes(fit2_class)
pvalues_class<-apply(fit2_class$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
class_limma<-cbind("TaxonomicLevel"=rep("Class",nrow(pvalues_class)),
	 fit2_class$coefficients,pvalues_class)
rownames(class_limma)<-rownames(fit2_class$coefficients)


phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

class_plot_data[,"Abundance"]<-as.numeric(class_plot_data[,"Abundance"])
class_plot_data[,"BMInorm"]<-as.numeric(class_plot_data[,"BMInorm"])
class_plot_data$Sample <- factor(class_plot_data$Sample, levels = class_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/ClassvsSample.pdf')              
p = ggplot(class_plot_data, aes(x = Sample, y = Abundance, fill = Class, show.legend=FALSE))
p = p + geom_bar(aes(color=Class, fill=Class),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") + 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
	geom_point(data=class_plot_data, aes(x = Sample, y = as.numeric(BMInorm),show.legend=FALSE),colour="black",
	shape = 2,  fill = "white", size = 2) 
print(p)
dev.off()

class_plot_data[,"Abundance"]<-as.numeric(class_plot_data[,"Abundance"])
class_plot_data$Sample <- factor(class_plot_data$Sample, levels = class_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/ClassvsGroup_facetwrap.pdf')              
p = ggplot(class_plot_data, aes(x = Sample, y = Abundance, fill = Class, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) +
	facet_wrap(~Class, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/SamplevsBMI.pdf')              
limits_BMI=seq(minBMI,maxBMI,by=5)
p1=ggplot(data=class_plot_data, aes(x = Sample, y = as.numeric(BMI),colour=Country,show.legend=FALSE))+geom_point()+ 	
	theme_bw() + ylab("BMI") +scale_fill_manual (values=c("red","blue"))+ 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +theme(text = element_text(size = 12)) + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
print(p1)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/ClassvsObese_facetwrap.pdf')              
p = ggplot(class_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Class, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()

# For the order level
#==================================================================================================================
norder<-sort(tapply(dd$Abundance,dd$Order,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
order_plot_data<-NULL
for (isample in nsample){
	for (iorder in names(norder)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_order<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_iorder<-isample_melt[which(isample_melt$Order==iorder),]
		count_iorder<-sum(as.numeric(isample_iorder[,"Abundance"]))
		abundance_iorder<-sum(as.numeric(isample_iorder[,"Abundance"]))/total_order*100
		if (is.null(order_plot_data)){
			order_plot_data<-c(isample,iorder,count_iorder,abundance_iorder,
				isample_melt[which(isample_melt$Order==iorder),"BMI"][1],
				isample_melt[which(isample_melt$Order==iorder),"Country"][1],
				isample_melt[which(isample_melt$Order==iorder),"obese"][1],
				isample_melt[which(isample_melt$Order==iorder),"AGE"][1],
				isample_melt[which(isample_melt$Order==iorder),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Order==iorder),"formicacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"aceticacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"propionicacid"][1],							
				isample_melt[which(isample_melt$Order==iorder),"butyricacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"isovalericacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"adiponectin"][1],			
				isample_melt[which(isample_melt$Order==iorder),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Order==iorder),"leptin"][1],
				isample_melt[which(isample_melt$Order==iorder),"insulin"][1],
				isample_melt[which(isample_melt$Order==iorder),"glucosemmol"][1],
				isample_melt[which(isample_melt$Order==iorder),"insulinmmol"][1],
				isample_melt[which(isample_melt$Order==iorder),"homa_ir"][1],
				isample_melt[which(isample_melt$Order==iorder),"ratio"][1])
		} else {
			order_plot_data<-rbind(order_plot_data,
				c(isample,iorder,count_iorder,abundance_iorder,
				isample_melt[which(isample_melt$Order==iorder),"BMI"][1],
				isample_melt[which(isample_melt$Order==iorder),"Country"][1],
				isample_melt[which(isample_melt$Order==iorder),"obese"][1],
				isample_melt[which(isample_melt$Order==iorder),"AGE"][1],
				isample_melt[which(isample_melt$Order==iorder),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Order==iorder),"formicacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"aceticacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"propionicacid"][1],							
				isample_melt[which(isample_melt$Order==iorder),"butyricacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"isovalericacid"][1],
				isample_melt[which(isample_melt$Order==iorder),"adiponectin"][1],			
				isample_melt[which(isample_melt$Order==iorder),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Order==iorder),"leptin"][1],
				isample_melt[which(isample_melt$Order==iorder),"insulin"][1],
				isample_melt[which(isample_melt$Order==iorder),"glucosemmol"][1],
				isample_melt[which(isample_melt$Order==iorder),"insulinmmol"][1],
				isample_melt[which(isample_melt$Order==iorder),"homa_ir"][1],
				isample_melt[which(isample_melt$Order==iorder),"ratio"][1]))
		}
	}
}
colnames(order_plot_data)<- c("Sample","Order","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(order_plot_data[,"BMI"])-min(as.numeric(order_plot_data[,"BMI"])))/
	(max(as.numeric(order_plot_data[,"BMI"]))-min(as.numeric(order_plot_data[,"BMI"])))*100
order_plot_data<-data.frame(cbind(order_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(order_plot_data)<- c("Sample","Order","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Order=order_plot_data$Order[1:length(unique(order_plot_data$Order))],
	Order_dummy=order_plot_data$Order[1:length(unique(order_plot_data$Order))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(order_plot_data$Counts),nrow=length(unique(order_plot_data$Order)),byrow=FALSE)
rownames(matrixData)<-order_plot_data$Order[1:length(unique(order_plot_data$Order))]
colnames(matrixData)<-order_plot_data$Sample[seq(1,nrow(order_plot_data),by=length(unique(order_plot_data$Order)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
order_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
order_metagenomeSeq_filter = filterData(order_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(order_metagenomeSeq_filter )
order_metagenomeSeq_filter <- cumNorm(order_metagenomeSeq_filter, p =p) 
order_metagenomeSeq_filter_nor = MRcounts(order_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(order_metagenomeSeq_filter))
pd <- pData(order_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(order_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(order_metagenomeSeq_filter)$obese
AGE = pData(order_metagenomeSeq_filter)$AGE 
normFactor = normFactors(order_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = order_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(order_metagenomeSeq_filter)$counts))
  
obese = pData(order_metagenomeSeq_filter)$obese
AGE = pData(order_metagenomeSeq_filter)$AGE 
Country= pData(order_metagenomeSeq_filter)$Country
normFactor = normFactors(order_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = order_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(order_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(order_metagenomeSeq_filter)$obese,
	pData(order_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age")
# fitting the ZIG model 
res = fitZig(obj = order_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	(obesity1_USA + obesity1_Ghana+obesity0_USA)/3 -  obesity0_Ghana,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_order<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
order_mG<-cbind("TaxonomicLevel"=rep("Order",nrow(pvalues_order)),
	 fit2$coefficients,pvalues_order)
rownames(order_mG)<-rownames(fit2$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
matrixData<-matrix(as.numeric(order_plot_data$Counts),nrow=length(unique(order_plot_data$Order)),byrow=FALSE)
rownames(matrixData)<-order_plot_data$Order[1:length(unique(order_plot_data$Order))]
colnames(matrixData)<-order_plot_data$Sample[seq(1,nrow(order_plot_data),by=length(unique(order_plot_data$Order)))]
data_deseq2<-DESeqDataSetFromMatrix(matrixData,
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="mean")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
order_DESeq2<-cbind("TaxonomicLevel"=rep("Order", nrow(res_final)), res_final)
rownames(order_DESeq2)<-rownames(matrixData)
 
# Fitting models with limma
dummy_matrix<-matrix(0,ncol=length(unique(order_plot_data[,"Order"])), 
	nrow=length(unique(order_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(order_plot_data[counter,"Order"]==unique(order_plot_data[,"Order"])[icol],
				order_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(order_plot_data[,"Order"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(order_plot_data[,"Order"])
rownames(dummy_matrix)<-unique(order_plot_data[,"Sample"])
  
design_data<-order_plot_data[!duplicated(order_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
order_matrix<-t(dummy_matrix)
fit_order<-lmFit(order_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_order<-contrasts.fit(fit_order, cont.matrix)
fit2_order<-eBayes(fit2_order)
pvalues_order<-apply(fit2_order$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
order_limma<-cbind("TaxonomicLevel"=rep("Order",nrow(pvalues_order)),
	 fit2_order$coefficients,pvalues_order)
rownames(order_limma)<-rownames(fit2_order$coefficients)


phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

order_plot_data[,"Abundance"]<-as.numeric(order_plot_data[,"Abundance"])
order_plot_data[,"BMInorm"]<-as.numeric(order_plot_data[,"BMInorm"])
order_plot_data$Sample <- factor(order_plot_data$Sample, levels = order_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/OrdervsSample.pdf')              
p = ggplot(order_plot_data, aes(x = Sample, y = Abundance, fill = Order, show.legend=FALSE))
p = p + geom_bar(aes(color=Order, fill=Order),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") + 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
	geom_point(data=order_plot_data, aes(x = Sample, y = as.numeric(BMInorm),show.legend=FALSE),colour="black",
	shape = 2,  fill = "white", size = 2) 
print(p)
dev.off()

order_plot_data[,"Abundance"]<-as.numeric(order_plot_data[,"Abundance"])
order_plot_data$Sample <- factor(order_plot_data$Sample, levels = order_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/OrdervsGroup_facetwrap.pdf')              
p = ggplot(order_plot_data, aes(x = Sample, y = Abundance, fill = Order, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) +
	facet_wrap(~Order, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/OrdervsObese_facetwrap.pdf')              
p = ggplot(order_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Order, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()

# For the family level
#==================================================================================================================
nfamily<-sort(tapply(dd$Abundance,dd$Family,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
family_plot_data<-NULL
for (isample in nsample){
	for (ifamily in names(nfamily)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_family<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_ifamily<-isample_melt[which(isample_melt$Family==ifamily),]
		count_ifamily<-sum(as.numeric(isample_ifamily[,"Abundance"]))
		abundance_ifamily<-sum(as.numeric(isample_ifamily[,"Abundance"]))/total_family*100
		if (is.null(family_plot_data)){
			family_plot_data<-c(isample,ifamily,count_ifamily,abundance_ifamily,
				isample_melt[which(isample_melt$Family==ifamily),"BMI"][1],
				isample_melt[which(isample_melt$Family==ifamily),"Country"][1],
				isample_melt[which(isample_melt$Family==ifamily),"obese"][1],
				isample_melt[which(isample_melt$Family==ifamily),"AGE"][1],
				isample_melt[which(isample_melt$Family==ifamily),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Family==ifamily),"formicacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"aceticacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"propionicacid"][1],							
				isample_melt[which(isample_melt$Family==ifamily),"butyricacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"isovalericacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"adiponectin"][1],			
				isample_melt[which(isample_melt$Family==ifamily),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Family==ifamily),"leptin"][1],
				isample_melt[which(isample_melt$Family==ifamily),"insulin"][1],
				isample_melt[which(isample_melt$Family==ifamily),"glucosemmol"][1],
				isample_melt[which(isample_melt$Family==ifamily),"insulinmmol"][1],
				isample_melt[which(isample_melt$Family==ifamily),"homa_ir"][1],
				isample_melt[which(isample_melt$Family==ifamily),"ratio"][1])
		} else {
			family_plot_data<-rbind(family_plot_data,
				c(isample,ifamily,count_ifamily,abundance_ifamily,
				isample_melt[which(isample_melt$Family==ifamily),"BMI"][1],
				isample_melt[which(isample_melt$Family==ifamily),"Country"][1],
				isample_melt[which(isample_melt$Family==ifamily),"obese"][1],
				isample_melt[which(isample_melt$Family==ifamily),"AGE"][1],
				isample_melt[which(isample_melt$Family==ifamily),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Family==ifamily),"formicacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"aceticacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"propionicacid"][1],							
				isample_melt[which(isample_melt$Family==ifamily),"butyricacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"isovalericacid"][1],
				isample_melt[which(isample_melt$Family==ifamily),"adiponectin"][1],			
				isample_melt[which(isample_melt$Family==ifamily),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Family==ifamily),"leptin"][1],
				isample_melt[which(isample_melt$Family==ifamily),"insulin"][1],
				isample_melt[which(isample_melt$Family==ifamily),"glucosemmol"][1],
				isample_melt[which(isample_melt$Family==ifamily),"insulinmmol"][1],
				isample_melt[which(isample_melt$Family==ifamily),"homa_ir"][1],
				isample_melt[which(isample_melt$Family==ifamily),"ratio"][1]))
		}
	}
}
colnames(family_plot_data)<- c("Sample","Family","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(family_plot_data[,"BMI"])-min(as.numeric(family_plot_data[,"BMI"])))/
	(max(as.numeric(family_plot_data[,"BMI"]))-min(as.numeric(family_plot_data[,"BMI"])))*100
family_plot_data<-data.frame(cbind(family_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(family_plot_data)<- c("Sample","Family","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Family=family_plot_data$Family[1:length(unique(family_plot_data$Family))],
	Family_dummy=family_plot_data$Family[1:length(unique(family_plot_data$Family))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(family_plot_data$Counts),nrow=length(unique(family_plot_data$Family)),byrow=FALSE)
rownames(matrixData)<-family_plot_data$Family[1:length(unique(family_plot_data$Family))]
colnames(matrixData)<-family_plot_data$Sample[seq(1,nrow(family_plot_data),by=length(unique(family_plot_data$Family)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
family_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
family_metagenomeSeq_filter = filterData(family_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(family_metagenomeSeq_filter )
family_metagenomeSeq_filter <- cumNorm(family_metagenomeSeq_filter, p =p) 
family_metagenomeSeq_filter_nor = MRcounts(family_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(family_metagenomeSeq_filter))
pd <- pData(family_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(family_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(family_metagenomeSeq_filter)$obese
AGE = pData(family_metagenomeSeq_filter)$AGE 
normFactor = normFactors(family_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = family_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(family_metagenomeSeq_filter)$counts))
  
obese = pData(family_metagenomeSeq_filter)$obese
AGE = pData(family_metagenomeSeq_filter)$AGE 
Country= pData(family_metagenomeSeq_filter)$Country
normFactor = normFactors(family_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = family_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(family_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(family_metagenomeSeq_filter)$obese,
	pData(family_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age")
# fitting the ZIG model 
res = fitZig(obj = family_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	(obesity1_USA + obesity1_Ghana+obesity0_USA)/3- obesity0_Ghana,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_family<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
family_mG<-cbind("TaxonomicLevel"=rep("Family",nrow(pvalues_family)),
	 fit2$coefficients,pvalues_family)
rownames(family_mG)<-rownames(fit2$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
matrixData<-matrix(as.numeric(family_plot_data$Counts),nrow=length(unique(family_plot_data$Family)),byrow=FALSE)
rownames(matrixData)<-family_plot_data$Family[1:length(unique(family_plot_data$Family))]
colnames(matrixData)<-family_plot_data$Sample[seq(1,nrow(family_plot_data),by=length(unique(family_plot_data$Family)))]
data_deseq2<-DESeqDataSetFromMatrix(matrixData,
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="mean")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
family_DESeq2<-cbind("TaxonomicLevel"=rep("Family", nrow(res_final)), res_final)
rownames(family_DESeq2)<-rownames(matrixData)
 
# Fitting models with limma
dummy_matrix<-matrix(0,ncol=length(unique(family_plot_data[,"Family"])), 
	nrow=length(unique(family_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(family_plot_data[counter,"Family"]==unique(family_plot_data[,"Family"])[icol],
				family_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(family_plot_data[,"Family"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(family_plot_data[,"Family"])
rownames(dummy_matrix)<-unique(family_plot_data[,"Sample"])
  
design_data<-family_plot_data[!duplicated(family_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
family_matrix<-t(dummy_matrix)
fit_family<-lmFit(family_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_family<-contrasts.fit(fit_family, cont.matrix)
fit2_family<-eBayes(fit2_family)
pvalues_family<-apply(fit2_family$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
family_limma<-cbind("TaxonomicLevel"=rep("Family",nrow(pvalues_family)),
	 fit2_family$coefficients,pvalues_family)
rownames(family_limma)<-rownames(fit2_family$coefficients)


phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

family_plot_data[,"Abundance"]<-as.numeric(family_plot_data[,"Abundance"])
family_plot_data[,"BMInorm"]<-as.numeric(family_plot_data[,"BMInorm"])
family_plot_data$Sample <- factor(family_plot_data$Sample, levels = family_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/FamilyvsSample.pdf')              
p = ggplot(family_plot_data, aes(x = Sample, y = Abundance, fill = Family, show.legend=FALSE))
p = p + geom_bar(aes(color=Family, fill=Family),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") + 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
print(p)
dev.off()

family_plot_data[,"Abundance"]<-as.numeric(family_plot_data[,"Abundance"])
family_plot_data$Sample <- factor(family_plot_data$Sample, levels = family_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/FamilyvsGroup_facetwrap.pdf')              
p = ggplot(family_plot_data, aes(x = Sample, y = Abundance, fill = Family, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) +
	facet_wrap(~Family, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/FamilyvsObese_facetwrap.pdf')              
p = ggplot(family_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 10)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Family, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()

# For the genus level
#==================================================================================================================
ngenus<-sort(tapply(dd$Abundance,dd$Genus,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
genus_plot_data<-NULL
for (isample in nsample){
	for (igenus in names(ngenus)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_genus<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_igenus<-isample_melt[which(isample_melt$Genus==igenus),]
		count_igenus<-sum(as.numeric(isample_igenus[,"Abundance"]))
		abundance_igenus<-sum(as.numeric(isample_igenus[,"Abundance"]))/total_genus*100
		if (is.null(genus_plot_data)){
			genus_plot_data<-c(isample,igenus,count_igenus,abundance_igenus,
				isample_melt[which(isample_melt$Genus==igenus),"BMI"][1],
				isample_melt[which(isample_melt$Genus==igenus),"Country"][1],
				isample_melt[which(isample_melt$Genus==igenus),"obese"][1],
				isample_melt[which(isample_melt$Genus==igenus),"AGE"][1],
				isample_melt[which(isample_melt$Genus==igenus),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Genus==igenus),"formicacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"aceticacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"propionicacid"][1],							
				isample_melt[which(isample_melt$Genus==igenus),"butyricacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"isovalericacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"adiponectin"][1],			
				isample_melt[which(isample_melt$Genus==igenus),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Genus==igenus),"leptin"][1],
				isample_melt[which(isample_melt$Genus==igenus),"insulin"][1],
				isample_melt[which(isample_melt$Genus==igenus),"glucosemmol"][1],
				isample_melt[which(isample_melt$Genus==igenus),"insulinmmol"][1],
				isample_melt[which(isample_melt$Genus==igenus),"homa_ir"][1],
				isample_melt[which(isample_melt$Genus==igenus),"ratio"][1])
		} else {
			genus_plot_data<-rbind(genus_plot_data,
				c(isample,igenus,count_igenus,abundance_igenus,
				isample_melt[which(isample_melt$Genus==igenus),"BMI"][1],
				isample_melt[which(isample_melt$Genus==igenus),"Country"][1],
				isample_melt[which(isample_melt$Genus==igenus),"obese"][1],
				isample_melt[which(isample_melt$Genus==igenus),"AGE"][1],
				isample_melt[which(isample_melt$Genus==igenus),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Genus==igenus),"formicacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"aceticacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"propionicacid"][1],							
				isample_melt[which(isample_melt$Genus==igenus),"butyricacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"isovalericacid"][1],
				isample_melt[which(isample_melt$Genus==igenus),"adiponectin"][1],			
				isample_melt[which(isample_melt$Genus==igenus),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Genus==igenus),"leptin"][1],
				isample_melt[which(isample_melt$Genus==igenus),"insulin"][1],
				isample_melt[which(isample_melt$Genus==igenus),"glucosemmol"][1],
				isample_melt[which(isample_melt$Genus==igenus),"insulinmmol"][1],
				isample_melt[which(isample_melt$Genus==igenus),"homa_ir"][1],
				isample_melt[which(isample_melt$Genus==igenus),"ratio"][1]))
		}
	}
}
colnames(genus_plot_data)<- c("Sample","Genus","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(genus_plot_data[,"BMI"])-min(as.numeric(genus_plot_data[,"BMI"])))/
	(max(as.numeric(genus_plot_data[,"BMI"]))-min(as.numeric(genus_plot_data[,"BMI"])))*100
genus_plot_data<-data.frame(cbind(genus_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(genus_plot_data)<- c("Sample","Genus","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Genus=genus_plot_data$Genus[1:length(unique(genus_plot_data$Genus))],
	Genus_dummy=genus_plot_data$Genus[1:length(unique(genus_plot_data$Genus))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(genus_plot_data$Counts),nrow=length(unique(genus_plot_data$Genus)),byrow=FALSE)
rownames(matrixData)<-genus_plot_data$Genus[1:length(unique(genus_plot_data$Genus))]
colnames(matrixData)<-genus_plot_data$Sample[seq(1,nrow(genus_plot_data),by=length(unique(genus_plot_data$Genus)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
genus_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
genus_metagenomeSeq_filter = filterData(genus_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(genus_metagenomeSeq_filter )
genus_metagenomeSeq_filter <- cumNorm(genus_metagenomeSeq_filter, p =p) 
genus_metagenomeSeq_filter_nor = MRcounts(genus_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(genus_metagenomeSeq_filter))
pd <- pData(genus_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(genus_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(genus_metagenomeSeq_filter)$obese
AGE = pData(genus_metagenomeSeq_filter)$AGE 
normFactor = normFactors(genus_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = genus_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(genus_metagenomeSeq_filter)$counts))
  
obese = pData(genus_metagenomeSeq_filter)$obese
AGE = pData(genus_metagenomeSeq_filter)$AGE 
Country= pData(genus_metagenomeSeq_filter)$Country
normFactor = normFactors(genus_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = genus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(genus_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(genus_metagenomeSeq_filter)$obese,
	pData(genus_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age")
# fitting the ZIG model 
res = fitZig(obj = genus_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	(obesity1_USA + obesity1_Ghana+obesity0_USA) - obesity0_Ghana,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_genus<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
genus_mG<-cbind("TaxonomicLevel"=rep("Genus",nrow(pvalues_genus)),
	 fit2$coefficients,pvalues_genus)
rownames(genus_mG)<-rownames(fit2$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
matrixData<-matrix(as.numeric(genus_plot_data$Counts),nrow=length(unique(genus_plot_data$Genus)),byrow=FALSE)
rownames(matrixData)<-genus_plot_data$Genus[1:length(unique(genus_plot_data$Genus))]
colnames(matrixData)<-genus_plot_data$Sample[seq(1,nrow(genus_plot_data),by=length(unique(genus_plot_data$Genus)))]
data_deseq2<-DESeqDataSetFromMatrix(matrixData,
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="local")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
genus_DESeq2<-cbind("TaxonomicLevel"=rep("Genus", nrow(res_final)), res_final)
rownames(genus_DESeq2)<-rownames(dds)
 
# Fitting models with limma
dummy_matrix<-matrix(0,ncol=length(unique(genus_plot_data[,"Genus"])), 
	nrow=length(unique(genus_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(genus_plot_data[counter,"Genus"]==unique(genus_plot_data[,"Genus"])[icol],
				genus_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(genus_plot_data[,"Genus"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(genus_plot_data[,"Genus"])
rownames(dummy_matrix)<-unique(genus_plot_data[,"Sample"])
  
design_data<-genus_plot_data[!duplicated(genus_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
genus_matrix<-t(dummy_matrix)
fit_genus<-lmFit(genus_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_genus<-contrasts.fit(fit_genus, cont.matrix)
fit2_genus<-eBayes(fit2_genus)
pvalues_genus<-apply(fit2_genus$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
genus_limma<-cbind("TaxonomicLevel"=rep("Genus",nrow(pvalues_genus)),
	 fit2_genus$coefficients,pvalues_genus)
rownames(genus_limma)<-rownames(fit2_genus$coefficients)


phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

genus_plot_data[,"Abundance"]<-as.numeric(genus_plot_data[,"Abundance"])
genus_plot_data[,"BMInorm"]<-as.numeric(genus_plot_data[,"BMInorm"])
genus_plot_data$Sample <- factor(genus_plot_data$Sample, levels = genus_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/GenusvsSample.pdf')              
p = ggplot(genus_plot_data, aes(x = Sample, y = Abundance, fill = Genus, show.legend=FALSE))
p = p + geom_bar(aes(color=Genus, fill=Genus),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") + 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
print(p)
dev.off()

genus_plot_data[,"Abundance"]<-as.numeric(genus_plot_data[,"Abundance"])
genus_plot_data$Sample <- factor(genus_plot_data$Sample, levels = genus_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/GenusvsGroup_facetwrap.pdf')              
p = ggplot(genus_plot_data, aes(x = Sample, y = Abundance, fill = Genus, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 10)) +
	facet_wrap(~Genus, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/GenusvsObese_facetwrap.pdf')              
p = ggplot(genus_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 10)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Genus, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()


# Species level- CANNOT TRUST THIS RESULTS
#---------------------------------------------------------
nspecies<-sort(tapply(dd$Abundance,dd$Species,sum),decreasing=TRUE)
nsample<-unique(dd$Sample)
species_plot_data<-NULL
for (isample in nsample){
	for (ispecies in names(nspecies)){
		isample_melt<-dd[which(dd$Sample==isample),]
		total_species<-sum(as.numeric(isample_melt[,"Abundance"]))
		isample_ispecies<-isample_melt[which(isample_melt$Species==ispecies),]
		count_ispecies<-sum(as.numeric(isample_ispecies[,"Abundance"]))
		abundance_ispecies<-sum(as.numeric(isample_ispecies[,"Abundance"]))/total_species*100
		if (is.null(species_plot_data)){
			species_plot_data<-c(isample,ispecies,count_ispecies,abundance_ispecies,
				isample_melt[which(isample_melt$Species==ispecies),"BMI"][1],
				isample_melt[which(isample_melt$Species==ispecies),"Country"][1],
				isample_melt[which(isample_melt$Species==ispecies),"obese"][1],
				isample_melt[which(isample_melt$Species==ispecies),"AGE"][1],
				isample_melt[which(isample_melt$Species==ispecies),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Species==ispecies),"formicacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"aceticacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"propionicacid"][1],							
				isample_melt[which(isample_melt$Species==ispecies),"butyricacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"isovalericacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"adiponectin"][1],			
				isample_melt[which(isample_melt$Species==ispecies),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Species==ispecies),"leptin"][1],
				isample_melt[which(isample_melt$Species==ispecies),"insulin"][1],
				isample_melt[which(isample_melt$Species==ispecies),"glucosemmol"][1],
				isample_melt[which(isample_melt$Species==ispecies),"insulinmmol"][1],
				isample_melt[which(isample_melt$Species==ispecies),"homa_ir"][1],
				isample_melt[which(isample_melt$Species==ispecies),"ratio"][1])
		} else {
			species_plot_data<-rbind(species_plot_data,
				c(isample,ispecies,count_ispecies,abundance_ispecies,
				isample_melt[which(isample_melt$Species==ispecies),"BMI"][1],
				isample_melt[which(isample_melt$Species==ispecies),"Country"][1],
				isample_melt[which(isample_melt$Species==ispecies),"obese"][1],
				isample_melt[which(isample_melt$Species==ispecies),"AGE"][1],
				isample_melt[which(isample_melt$Species==ispecies),"Percentbodyfatmass"][1],
				isample_melt[which(isample_melt$Species==ispecies),"formicacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"aceticacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"propionicacid"][1],							
				isample_melt[which(isample_melt$Species==ispecies),"butyricacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"isovalericacid"][1],
				isample_melt[which(isample_melt$Species==ispecies),"adiponectin"][1],			
				isample_melt[which(isample_melt$Species==ispecies),"GlucoseResultGLUCRSLT"][1],
				isample_melt[which(isample_melt$Species==ispecies),"leptin"][1],
				isample_melt[which(isample_melt$Species==ispecies),"insulin"][1],
				isample_melt[which(isample_melt$Species==ispecies),"glucosemmol"][1],
				isample_melt[which(isample_melt$Species==ispecies),"insulinmmol"][1],
				isample_melt[which(isample_melt$Species==ispecies),"homa_ir"][1],
				isample_melt[which(isample_melt$Species==ispecies),"ratio"][1]))
		}
	}
}
colnames(species_plot_data)<- c("Sample","Species","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT", "adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir","ratio")
	
BMInorm=(as.numeric(species_plot_data[,"BMI"])-min(as.numeric(species_plot_data[,"BMI"])))/
	(max(as.numeric(species_plot_data[,"BMI"]))-min(as.numeric(species_plot_data[,"BMI"])))*100
species_plot_data<-data.frame(cbind(species_plot_data,BMInorm), stringsAsFactors=FALSE)
colnames(species_plot_data)<- c("Sample","Species","Counts","Abundance","BMI","Country","obese",
	"AGE","Percentbodyfatmass","formicacid","aceticacid","propionicacid","butyricacid","isovalericacid",
	"GlucoseResultGLUCRSLT","adiponectin","leptin","insulin","glucosemmol","insulinmmol","homa_ir",
	"ratio","BMInorm")

# Fitting linear models with metagenomeSeq
featureData =data.frame(Species=species_plot_data$Species[1:length(unique(species_plot_data$Species))],
	Species_dummy=species_plot_data$Species[1:length(unique(species_plot_data$Species))])
rownames(featureData)<-featureData[,1]
matrixData<-matrix(as.numeric(species_plot_data$Counts),nrow=length(unique(species_plot_data$Species)),byrow=FALSE)
rownames(matrixData)<-species_plot_data$Species[1:length(unique(species_plot_data$Species))]
colnames(matrixData)<-species_plot_data$Sample[seq(1,nrow(species_plot_data),by=length(unique(species_plot_data$Species)))]
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
species_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
species_metagenomeSeq_filter = filterData(species_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(species_metagenomeSeq_filter )
species_metagenomeSeq_filter <- cumNorm(species_metagenomeSeq_filter, p =p) 
species_metagenomeSeq_filter_nor = MRcounts(species_metagenomeSeq_filter , norm = TRUE, log = FALSE)
head(normFactors(species_metagenomeSeq_filter))
pd <- pData(species_metagenomeSeq_filter) 
mod_obese <- model.matrix(~1 + obese, data = pd)
obese_res = fitFeatureModel(species_metagenomeSeq_filter, mod_obese) 
head(MRcoefs(obese_res ))

obese = pData(species_metagenomeSeq_filter)$obese
AGE = pData(species_metagenomeSeq_filter)$AGE 
normFactor = normFactors(species_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age= model.matrix(~obese + AGE + normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age = fitZig(obj = species_metagenomeSeq_filter, mod = mod_obese_age, useCSSoffset =TRUE, control = settings)
results_obese_age<-MRfulltable(fit_obese_age, number = nrow(assayData(species_metagenomeSeq_filter)$counts))
  
obese = pData(species_metagenomeSeq_filter)$obese
AGE = pData(species_metagenomeSeq_filter)$AGE 
Country= pData(species_metagenomeSeq_filter)$Country
normFactor = normFactors(species_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
mod_obese_age_country= model.matrix(~obese + AGE + Country+normFactor) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = species_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(species_metagenomeSeq_filter)$counts))

settings = zigControl(maxit = 100, verbose = TRUE) 
obese_country<-paste("obesity",
	paste(pData(species_metagenomeSeq_filter)$obese,
	pData(species_metagenomeSeq_filter)$Country,sep="_"),
	sep="")
mod = model.matrix(~0+obese_country+AGE) 
colnames(mod)<-c("obesity0_Ghana","obesity0_USA","obesity1_Ghana","obesity1_USA","Age","normFactor")
# fitting the ZIG model 
res = fitZig(obj = species_metagenomeSeq_filter, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = res$fit 
finalMod = res$fit$design 
contrast.matrix = makeContrasts(obesity1_Ghana - obesity0_Ghana, obesity1_USA - obesity1_Ghana, 
	obesity1_USA - obesity0_USA, obesity0_USA - obesity0_Ghana,
	(obesity1_USA + obesity0_USA) - (obesity1_Ghana - obesity0_Ghana), 
	(obesity1_USA + obesity1_Ghana) - (obesity0_USA - obesity0_Ghana),
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues_species<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
species_mG<-cbind("TaxonomicLevel"=rep("Species",nrow(pvalues_species)),
	 fit2$coefficients,pvalues_species)
rownames(species_mG)<-rownames(fit2$coefficients)

# Fitting linear models with DESeq2
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered)))]
colData<-data.frame(metadata, obese_country=obese_country)
Ghana0=matrix(0,ncol=1,nrow=nrow(colData))
Ghana0[which(colData$obese_country=="obesity0_Ghana")]<-1
Ghana1=matrix(0,ncol=1,nrow=nrow(colData))
Ghana1[which(colData$obese_country=="obesity1_Ghana")]<-1
USA0=matrix(0,ncol=1,nrow=nrow(colData))
USA0[which(colData$obese_country=="obesity0_USA")]<-1
USA1=matrix(0,ncol=1,nrow=nrow(colData))
USA1[which(colData$obese_country=="obesity1_USA")]<-1
colData<-data.frame(colData, Ghana0, Ghana1, USA0, USA1)
matrixData<-matrix(as.numeric(species_plot_data$Counts),nrow=length(unique(species_plot_data$Species)),byrow=FALSE)
rownames(matrixData)<-species_plot_data$Species[1:length(unique(species_plot_data$Species))]
colnames(matrixData)<-species_plot_data$Sample[seq(1,nrow(species_plot_data),by=length(unique(species_plot_data$Species)))]
data_deseq2<-DESeqDataSetFromMatrix(matrixData,
	colData=colData, ~AGE+ obese_country, tidy = FALSE, ignoreRank = FALSE)
dds <- DESeq(data_deseq2, fitType="mean")
res1 <- results(dds,contrast=c("obese_country", "obesity1_Ghana","obesity0_Ghana"))
# Use log2FoldChange and padj
res2 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity1_Ghana"))
res3 <- results(dds,contrast=c("obese_country", "obesity1_USA","obesity0_USA"))
res6 <- results(dds,contrast=c("obese_country", "obesity0_USA","obesity0_Ghana"))
res4 <- results(dds,contrast= c(0,0,-1,1,-1,1))
res5 <- results(dds,contrast=c(0,0,-1,-1,1,1))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_coeff=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"])
species_DESeq2<-cbind("TaxonomicLevel"=rep("Species", nrow(res_final)), res_final)
rownames(species_DESeq2)<-rownames(dds)
 
# Fitting models with limma
dummy_matrix<-matrix(0,ncol=length(unique(species_plot_data[,"Species"])), 
	nrow=length(unique(species_plot_data[,"Sample"])))
	counter=1
	for (irow in 1:nrow(dummy_matrix)){
		for (icol in 1: ncol(dummy_matrix)){	
			dummy_matrix[irow, icol]<-ifelse(species_plot_data[counter,"Species"]==unique(species_plot_data[,"Species"])[icol],
				species_plot_data[counter,"Abundance"],0)
			counter=counter+1
		}
	} 
colnames(dummy_matrix)<-unique(species_plot_data[,"Species"])
dummy_matrix<-matrix(as.numeric(dummy_matrix),ncol=ncol(dummy_matrix),
	nrow=nrow(dummy_matrix))
colnames(dummy_matrix)<-unique(species_plot_data[,"Species"])
rownames(dummy_matrix)<-unique(species_plot_data[,"Sample"])
  
design_data<-species_plot_data[!duplicated(species_plot_data[,c("Sample","Country","BMI","obese","AGE")]),
	c("Country","BMI","obese","AGE")]

design<-model.matrix(~design_data[,"obese"]*design_data[,"Country"])
colnames(design)<-c("Ghanalean","obesevsleanGhana","GhanavsUSAlean","interaction")
species_matrix<-t(dummy_matrix)
fit_species<-lmFit(species_matrix,design)
cont.matrix<-cbind(obesity=c(0,2,0,1), country=c(0,0,2,1),
	obesityinGhana=c(0,1,0,0),obesityUSA=c(0,1,0,1),
	USAvsGhanainlean=c(0,1,0,0),
	interaction=c(0,0,0,1))
fit2_species<-contrasts.fit(fit_species, cont.matrix)
fit2_species<-eBayes(fit2_species)
pvalues_species<-apply(fit2_species$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
species_limma<-cbind("TaxonomicLevel"=rep("Species",nrow(pvalues_species)),
	 fit2_species$coefficients,pvalues_species)
rownames(species_limma)<-rownames(fit2_species$coefficients)


phy <- c(Actinobacteria="#ff7f00", Bacteroidetes="#6a3d9a", 
         Cyanobacteria ="#b15928", 
         Firmicutes="#33A02C", Fusobacteria="#fb9a99",
        Lentisphaerae="#fdbf6f", Proteobacteria="#A6CEE3",
         Spirochaetes="#1F78B4", Synergistetes="#B2DF8A", 
         Tenericutes ="#e31a1c",TM7="#4D4D4D",
         Verrucomicrobia="#ffff99") 

species_plot_data[,"Abundance"]<-as.numeric(species_plot_data[,"Abundance"])
species_plot_data[,"BMInorm"]<-as.numeric(species_plot_data[,"BMInorm"])
species_plot_data$Sample <- factor(species_plot_data$Sample, levels = species_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/SpeciesvsSample.pdf')              
p = ggplot(species_plot_data, aes(x = Sample, y = Abundance, fill = Species, show.legend=FALSE))
p = p + geom_bar(aes(color=Species, fill=Species),stat="identity", position="stack",colour="white")+		
	theme_bw() + ylab("Proportions") + 
	guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) 
p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
print(p)
dev.off()

species_plot_data[,"Abundance"]<-as.numeric(species_plot_data[,"Abundance"])
species_plot_data$Sample <- factor(species_plot_data$Sample, levels = species_plot_data$Sample)
pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/SpeciesvsGroup_facetwrap.pdf')              
p = ggplot(species_plot_data, aes(x = Sample, y = Abundance, fill = Species, show.legend=FALSE))
p = p + geom_bar(aes(color=obese, fill=obese),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 12)) +
	facet_wrap(~Species, scale="free") +theme(legend.position="none")
p2=p
print(p2)
dev.off()

pdf(file = '/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/SpeciesvsObese_facetwrap.pdf')              
p = ggplot(species_plot_data, aes(x = BMI, y = Abundance, fill = Country, show.legend=FALSE))
p = p + geom_bar(aes(color=Country, fill=Country),stat="identity", position="stack")+		
	theme_bw() + ylab("Proportions") 
p = p + theme(text = element_text(size = 10)) + theme(axis.text.x=element_blank()) +
	facet_wrap(~Species, scale="free")+theme(legend.position="none")
p2=p
print(p2)
dev.off()
	
	
# Summary of the data
mG_taxa_sig<-rbind(phylum_mG,class_mG,order_mG,family_mG,genus_mG)
DESeq2_taxa_sig<-rbind(phylum_DESeq2,class_DESeq2,order_DESeq2,family_DESeq2,genus_DESeq2)
limma_taxa_sig<-rbind(phylum_limma,class_limma,order_limma,family_limma,genus_limma)
write.table(mG_taxa_sig, file='/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/mG_taxa_sig.txt', 
	sep="\t", row.names=FALSE)
write.table(DESeq2_taxa_sig, file='/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/mG_taxa_sig.txt', 
	sep="\t", row.names=FALSE)
	write.table(limma_taxa_sig, file='/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/mG_taxa_sig.txt', 
	sep="\t", row.names=FALSE)
	
# Plot richeness as a function of different values
#==================================================================================

dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref")

# Normalized data
featureData =data.frame(tax_table(final_filtered))
matrixData<-matrix(otu_table(final_filtered),ncol=ncol(otu_table(final_filtered)))
rownames(matrixData)<-rownames(otu_table(final_filtered))
colnames(matrixData)<-colnames(otu_table(final_filtered))
metadata<-sample_data(final_filtered)[match(colnames(matrixData),rownames(sample_data(final_filtered))), ]
otus_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
	
otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeq_filter )
otus_metagenomeSeq_filter <- cumNorm(otus_metagenomeSeq_filter, p =p) 
otus_metagenomeSeq_filter_nor = MRcounts(otus_metagenomeSeq_filter , norm = TRUE, log = FALSE)
matrixData_filter<-otus_metagenomeSeq_filter_nor
matrixData_filter<-apply(matrixData_filter,c(1,2),function(x){round(x,0)})
metadata_filter<-metadata[ match(colnames(matrixData_filter),rownames(metadata)),]

final_mGnorm<-phyloseq(otu_table(matrixData_filter,taxa_are_rows=TRUE), sample_data(metadata_filter),tax_table(final_filtered), tree)
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
richness=estimate_richness(final_mGnorm, split = TRUE, measures = c("Shannon", "InvSimpson", "Chao1"))
richness=data.frame(sample_data(final_mGnorm), country_obese=country_obese, richness,stringsAsFactors=FALSE)
richness_pvalues<-c()
aa<-summary(lm(Chao1~ AGE+country_obese, data=richness))
richness_pvalues<-cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients)

aa<-summary(lm(Shannon~ AGE+country_obese, data=richness))
richness_pvalues<-rbind(richness_pvalues,
	cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))
	
aa<-summary(lm(InvSimpson~ AGE+country_obese, data=richness))
richness_pvalues<-rbind(richness_pvalues,
	cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+obese, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+obese, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+obese, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+BMI+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+BMI+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+BMI+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+WAIST+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+WAIST+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+WAIST+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+Percentbodyfatmass+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+Percentbodyfatmass+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+Percentbodyfatmass+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+meansfa+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+meansfa+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+meansfa+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+aceticacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+aceticacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+aceticacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+propionicacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+propionicacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+propionicacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+butyricacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+butyricacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+butyricacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+isovalericacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+isovalericacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+isovalericacid+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+adiponectin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+adiponectin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~AGE+adiponectin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+leptin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+leptin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~AGE+leptin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+insulin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+insulin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~AGE+insulin+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+glucosemmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+glucosemmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~AGE+glucosemmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~AGE+ insulinmmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+insulinmmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~AGE+insulinmmol+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Chao1~ AGE+homa_ir+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Chao1",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(Shannon~ AGE+homa_ir+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("Shannon",nrow(aa$coefficients)),aa$coefficients))

# aa<-summary(lm(InvSimpson~ AGE+homa_ir+obese+Country, data=richness))
# richness_pvalues<-rbind(richness_pvalues,
	# cbind(Method=rep("InvSimpson",nrow(aa$coefficients)),aa$coefficients))

write.table(richness_pvalues, 
	file='/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/richness_pvalues.txt', 
	sep="\t", row.names=FALSE)

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obese_country_Chao1.pdf")
p = ggplot(data=richness, aes(x = country_obese, y = Chao1, fill = AGE, show.legend=FALSE))
p = p + geom_point(aes(color=AGE, fill=AGE),stat="identity")+		
	theme_bw() + ylab("Chao1") 
 p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))+geom_boxplot(colour=c("red","blue","purple","darkgreen"),alpha=0.05,fill=c("red","blue","purple","darkgreen"))
 print(p)
 dev.off()
 
 pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obese_country_Shannon.pdf")
p = ggplot(data=richness, aes(x = country_obese, y = Shannon, fill = AGE, show.legend=FALSE))
p = p + geom_point(aes(color=AGE, fill=AGE),stat="identity")+		
	theme_bw() + ylab("Shannon") 
 p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))+geom_boxplot(colour=c("red","blue","purple","darkgreen"),alpha=0.05,fill=c("red","blue","purple","darkgreen"))
 print(p)
 dev.off()

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obese_country_InvSimpson.pdf")
p = ggplot(data=richness, aes(x = country_obese, y = InvSimpson, fill = AGE, show.legend=FALSE))
p = p + geom_point(aes(color=AGE, fill=AGE),stat="identity")+		
	theme_bw() + ylab("InvSimpson") 
 p = p + theme(text = element_text(size = 12)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))+geom_boxplot(colour=c("red","blue","purple","darkgreen"),alpha=0.05,fill=c("red","blue","purple","darkgreen"))
  print(p)
 dev.off()


# dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref")
# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/Countrybyobese.pdf")
# p1 = plot_richness(final_filtered, x="Country", color="obese", measures = c("Shannon", "InvSimpson", "Chao1"))
# p1 + theme(text = element_text(size = 8)) +
     # theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obeseyscfa_tot.pdf")
# p1 = plot_richness(final_filtered, x="Country", color="Percentbodyfatmass", measures = c("Shannon", "InvSimpson", "Chao1"))
# p2=p1 + theme(text = element_text(size = 8)) + scale_colour_gradient(low="darkblue", high="red")+
     # theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
    # print(p2)
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obesebyaceticac.pdf")
# p1 = plot_richness(final_filtered, x="obese", color="aceticacid", shape="Country", measures = c("Shannon", "InvSimpson", "Chao1"))
# p1 + theme(text = element_text(size = 8)) + scale_colour_gradient(low="darkblue", high="red")+
     # theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/RichnessPlots_openref/obesebybutyricac.pdf")
# p1 = plot_richness(final_filtered, x="obese", color="butyricacid", shape="Country", measures = c("Shannon", "InvSimpson", "Chao1"))
# p1 + theme(text = element_text(size = 8)) + scale_colour_gradient(low="darkblue", high="red")+
     # theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
# dev.off()

# Calculate Unifrac distances  
#===================================================================================    

#Rarefaction
#===================================================================================
S <- specnumber(t(otu_table(final_filtered))) # observed number of species
raremax <- min(colSums(otu_table(final_filtered)))
Srare <- rarefy(t(otu_table(final_filtered)), raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(t(otu_table(final_filtered)), step = 20, sample = raremax, col = "blue", cex = 0.6,label = TRUE)
final_even500= rarefy_even_depth(final, replace=TRUE,trimOTUs=TRUE,sample.size=500,rngseed=711)

unifrac_weighted_normalized=UniFrac(final_even500, weighted=TRUE, normalized=TRUE, parallel=FALSE, 
	fast=TRUE)
unifrac_weighted_nonnormalized=UniFrac(final_even500, weighted=TRUE, normalized=FALSE, parallel=FALSE, 
	fast=TRUE)
unifrac_unweighted_nonnormalized=UniFrac(final_even500, weighted=FALSE, normalized=FALSE, parallel=FALSE, 
	fast=TRUE)           
unifrac_unweighted_normalized=UniFrac(final_even500, weighted=FALSE, normalized=TRUE, parallel=FALSE, 
	fast=TRUE)  

#PCoA of beta-diversity (unweighted unifrac)
#===================================================================================
dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/")
aaa<-data.frame(Country=sample_data(final_even500)[,"Country"],obese=sample_data(final_even500)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")

metadata<-sample_data(data.frame(sample_data(final_even500), country_obese=aaa))
final_even500=phyloseq(otu_table(final_even500),metadata,tax_table(final_even500),tree)

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_Unweigthed_NonNormalized.pdf")
unweighted_u.pco2 = ordinate(final_even500, method = "PCoA", 
	distance = unifrac_unweighted_nonnormalized)
p2.2 = plot_ordination(final_even500, unweighted_u.pco2, type = "samples", 
	color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + 
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
dev.off()
 
pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_Unweigthed_Normalized.pdf")     
unweighted_n.pco2 = ordinate(final_even500, method = "PCoA", distance = unifrac_unweighted_normalized)
p2.2 = plot_ordination(final_even500, unweighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + 
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
dev.off()
        
pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_Weigthed_NonNormalized.pdf")
weighted_u.pco2 = ordinate(final_even500, method = "PCoA", distance = unifrac_weighted_nonnormalized)
p2.2 = plot_ordination(final_even500, weighted_u.pco2, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +          guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
dev.off()
 
pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_Weigthed_Normalized.pdf")     
weighted_n.pco2 = ordinate(final_even500, method = "PCoA", distance = unifrac_weighted_normalized)
p2.2 = plot_ordination(final_even500, weighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10)))   
dev.off()
         
#NMDS of beta-diversity (unweighted unifrac)
#===================================================================================
pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_NMDS_Unweigthed_NonNormalized.pdf")     
unweighted_u.nmds = ordinate(final_even500, method = "NMDS", distance = unifrac_unweighted_nonnormalized, try=50,trymax=100)
p2.2 = plot_ordination(final_even500, unweighted_u.nmds, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10)))   
dev.off()

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_NMDS_Unweigthed_Normalized.pdf")     
unweighted_n.nmds = ordinate(final_even500, method = "NMDS", distance = unifrac_unweighted_normalized)
p2.2 = plot_ordination(final_even500, unweighted_n.nmds, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10)))   
dev.off()

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_NMDS_Weigthed_NonNormalized.pdf")     
weighted_u.nmds = ordinate(final_even500, method = "NMDS", distance = unifrac_weighted_nonnormalized)
p2.2 = plot_ordination(final_even500, weighted_u.nmds, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10)))   
dev.off()

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/UniFrac_NMDS_Weigthed_Normalized.pdf")     
weighted_n.nmds = ordinate(final_even500, method = "NMDS", distance = unifrac_unweighted_normalized)
p2.2 = plot_ordination(final_even500, weighted_n.nmds, type = "samples", color = "obese", shape = "Country") 
p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +          guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10)))   
dev.off()
        

#PCoA of beta-diversity (bray-curtis)
#===================================================================================
pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA.pdf")
bray.pcoa = ordinate(final_mGnorm, method = "PCoA", distance = "bray")
p3.1 = plot_ordination(final_mGnorm, bray.pcoa, type = "sample", color = "obese",shape = "Country")
p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")+
        guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
dev.off()

pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA.pdf")
bray.pcoa = ordinate(final_mGnorm, method = "NMDS", distance = "bray")
p3.1 = plot_ordination(final_mGnorm, bray.pcoa, type = "sample", color = "obese",shape = "Country")
p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")+ guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
     
dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_DCA.pdf")
# bray.DCA = ordinate(final_even500, method = "DCA", distance = "bray")
# p3.1 = plot_ordination(final_even500, bray.DCA, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_CCA.pdf")
# bray.CCA = ordinate(final_even500, method = "CCA", distance = "bray")
# p3.1 = plot_ordination(final_even500, bray.CCA, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_RDA.pdf")
# bray.RDA = ordinate(final_filtered, method = "RDA", distance = "bray")
# p3.1 = plot_ordination(final_filtered, bray.RDA, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_CAP.pdf")
# bray.CAP = ordinate(final_even500, method = "CAP", distance = "bray")
# p3.1 = plot_ordination(final_even500, bray.CAP, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_DPCoA.pdf")
# bray.DPCoA = ordinate(final_even500, method = "DPCoA", distance = "bray")
# p3.1 = plot_ordination(final_even500, bray.DPCoA, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_MDS.pdf")
# bray.MDS = ordinate(final_even500, method = "MDS", distance = "bray")
# p3.1 = plot_ordination(final_even500, bray.MDS, type = "sample", color = "obese",shape = "Country")
# p3.1 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) +  scale_colour_gradient(low="red", high="blue")
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

 
# #PCoA+NMDS on top 10 phyla across samples using bray distance
# #==================================================================================

# otu_table_norm<-otu_table(final_filtered)
 # otu_table_norm=otu_table_norm/colSums(otu_table_norm)
 # final_norm<-merge_phyloseq( otu_table_norm, sample_data_order_phy,tax_split,tree)
# class.sum = tapply(taxa_sums( final_norm), tax_table( final_norm)[, "Phylum"], sum, na.rm = TRUE)
# top10phyla = names(sort(class.sum, TRUE))[1:10]
# final_top10phyla= prune_taxa((tax_table( final_norm)[, "Phylum"] %in% top10phyla),  final_norm)

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_NMDS_ByPhylum.pdf")
# final_top10phyla.ord.NMDS = ordinate(final_top10phyla, "NMDS", distance = "bray")
# p4 = plot_ordination(final_top10phyla, final_top10phyla.ord.NMDS, type = "taxa", color = "Phylum")
# p5=p4 + facet_wrap(~Phylum, 5) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA_ByPhylum.pdf")
# final_top10phyla.ord.PCoA= ordinate(final_top10phyla, "PCoA", distance = "bray")
# p4 = plot_ordination(final_top10phyla, final_top10phyla.ord.PCoA, type = "taxa", color = "Phylum")
# p4 + facet_wrap(~Phylum, 5) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_MDS_ByPhylum.pdf")
# final_top10phyla.ord.MDS= ordinate(final_top10phyla, "MDS", distance = "bray")
# p4 = plot_ordination(final_top10phyla, final_top10phyla.ord.MDS, type = "taxa", color = "Phylum")
# p4 + facet_wrap(~Phylum, 5) 
# dev.off()

# #PCoA+NMDS on top 10 genus across samples using bray-curtis distance
# #===================================================================================

# genus.sum = tapply(taxa_sums( final_norm), tax_table( final_norm)[, "Genus"], sum, na.rm = TRUE)
# top10genus = names(sort(genus.sum, TRUE))[1:10]
# final_top10genus= prune_taxa((tax_table( final_norm)[, "Genus"] %in% top10genus),  final_norm)

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_NMDS_ByGenus.pdf")
# final_top10genus.ord.NMDS = ordinate(final_top10genus, "NMDS", distance = "bray")
# p4 = plot_ordination(final_top10genus, final_top10genus.ord.NMDS, type = "taxa", color = "Genus")
# p4 + facet_wrap(~Genus, 5)
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA_ByGenus.pdf")
# final_top10genus.ord.PCoA = ordinate(final_top10genus, "PCoA", distance = "bray")
# p4 = plot_ordination(final_top10genus, final_top10genus.ord.PCoA, type = "taxa", color = "Genus")
# p4 + facet_wrap(~Genus, 5)
# dev.off()

# #Split Graphic
# #===================================================================================
# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA_TopPhylum_Country.pdf")
# p5 = plot_ordination(final_top10phyla,final_top10phyla.ord.PCoA, type = "split", color = "Phylum", shape = "Country")
# p5+scale_shape_manual(values=seq(0,15))
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA_TopPhylum_obese.pdf")
# p5 = plot_ordination(final_top10phyla,final_top10phyla.ord.PCoA, type = "split", color = "obese", 
	# shape = "Phylum")
# p5+scale_shape_manual(values=seq(0,15))
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity_openref/Bray_PCoA_TopGenus_Country.pdf")
# p5 = plot_ordination(final_top10genus,final_top10genus.ord.PCoA, type = "split", color = "Genus", shape = "Country")
# p5+scale_shape_manual(values=seq(0,15))
# dev.off()

# #Bar graph of OTUs by env_type
# #===================================================================================

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/Ruminococcaceae_Country_Obese_dis.pdf")
# final.Bacteroidaceae= subset_taxa(final_filtered, Family == "Bacteroidaceae")
# p6 = plot_bar(final.Bacteroidaceae, x="BMI", y="Abundance", fill="Genus", facet_grid=Country~obese)
# p6 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + 
	# theme(text = element_text(size = 14)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
# dev.off()
        
# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/Ruminococcaceae_Country_Obese_dis.pdf")
# final.Ruminococcaceae= subset_taxa(final_filtered, Family == "Ruminococcaceae")
# p6 = plot_bar(final.Ruminococcaceae, x="BMI", y="Abundance", fill="Genus", facet_grid=Country~obese)
# p6 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + 
	# theme(text = element_text(size = 14)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/Prevotellaceae_Country_Obese_dis.pdf")
# final.Prevotellaceae= subset_taxa(final_filtered, Family == "Prevotellaceae")
# p6 = plot_bar(final.Prevotellaceae, x="BMI", y="Abundance", fill="Genus", facet_grid=Country~obese)
# p6 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + 
	# theme(text = element_text(size = 14)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/AbundancePlots_openref/Clostridiaceae_Country_Obese_dis.pdf")
# final.Clostridiaceae= subset_taxa(final_filtered, Family == "Clostridiaceae")
# p6 = plot_bar(final.Clostridiaceae, x="BMI", y="Abundance", fill="Genus", facet_grid=Country~obese)
# p6 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + 
	# theme(text = element_text(size = 14)) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
# dev.off()


# #Network analysis
# #===================================================================================

# set.seed(711L)
# plot_net(final_normalized, maxdist=0.4, point_label = "Country", color ="obese", shape="Country")
# # set.seed(711L)
# # p5=plot_net(final_normalized, maxdist=0.4, point_label = "CaseName", color ="PMI_dis", 
	# # shape = "Organ")
# # p5+scale_colour_gradient(low="red", high="blue")

# # ig = make_network(final_normalized, "samples", max.dist=0.4, dist.fun="bray")
# # plot_network(ig,final_normalized, point_label = "CaseName",color = "PMI_dis", shape = "Organ", line_weight = 0.3)

# # ig = make_network(final_normalized, "samples", max.dist=0.4, dist.fun="bray")
# # plot_network(ig,final_normalized, color = "AmbientTemperatureCelsius_dis", shape = "Organ", line_weight = 0.3)

# # Specific generalized linear models for certain genus and phyla
# #===================================================================================

# # # Determine the factors the mostly influenced the presence of certain species
# # otu_table_fuso<-otu_table(final_normalized)
# # Fusobacterium<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Genus"]=="Fusobacterium"),]
# # otu_table_fuso<-otu_table_fuso[na.omit(match(rownames(Fusobacterium),rownames(otu_table_fuso))),]
# # data_glm<-cbind(sample_data(final_normalized),Fusobacterium=colSums(otu_table_fuso))

 # # summary(glm(Fusobacterium~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age, species=quasibinomial,
	# # data=data_glm))

# # # Call:
# # # glm(formula = Fusobacterium ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.27787  -0.06635  -0.03090   0.02217   0.19608  

# # # Coefficients:
                                        # # # Estimate Std. Error t value Pr(>|t|)    
# # # (Intercept)                            -2.022795   2.088514  -0.969 0.339067    
# # # CauseofDeathCar Accident               -3.516718   3.178778  -1.106 0.275733    
# # # CauseofDeathChest Trauma               -7.050699   5.192889  -1.358 0.182763    
# # # CauseofDeathComplication after surgery -5.299496   7.052911  -0.751 0.457170    
# # # CauseofDeathCoronary Heart Disease     -2.695781   0.687697  -3.920 0.000369 ***
# # # CauseofDeathDrowning                   -3.552176   1.446247  -2.456 0.018860 *  
# # # CauseofDeathDrug Overdose              -0.482946   0.895439  -0.539 0.592884    
# # # CauseofDeathGSW Chest                  -1.675896   1.137144  -1.474 0.149000    
# # # CauseofDeathGSW Head                   -5.789909   1.722766  -3.361 0.001815 ** 
# # # CauseofDeathMultiple gunshot wounds    -3.038910   1.189097  -2.556 0.014841 *  
# # # CauseofDeathOverdose                   -3.923648   0.798140  -4.916 1.83e-05 ***
# # # CauseofDeathPending Further Studies    -4.298038   1.086535  -3.956 0.000332 ***
# # # CauseofDeathShotgun Wound of Head      -6.465149   2.146223  -3.012 0.004656 ** 
# # # OrganBrain                             -0.027794   0.583323  -0.048 0.962253    
# # # OrganHeart                             -0.177113   0.603116  -0.294 0.770656    
# # # OrganLiver                             -0.835040   0.705194  -1.184 0.243914    
# # # OrganMouthSwab                          2.024678   0.793594   2.551 0.015000 *  
# # # OrganSpleen                            -0.872479   0.693142  -1.259 0.216012    
# # # AmbientTemperatureCelsius              -0.047924   0.068607  -0.699 0.489213    
# # # PMI                                     0.001389   0.013368   0.104 0.917830    
# # # GenderM                                 1.168422   0.494477   2.363 0.023498 *  
# # # Age                                     0.003009   0.026175   0.115 0.909112    
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.01239359)

    # # # Null deviance: 2.21607  on 58  degrees of freedom
# # # Residual deviance: 0.48113  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 10


# # otu_table_prevotella<-otu_table(final_normalized)
# # Prevotella<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Genus"]=="Prevotella"),]
# # otu_table_prevotella<-otu_table_prevotella[na.omit(match(rownames(Prevotella),rownames(otu_table_prevotella))),]
# # data_glm<-cbind(sample_data(final_normalized),Prevotella=colSums(otu_table_prevotella))

 # # summary(glm(Prevotella~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
# # # Call:
# # # glm(formula = Prevotella ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.63111  -0.10544  -0.03623   0.02860   0.39797  

# # # Coefficients:
                                         # # # Estimate Std. Error t value Pr(>|t|)   
# # # (Intercept)                            -7.070e+00  4.074e+00  -1.735  0.09104 . 
# # # CauseofDeathCar Accident               -1.158e+01  2.547e+03  -0.005  0.99640   
# # # CauseofDeathChest Trauma               -2.554e+00  6.095e+00  -0.419  0.67761   
# # # CauseofDeathComplication after surgery -7.100e-01  1.119e+01  -0.063  0.94976   
# # # CauseofDeathCoronary Heart Disease      2.148e+00  3.605e+00   0.596  0.55493   
# # # CauseofDeathDrowning                    2.736e+00  3.765e+00   0.727  0.47199   
# # # CauseofDeathDrug Overdose               5.047e+00  3.743e+00   1.348  0.18570   
# # # CauseofDeathGSW Chest                   3.495e+00  3.711e+00   0.942  0.35234   
# # # CauseofDeathGSW Head                    1.590e+00  3.760e+00   0.423  0.67481   
# # # CauseofDeathMultiple gunshot wounds     3.950e+00  3.690e+00   1.071  0.29130   
# # # CauseofDeathOverdose                    2.578e+00  3.580e+00   0.720  0.47604   
# # # CauseofDeathPending Further Studies     2.304e+00  3.649e+00   0.631  0.53169   
# # # CauseofDeathShotgun Wound of Head      -6.969e-01  4.476e+00  -0.156  0.87711   
# # # OrganBrain                              2.345e-01  1.075e+00   0.218  0.82856   
# # # OrganHeart                             -1.812e-01  1.166e+00  -0.155  0.87738   
# # # OrganLiver                             -1.205e+00  1.531e+00  -0.787  0.43630   
# # # OrganMouthSwab                          3.026e+00  1.065e+00   2.841  0.00727 **
# # # OrganSpleen                            -7.607e-01  1.214e+00  -0.626  0.53485   
# # # AmbientTemperatureCelsius               1.077e-03  6.715e-02   0.016  0.98729   
# # # PMI                                     1.824e-03  1.472e-02   0.124  0.90205   
# # # GenderM                                -4.618e-01  4.479e-01  -1.031  0.30922   
# # # Age                                     6.474e-03  2.383e-02   0.272  0.78740   
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.05610891)

    # # # Null deviance: 8.3341  on 58  degrees of freedom
# # # Residual deviance: 2.1538  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 18

# # otu_table_Bacteroidetes<-otu_table(final_normalized)
# # Bacteroidetes<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Bacteroidetes"),]
# # otu_table_Bacteroidetes<-otu_table_Bacteroidetes[na.omit(match(rownames(Bacteroidetes),
	# # rownames(otu_table_Bacteroidetes))),]
# # data_glm<-cbind(sample_data(final_normalized),Bacteroidetes=colSums(otu_table_Bacteroidetes))

 # # summary(glm(Bacteroidetes~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))

# # # Call:
# # # glm(formula = Bacteroidetes ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.63791  -0.14073  -0.05347   0.05329   0.63381  

# # # Coefficients:
                                        # # # Estimate Std. Error t value Pr(>|t|)    
# # # (Intercept)                             1.263664   2.093437   0.604 0.549772    
# # # CauseofDeathCar Accident               -6.976980   7.509875  -0.929 0.358890    
# # # CauseofDeathChest Trauma               -7.163903   5.626183  -1.273 0.210848    
# # # CauseofDeathComplication after surgery -6.640151  13.750421  -0.483 0.632008    
# # # CauseofDeathCoronary Heart Disease     -2.994874   0.728200  -4.113 0.000209 ***
# # # CauseofDeathDrowning                   -2.251440   1.321718  -1.703 0.096877 .  
# # # CauseofDeathDrug Overdose              -0.208048   1.107237  -0.188 0.851983    
# # # CauseofDeathGSW Chest                  -2.073624   1.234138  -1.680 0.101340    
# # # CauseofDeathGSW Head                   -5.126847   1.509636  -3.396 0.001646 ** 
# # # CauseofDeathMultiple gunshot wounds    -1.757232   1.080082  -1.627 0.112237    
# # # CauseofDeathOverdose                   -3.471189   0.685333  -5.065 1.15e-05 ***
# # # CauseofDeathPending Further Studies    -4.173962   1.132856  -3.684 0.000729 ***
# # # CauseofDeathShotgun Wound of Head      -6.935960   3.427733  -2.023 0.050292 .  
# # # OrganBrain                              0.482460   0.674379   0.715 0.478843    
# # # OrganHeart                             -0.234892   0.733639  -0.320 0.750638    
# # # OrganLiver                              0.061309   0.707323   0.087 0.931395    
# # # OrganMouthSwab                          2.323359   0.802840   2.894 0.006344 ** 
# # # OrganSpleen                            -1.085117   0.810055  -1.340 0.188559    
# # # AmbientTemperatureCelsius              -0.054759   0.068205  -0.803 0.427179    
# # # PMI                                     0.002979   0.015974   0.187 0.853062    
# # # GenderM                                 0.237724   0.485607   0.490 0.627350    
# # # Age                                    -0.018971   0.026743  -0.709 0.482536    
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.09492104)

    # # # Null deviance: 12.2295  on 58  degrees of freedom
# # # Residual deviance:  3.3462  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 9

# # otu_table_Firmicutes<-otu_table(final_normalized)
# # Firmicutes<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Firmicutes"),]
# # otu_table_Firmicutes<-otu_table_Firmicutes[na.omit(match(rownames(Firmicutes),
	# # rownames(otu_table_Firmicutes))),]
# # data_glm<-cbind(sample_data(final_normalized),Firmicutes=colSums(otu_table_Firmicutes))

 # # summary(glm(Firmicutes~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
 # # # Call:
# # # glm(formula = Firmicutes ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -1.29959  -0.35687   0.03066   0.49352   1.24040  

# # # Coefficients:
                                       # # # Estimate Std. Error t value Pr(>|t|)  
# # # (Intercept)                             2.39667    2.75664   0.869   0.3902  
# # # CauseofDeathCar Accident               -4.54595    5.16438  -0.880   0.3844  
# # # CauseofDeathChest Trauma                4.26367    3.45894   1.233   0.2255  
# # # CauseofDeathComplication after surgery  8.74925   25.35857   0.345   0.7320  
# # # CauseofDeathCoronary Heart Disease      1.85291    1.25893   1.472   0.1495  
# # # CauseofDeathDrowning                    2.03255    2.03597   0.998   0.3246  
# # # CauseofDeathDrug Overdose              -1.14170    2.01537  -0.566   0.5745  
# # # CauseofDeathGSW Chest                   0.19214    2.27163   0.085   0.9331  
# # # CauseofDeathGSW Head                   -2.49470    2.51055  -0.994   0.3268  
# # # CauseofDeathMultiple gunshot wounds     0.21653    1.43940   0.150   0.8812  
# # # CauseofDeathOverdose                    1.35808    0.98453   1.379   0.1760  
# # # CauseofDeathPending Further Studies    -0.76228    1.92635  -0.396   0.6946  
# # # CauseofDeathShotgun Wound of Head       2.11820    3.26879   0.648   0.5210  
# # # OrganBrain                             -0.76534    1.13321  -0.675   0.5036  
# # # OrganHeart                             -0.42321    1.13586  -0.373   0.7116  
# # # OrganLiver                             -0.47976    1.13471  -0.423   0.6749  
# # # OrganMouthSwab                         -0.39949    1.27782  -0.313   0.7563  
# # # OrganSpleen                            -0.26383    1.10860  -0.238   0.8132  
# # # AmbientTemperatureCelsius               0.01643    0.06018   0.273   0.7864  
# # # PMI                                     0.02835    0.02521   1.124   0.2681  
# # # GenderM                                 1.79140    0.80790   2.217   0.0328 *
# # # Age                                    -0.09153    0.04937  -1.854   0.0717 .
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.6414088)

    # # # Null deviance: 38.735  on 58  degrees of freedom
# # # Residual deviance: 22.063  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 9

# # otu_table_Actinobacteria<-otu_table(final_normalized)
# # Actinobacteria<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Actinobacteria"),]
# # otu_table_Actinobacteria<-otu_table_Actinobacteria[na.omit(match(rownames(Actinobacteria),
	# # rownames(otu_table_Actinobacteria))),]
# # data_glm<-cbind(sample_data(final_normalized),Actinobacteria=colSums(otu_table_Actinobacteria))

 # # summary(glm(Actinobacteria~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
 	
# # # Call:
# # # glm(formula = Actinobacteria ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.15607  -0.05525  -0.01629   0.01021   0.22910  

# # # Coefficients:
                                         # # # Estimate Std. Error t value Pr(>|t|)   
# # # (Intercept)                              -6.68904    2.20879  -3.028  0.00429 **
# # # CauseofDeathCar Accident                 -2.57081    3.83831  -0.670  0.50685   
# # # CauseofDeathChest Trauma                 -2.04617    1.74395  -1.173  0.24761   
# # # CauseofDeathComplication after surgery  -15.76642 2232.04196  -0.007  0.99440   
# # # CauseofDeathCoronary Heart Disease       -1.26160    0.93378  -1.351  0.18427   
# # # CauseofDeathDrowning                     -0.60216    1.40022  -0.430  0.66947   
# # # CauseofDeathDrug Overdose                 2.59047    1.43344   1.807  0.07826 . 
# # # CauseofDeathGSW Chest                    -1.16482    1.96149  -0.594  0.55596   
# # # CauseofDeathGSW Head                      4.45157    1.28073   3.476  0.00124 **
# # # CauseofDeathMultiple gunshot wounds       0.78258    1.26156   0.620  0.53856   
# # # CauseofDeathOverdose                     -1.83950    0.98944  -1.859  0.07038 . 
# # # CauseofDeathPending Further Studies       0.87766    1.09319   0.803  0.42681   
# # # CauseofDeathShotgun Wound of Head         0.06517    1.20343   0.054  0.95708   
# # # CauseofDeathStroke                       -7.73044    4.27452  -1.808  0.07805 . 
# # # OrganBrain                               -0.46968    1.17057  -0.401  0.69038   
# # # OrganHeart                               -0.25743    1.18755  -0.217  0.82949   
# # # OrganLiver                                0.21628    1.09528   0.197  0.84447   
# # # OrganMouthSwab                            1.83403    1.10593   1.658  0.10507   
# # # OrganSpleen                               0.16271    1.07663   0.151  0.88063   
# # # AmbientTemperatureCelsius                -0.02384    0.06302  -0.378  0.70724   
# # # PMI                                       0.02064    0.01291   1.599  0.11776   
# # # GenderM                                  -1.10341    0.53125  -2.077  0.04427 * 
# # # Age                                       0.04849    0.02979   1.628  0.11143   
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.0158478)

    # # # Null deviance: 5.23106  on 62  degrees of freedom
# # # Residual deviance: 0.39465  on 40  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 19

# # otu_table_Acidobacteria<-otu_table(final_normalized)
# # Acidobacteria<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Acidobacteria"),]
# # otu_table_Acidobacteria<-otu_table_Acidobacteria[na.omit(match(rownames(Acidobacteria),
	# # rownames(otu_table_Acidobacteria))),]
# # data_glm<-cbind(sample_data(final_normalized),Acidobacteria=colSums(otu_table_Acidobacteria))

 # # summary(glm(Acidobacteria~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
# # # Call:
# # # glm(formula = Acidobacteria ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
      # # # Min         1Q     Median         3Q        Max  
# # # -0.061196  -0.006460  -0.000009   0.000000   0.097200  

# # # Coefficients:
                                         # # # Estimate Std. Error t value Pr(>|t|)   
# # # (Intercept)                             1.481e+01  8.507e+00   1.742  0.08990 . 
# # # CauseofDeathCar Accident               -4.344e+01  1.411e+04  -0.003  0.99756   
# # # CauseofDeathChest Trauma                1.239e+01  3.860e+03   0.003  0.99746   
# # # CauseofDeathComplication after surgery -1.095e+01  1.411e+04  -0.001  0.99939   
# # # CauseofDeathCoronary Heart Disease      3.969e-02  3.537e+00   0.011  0.99111   
# # # CauseofDeathDrowning                    8.874e-01  5.368e+00   0.165  0.86959   
# # # CauseofDeathDrug Overdose               5.065e+00  4.546e+00   1.114  0.27240   
# # # CauseofDeathGSW Chest                   5.047e+00  1.463e+04   0.000  0.99973   
# # # CauseofDeathGSW Head                   -1.532e+01  1.463e+04  -0.001  0.99917   
# # # CauseofDeathMultiple gunshot wounds     2.488e+00  4.826e+00   0.516  0.60923   
# # # CauseofDeathOverdose                   -6.258e+00  2.769e+00  -2.260  0.02978 * 
# # # CauseofDeathPending Further Studies    -2.667e+01  7.715e+03  -0.003  0.99726   
# # # CauseofDeathShotgun Wound of Head      -1.463e+01  1.463e+04  -0.001  0.99921   
# # # OrganBrain                             -5.049e-01  2.326e+00  -0.217  0.82933   
# # # OrganHeart                             -2.566e+00  2.698e+00  -0.951  0.34786   
# # # OrganLiver                              1.342e+00  2.276e+00   0.589  0.55917   
# # # OrganMouthSwab                         -1.838e+01  3.860e+03  -0.005  0.99623   
# # # OrganSpleen                            -2.245e+00  2.584e+00  -0.869  0.39049   
# # # AmbientTemperatureCelsius              -8.138e-01  2.598e-01  -3.133  0.00338 **
# # # PMI                                     2.600e-01  9.145e-02   2.844  0.00722 **
# # # GenderM                                 5.593e+00  2.816e+00   1.986  0.05444 . 
# # # Age                                    -2.798e-01  1.105e-01  -2.532  0.01573 * 
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.004268548)

    # # # Null deviance: 0.730517  on 58  degrees of freedom
# # # Residual deviance: 0.039439  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 24


# # # With just the significant variables become irrelevant
# # otu_table_Acidobacteria<-otu_table(final_normalized)
# # Acidobacteria<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Acidobacteria"),]
# # otu_table_Acidobacteria<-otu_table_Acidobacteria[na.omit(match(rownames(Acidobacteria),
	# # rownames(otu_table_Acidobacteria))),]
# # data_glm<-cbind(sample_data(final_normalized),Acidobacteria=colSums(otu_table_Acidobacteria))


# # otu_table_Fusobacteria<-otu_table(final_normalized)
# # Fusobacteria<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Fusobacteria"),]
# # otu_table_Fusobacteria<-otu_table_Fusobacteria[na.omit(match(rownames(Fusobacteria),
	# # rownames(otu_table_Fusobacteria))),]
# # data_glm<-cbind(sample_data(final_normalized),Fusobacteria=colSums(otu_table_Fusobacteria))

 # # summary(glm(Fusobacteria~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
# # # Call:
# # # glm(formula = Fusobacteria ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.26621  -0.08437  -0.02912   0.02152   0.19564  

# # # Coefficients:
                                         # # # Estimate Std. Error t value Pr(>|t|)    
# # # (Intercept)                             0.2777146  2.2557668   0.123 0.902684    
# # # CauseofDeathCar Accident               -4.6301648  3.7094533  -1.248 0.219797    
# # # CauseofDeathChest Trauma               -4.9982512  4.2922590  -1.164 0.251683    
# # # CauseofDeathComplication after surgery -4.9363041  8.2009187  -0.602 0.550898    
# # # CauseofDeathCoronary Heart Disease     -1.6177283  0.5937486  -2.725 0.009772 ** 
# # # CauseofDeathDrowning                   -1.0443277  1.1515031  -0.907 0.370315    
# # # CauseofDeathDrug Overdose               0.4108472  0.9330066   0.440 0.662247    
# # # CauseofDeathGSW Chest                  -0.5023840  1.0433915  -0.481 0.633003    
# # # CauseofDeathGSW Head                   -6.0046637  1.9733379  -3.043 0.004295 ** 
# # # CauseofDeathMultiple gunshot wounds    -1.8149446  1.1935047  -1.521 0.136840    
# # # CauseofDeathOverdose                   -3.5454906  0.8203351  -4.322 0.000112 ***
# # # CauseofDeathPending Further Studies    -4.2844668  1.2322500  -3.477 0.001313 ** 
# # # CauseofDeathShotgun Wound of Head      -5.7448224  2.3672409  -2.427 0.020221 *  
# # # OrganBrain                              0.3649267  0.6264785   0.583 0.563760    
# # # OrganHeart                             -0.1323569  0.6934961  -0.191 0.849683    
# # # OrganLiver                             -0.7087045  0.7927889  -0.894 0.377133    
# # # OrganMouthSwab                          1.7782226  0.7223642   2.462 0.018613 *  
# # # OrganSpleen                            -0.9332288  0.7835507  -1.191 0.241228    
# # # AmbientTemperatureCelsius              -0.1041644  0.0752376  -1.384 0.174507    
# # # PMI                                    -0.0000713  0.0146425  -0.005 0.996141    
# # # GenderM                                 0.9350903  0.4854873   1.926 0.061803 .  
# # # Age                                    -0.0289532  0.0255328  -1.134 0.264103    
# # # ---
# # # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # # (Dispersion parameter for quasibinomial family taken to be 0.01676009)

    # # # Null deviance: 2.58655  on 58  degrees of freedom
# # # Residual deviance: 0.64716  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 10

 	
# # otu_table_Proteobacteria<-otu_table(final_normalized)
# # Proteobacteria<-tax_table(final_normalized)[which(tax_table(final_normalized)[,"Phylum"]=="Proteobacteria"),]
# # otu_table_Proteobacteria<-otu_table_Proteobacteria[na.omit(match(rownames(Proteobacteria),
	# # rownames(otu_table_Proteobacteria))),]
# # data_glm<-cbind(sample_data(final_normalized),Proteobacteria=colSums(otu_table_Proteobacteria))

 # # summary(glm(Proteobacteria~CauseofDeath+Organ+AmbientTemperatureCelsius+PMI+Gender+Age,
 	# # family=quasibinomial,data=data_glm))
 	
 	
# # # Call:
# # # glm(formula = Proteobacteria ~ CauseofDeath + Organ + AmbientTemperatureCelsius + 
    # # # PMI + Gender + Age, family = quasibinomial, data = data_glm)

# # # Deviance Residuals: 
     # # # Min        1Q    Median        3Q       Max  
# # # -0.85592  -0.36562  -0.00368   0.27422   1.66468  

# # # Coefficients:
                                        # # # Estimate Std. Error t value Pr(>|t|)
# # # (Intercept)                             -8.16598   11.45404  -0.713    0.480
# # # CauseofDeathCar Accident                 8.88636   14.69717   0.605    0.549
# # # CauseofDeathChest Trauma                -3.03200   16.34218  -0.186    0.854
# # # CauseofDeathComplication after surgery  -9.40677  146.71377  -0.064    0.949
# # # CauseofDeathCoronary Heart Disease      -1.04284    5.08382  -0.205    0.839
# # # CauseofDeathDrowning                    -2.47094    9.56697  -0.258    0.798
# # # CauseofDeathDrug Overdose                2.07358    7.26278   0.286    0.777
# # # CauseofDeathGSW Chest                   -0.53486   16.83894  -0.032    0.975
# # # CauseofDeathGSW Head                     5.67234   13.42673   0.422    0.675
# # # CauseofDeathMultiple gunshot wounds      1.45397    4.86433   0.299    0.767
# # # CauseofDeathOverdose                    -0.16193    3.25421  -0.050    0.961
# # # CauseofDeathPending Further Studies      4.09030    7.74703   0.528    0.601
# # # CauseofDeathShotgun Wound of Head        0.21071   27.07865   0.008    0.994
# # # OrganBrain                               0.40432    3.77873   0.107    0.915
# # # OrganHeart                               0.30221    3.76020   0.080    0.936
# # # OrganLiver                               0.17213    3.78837   0.045    0.964
# # # OrganMouthSwab                          -2.19970    5.36027  -0.410    0.684
# # # OrganSpleen                              0.34441    3.71467   0.093    0.927
# # # AmbientTemperatureCelsius                0.03893    0.20333   0.191    0.849
# # # PMI                                     -0.05018    0.10196  -0.492    0.626
# # # GenderM                                 -3.12081    3.20490  -0.974    0.336
# # # Age                                      0.17960    0.22887   0.785    0.438

# # # (Dispersion parameter for quasibinomial family taken to be 5.385577)

    # # # Null deviance: 37.600  on 58  degrees of freedom
# # # Residual deviance: 16.791  on 37  degrees of freedom
# # # AIC: NA

# # # Number of Fisher Scoring iterations: 10 	
 	
# Using mG to predict obesity
#===================================================================================================
dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref")
featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
metadata<-updated_sampledata[match(colnames(matrixData),rownames(updated_sampledata)), ]
otus_metagenomeSeq<-newMRexperiment(matrixData, 
	phenoData = AnnotatedDataFrame(metadata), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeq_filter )
otus_metagenomeSeq_filter <- cumNorm(otus_metagenomeSeq_filter, p =p) 
otus_metagenomeSeq_filter_nor = MRcounts(otus_metagenomeSeq_filter , norm = TRUE, log = FALSE)
normFactor = normFactors(otus_metagenomeSeq_filter) 
normFactor = log2(normFactor/median(normFactor) + 1) 
AGE=pData(otus_metagenomeSeq_filter)$AGE
mod_obese_age_country= model.matrix(~0+country_obese + AGE) 
colnames(mod_obese_age_country)<-c("Ghana0","Ghana1","USA0","USA1","Age")
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(Ghana1 - Ghana0, USA1 - Ghana1, 
	USA1 - USA0, USA0 - Ghana0,
	(USA1 +USA0) - (Ghana1 - Ghana0), 
	(USA1 +Ghana1) - (USA0 - Ghana0),
	(USA1+USA0+Ghana1)/3-Ghana0,
	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(any(as.numeric(x)<0.001),TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_countryobese.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_countryobese.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
most_sig_OTUs<-cbind(rownames(tax_table_sig), rep("AGE_countryobese",nrow(tax_table_sig)))


AGE=pData(otus_metagenomeSeq_filter)$AGE
Country=pData(otus_metagenomeSeq_filter)$Country
BMI=pData(otus_metagenomeSeq_filter)$BMI
mod_obese_age_country= model.matrix(~BMI+AGE+Country) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(BMI, 	CountryUSA,	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_BMI.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_BMI.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]	
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_BMI",nrow(tax_table_sig))))
}		

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
updated_sampledata<-data.frame(sample_data(final_mGnorm),stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$Percentbodyfatmass),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$Percentbodyfatmass)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
Country=pData(otus_metagenomeSeqv1_filter)$Country
Percentbodyfatmass=pData(otus_metagenomeSeqv1_filter)$Percentbodyfatmass
mod_obese_age_country= model.matrix(~AGE+Country+Percentbodyfatmass) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(Percentbodyfatmass, CountryUSA,	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1]<0.001),TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_Percentbodyfatmass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_Percentbodyfatmass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]	
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig),rep("AGE_Percentbodyfatmass",nrow(tax_table_sig))))
}	
	
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
formicacid=pData(otus_metagenomeSeq_filter)$formicacid
mod_obese_age_country= model.matrix(~AGE+country_obese+formicacid) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(formicacid, country_obeseUSA0, levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_formicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_formicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]	
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig),rep("AGE_countryobese_formicacid",nrow(tax_table_sig))))
}	
	
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
aceticacid=pData(otus_metagenomeSeq_filter)$aceticacid
mod_obese_age_country= model.matrix(~AGE+country_obese+aceticacid) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(aceticacid,country_obeseUSA0,	levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_aceticacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_aceticacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_aceticacid",nrow(tax_table_sig))))
}	
	
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
propionicacid=pData(otus_metagenomeSeq_filter)$propionicacid
mod_obese_age_country= model.matrix(~AGE+country_obese+propionicacid) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(propionicacid, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_propionicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_propionicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
	cbind(rownames(tax_table_sig), rep("AGE_countryobese_propionicacid",nrow(tax_table_sig))))
}	
	
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
butyricacid=pData(otus_metagenomeSeq_filter)$butyricacid
mod_obese_age_country= model.matrix(~AGE+country_obese+butyricacid) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(butyricacid,country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_butyricacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_butyricacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
cbind(rownames(tax_table_sig), rep("AGE_countryobese_byturicacid",nrow(tax_table_sig))))
}	
	
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
isovalericacid=pData(otus_metagenomeSeq_filter)$isovalericacid
mod_obese_age_country= model.matrix(~AGE+country_obese+isovalericacid) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(isovalericacid,country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_isovalericacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_isovalericacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_isovalericacid",nrow(tax_table_sig))))
}	
	
 GlucoseResultGLUCRSLT adiponectin  leptin insulin glucosemmol insulinmmol   homa_ir    ratio
AGE=pData(otus_metagenomeSeq_filter)$AGE
country_obese=pData(otus_metagenomeSeq_filter)$country_obese
totscfa=pData(otus_metagenomeSeq_filter)$totscfa
mod_obese_age_country= model.matrix(~totscfa+AGE+country_obese) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeq_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(totscfa, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_totscfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_totscfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_totscfa",nrow(tax_table_sig))))
}		
	
adiponectin  leptin insulin glucosemmol insulinmmol   homa_ir    ratio
featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$GlucoseResultGLUCRSLT),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$GlucoseResultGLUCRSLT)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
GlucoseResultGLUCRSLT=pData(otus_metagenomeSeqv1_filter)$GlucoseResultGLUCRSLT
mod_obese_age_country= model.matrix(~AGE+country_obese+GlucoseResultGLUCRSLT) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(GlucoseResultGLUCRSLT, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_GlucoseResultGLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_GlucoseResultGLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_GlucoseResultGLUCRSLT",nrow(tax_table_sig))))
}	
	

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$adiponectin ),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$adiponectin )]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
adiponectin =pData(otus_metagenomeSeqv1_filter)$adiponectin 
mod_obese_age_country= model.matrix(~AGE+country_obese+adiponectin ) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(adiponectin , country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_adiponectin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_adiponectin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_adiponectin",nrow(tax_table_sig))))
}	
	

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$leptin),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$leptin)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
leptin=pData(otus_metagenomeSeqv1_filter)$leptin 
mod_obese_age_country= model.matrix(~AGE+country_obese+leptin ) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(leptin , country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_leptin",nrow(tax_table_sig))))
}	
	

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$insulin),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$insulin)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
insulin=pData(otus_metagenomeSeqv1_filter)$insulin 
mod_obese_age_country= model.matrix(~AGE+country_obese+insulin) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(insulin, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_insulin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_insulin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig),rep("AGE_countryobese_insulin",nrow(tax_table_sig))))
}	

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$glucosemmol),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$glucosemmol)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
glucosemmol=pData(otus_metagenomeSeqv1_filter)$glucosemmol
mod_obese_age_country= model.matrix(~AGE+country_obese+glucosemmol) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(glucosemmol, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_glucosemmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_glucosemmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_glucosemmol",nrow(tax_table_sig))))
}	
featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$insulinmmol),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$insulinmmol)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
insulinmmol=pData(otus_metagenomeSeqv1_filter)$insulinmmol
mod_obese_age_country= model.matrix(~AGE+country_obese+insulinmmol) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(insulinmmol, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_insulinmmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_insulinmmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_insulinmmol",nrow(tax_table_sig))))
}	
	
featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$homa_ir ),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$homa_ir)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
homa_ir=pData(otus_metagenomeSeqv1_filter)$homa_ir
mod_obese_age_country= model.matrix(~AGE+country_obese+homa_ir) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(homa_ir, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_homa_ir.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_homa_ir.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_homa_ir",nrow(tax_table_sig))))
}	

featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
updated_sampledata_v1<-updated_sampledata[!is.na(updated_sampledata$ratio),]
matrixData_v1<-matrixData[,!is.na(updated_sampledata$ratio)]
metadatav1<-updated_sampledata_v1[match(colnames(matrixData_v1),rownames(updated_sampledata_v1)), ]
otus_metagenomeSeqv1<-newMRexperiment(matrixData_v1, 
	phenoData = AnnotatedDataFrame(metadatav1), 
	featureData = AnnotatedDataFrame(featureData))
otus_metagenomeSeqv1_filter = filterData(otus_metagenomeSeqv1, present = round(97/4*0.5), depth = 1) 
p = cumNormStatFast(otus_metagenomeSeqv1_filter )
otus_metagenomeSeqv1_filter <- cumNorm(otus_metagenomeSeqv1_filter, p =p) 
AGE=pData(otus_metagenomeSeqv1_filter)$AGE
country_obese=pData(otus_metagenomeSeqv1_filter)$country_obese
ratio=pData(otus_metagenomeSeqv1_filter)$ratio
mod_obese_age_country= model.matrix(~AGE+country_obese+ratio) 
settings = zigControl(maxit = 100, verbose = TRUE, pvalMethod = "bootstrap") 
fit_obese_age_country = fitZig(obj = otus_metagenomeSeqv1_filter, mod = mod_obese_age_country, useCSSoffset =TRUE, control = settings)
results_obese_age_country<-MRfulltable(fit_obese_age_country, number = nrow(assayData(otus_metagenomeSeq_filter)$counts))
# The output of fitZig contains a list of various useful # items. hint: names(res). Probably the most useful is the # limma 'MLArrayLM' object called fit. 
zigFit = fit_obese_age_country$fit 
finalMod = fit_obese_age_country$fit$design 
contrast.matrix = makeContrasts(ratio, country_obeseUSA0,levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix) 
fit2 = eBayes(fit2) 
topTable(fit2)
pvalues<-apply(fit2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit2$coefficients,pvalues)
rownames(otus_mG)<-rownames(fit2$coefficients)
otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2+1):ncol(otus_mG)],1,function(x){ifelse(as.numeric(x[1])<0.001,TRUE,FALSE)}),]
write.table(otus_mG_filtered,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/otus_mG_filtered_Country_ratio.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_sig<-merge(featureData[match(rownames(otus_mG_filtered),rownames(featureData)),],otus_mG_filtered,by.x=0,by.y=0)
write.table(tax_table_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/tax_table_sig_Country_ratio.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
rownames(tax_table_sig)<-tax_table_sig[,1]
if (!is.null(dim(tax_table_sig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_sig), rep("AGE_countryobese_ratio",nrow(tax_table_sig))))
}	
write.table(most_sig_OTUs,file="/Users/beatrizp/Dropbox/SCFA_Obesity/mG_openref/most_sig_OTUs.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	

# Using DESeq to predict obesity
#===================================================================================================
dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref")
aaa<-data.frame(Country=sample_data(final_filtered)[,"Country"],obese=sample_data(final_filtered)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")

metadata<-sample_data(data.frame(sample_data(final_filtered), country_obese=country_obese))
final_filtered=phyloseq(otu_table(final_filtered),metadata,tax_table(final_filtered),tree)


final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ AGE+country_obese)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
res1 <- results(final_filtered2DESeq_res,contrast=c("country_obese", "Ghana1","Ghana0"))
# Use log2FoldChange and padj
res2 <- results(final_filtered2DESeq_res,contrast=c("country_obese", "USA1","Ghana1"))
res3 <- results(final_filtered2DESeq_res,contrast=c("country_obese", "USA1","USA0"))
res6 <- results(final_filtered2DESeq_res,contrast=c("country_obese", "USA0","Ghana0"))
res4 <- results(final_filtered2DESeq_res,contrast= c(0,0,-1,1,-1,1))
res5 <- results(final_filtered2DESeq_res,contrast=c(0,0,-1,-1,1,1))
res7<-results(final_filtered2DESeq_res,contrast=c(0,0,-1,1/3,1/3,1/3))
res_final<-cbind(ghana_coeff=res1[, "log2FoldChange"], obese_coeff=res2[,"log2FoldChange"],
	usa_coeff=res3[,"log2FoldChange"], lean_coeff=res6[,"log2FoldChange"], country_coeff=res4[,"log2FoldChange"],
	obesityglobal_coeff=res5[,"log2FoldChange"], allvsleanghana_coeff=res7[,"log2FoldChange"],ghana_padj=res1[, "padj"], obese_padj=res2[,"padj"], 
	usa_padj=res3[,"padj"], lean_padj=res6[,"padj"], country_padj=res4[,"padj"],obesityglobal_padj=res5[,"padj"],
	allvsleanghana_padj=res7[,"padj"])



results_all<-cbind(data.frame(mcols(final_filtered2DESeq_res), stringAsFactors=FALSE) ,res_final)
rownames(results_all)<-rownames(res1)
#rownames(results_all)<-rownames(results_final)
aaa<-results_all[,grep("_padj", colnames(results_all))]
results_all_filtered<-results_all[apply(aaa,1,function(x){
	ifelse(any(as.numeric(x[!is.na(x)])<=0.001),TRUE,FALSE)}),]
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[match(rownames(results_psig),rownames(otu_table_psig)),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_obeseAGE.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]

tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_obeseAGE.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
most_sig_OTUs<-c()
compids<-grep("padj", colnames(tax_table_psig))
if (!is.null(dim(tax_table_psig))){
	for (ids in compids){
		if (is.null(most_sig_OTUs))
			most_sig_OTUs<-cbind(OTUs=rownames(tax_table_psig), padj=tax_table_psig [, ids],
			variable=rep(strsplit(colnames(tax_table_psig)[ids],"_padj",fixed=TRUE)[[1]][1],nrow(tax_table_psig)))
		else {
			most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig [, ids],
					variable=rep(strsplit(colnames(tax_table_psig)[ids],"_padj",fixed=TRUE)[[1]][1],nrow(tax_table_psig))))
		}
	}	
}



# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ AGE+BMI)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
# results_final<-results(final_filtered2DESeq_res)
# results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
# #rownames(results_all)<-rownames(results_final)

# results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_BMIAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_BMIAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)

# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_BMI",nrow(tax_table_psig))))
# }	
	
# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
# results_final<-results(final_filtered2DESeq_res)
# results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
# #rownames(results_all)<-rownames(results_final)

# results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_Country.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_Country.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("us",nrow(tax_table_psig))))
# }	
	
# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us +AGE+BMI)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric",minReplicatesForReplace=3)
# results_final<-results(final_filtered2DESeq_res)
# results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
# #rownames(results_all)<-rownames(results_final)

# results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
# testOTUs<-rownames(results_all_filtered)

# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_CountryobeseAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_CountryobeseAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("us_AGE_BMI",nrow(tax_table_psig))))
# }


# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# data_psig0<-data_psig
# for (iOTU in rownames(results_psig)){
	# temp<-data_psig[,iOTU]-exp(tax_table_psig[which(rownames(tax_table_psig)==iOTU),"Intercept"]+
	# data_psig[,"us"]*tax_table_psig[which(rownames(tax_table_psig)==iOTU),"us"]+
	# data_psig[,"AGE"]*tax_table_psig[which(rownames(tax_table_psig)==iOTU),"AGE"])
	# data_psig<-cbind(data_psig,temp)
# }
# colnames(data_psig)<-c(colnames(data_psig0),paste(rownames(results_psig),"_adj",sep=""))
# data_psig_corr<-data_psig


# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ Country +AGE+BMI)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
# results_final<-results(final_filtered2DESeq_res)
# results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
# #rownames(results_all)<-rownames(results_final)

# results_all_filtered<-results_all[which((results_final$padj<=0.1)==TRUE),]	
# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_CountryBMIAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_CountryBMIAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("Country_AGE_BMI",nrow(tax_table_psig))))
# }

# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ Country +AGE+obese)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
# results_final<-results(final_filtered2DESeq_res)
# results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
# #rownames(results_all)<-rownames(results_final)

# results_all_filtered<-results_all[which((results_final$padj<=0.1)==TRUE),]	
# results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_CountryobeseAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig,results_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_CountryobeseAGE.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("Country_AGE_obese",nrow(tax_table_psig))))
# }

# sample_new<-sample_data(final_filtered);
# Country_obese<-paste(sample_new$Country,sample_new$obese,collapsed="",sep="")
# sample_new<-cbind(sample_new,Country_obese)
# sample_new<-sample_data(sample_new);
# final_filtered = merge_phyloseq(otu_table_final_filtersamples_OTUs, sample_new,tax_split,tree)
# final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ AGE+Country_obese)
# gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
# geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
# final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
# final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
# results_final_ghana<-results(final_filtered2DESeq_res,
	# contrast=c("Country_obese","Ghana1", "Ghana0"))
# results_final_usa<-results(final_filtered2DESeq_res,
	# contrast=c("Country_obese","USA1", "USA0"))
# results_final_obese<-results(final_filtered2DESeq_res,
	# contrast=c("Country_obese","USA1", "Ghana1"))
# results_final_lean<-results(final_filtered2DESeq_res,
	# contrast=c("Country_obese","USA0", "Ghana0"))
# results_final_Ghanalean<-results(final_filtered2DESeq_res,
	# contrast=list("Country_obeseGhana0",c("Country_obeseGhana1","Country_obeseUSA0","Country_obeseUSA1")),
	# listValues=c(1, -1/3))


# results_all_ghana<-cbind(mcols(final_filtered2DESeq_res), results_final_ghana)
# #rownames(results_all)<-rownames(results_final)

# results_all_ghana_filtered<-results_all_ghana[which((results_final_ghana$padj<=0.1)==TRUE),]	
# results_all_ghana_psig<-data.frame(results_all_ghana_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_all_ghana_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_AGEGhana.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig, results_all_ghana_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_AGEGhana.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_Ghana",nrow(tax_table_psig))))
# }

# results_all_USA<-cbind(mcols(final_filtered2DESeq_res), results_final_USA)
# results_all_USA_filtered<-results_all_USA[which((results_final_USA$padj<=0.1)==TRUE),]	
# results_all_USA_psig<-data.frame(results_all_USA_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_all_USA_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_AGEUSA.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig, results_all_USA_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_AGEUSA.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_USA",nrow(tax_table_psig))))
# }

# results_all_obese<-cbind(mcols(final_filtered2DESeq_res), results_final_obese)
# results_all_obese_filtered<-results_all_obese[which((results_final_obese$padj<=0.1)==TRUE),]	
# results_all_obese_psig<-data.frame(results_all_obese_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_all_obese_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_AGE_obeseBothCountries.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig, results_all_obese_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_AGE_obeseBothCountries.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_obeseBothCountries",nrow(tax_table_psig))))
# }

# results_all_lean<-cbind(mcols(final_filtered2DESeq_res), results_final_lean)
# results_all_lean_filtered<-results_all_lean[which((results_final_lean$padj<=0.1)==TRUE),]	
# results_all_lean_psig<-data.frame(results_all_lean_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_all_lean_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_AGE_leanBothCountries.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig, results_all_lean_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_AGE_leanBothCountries.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_leanBothCountries",nrow(tax_table_psig))))
# }

# results_all_Ghanalean<-cbind(mcols(final_filtered2DESeq_res), results_final_Ghanalean)
# results_all_Ghanalean_filtered<-results_all_Ghanalean[which((results_final_Ghanalean$padj<=0.1)==TRUE),]	
# results_all_Ghanalean_psig<-data.frame(results_all_Ghanalean_filtered,stringsAsFactors=FALSE)
# otu_table_psig<-otu_table(final_filtered)
# otu_table_psig<-otu_table_psig[rownames(results_all_Ghanalean_psig),]
# sample_data_psig<-sample_data(final_filtered)
# data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
# write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_AGEGhanalean.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# tax_table_psig<-tax_table(final_filtered)
# tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	# rownames(tax_table(final_filtered)))),]
# tax_table_psig<-cbind(tax_table_psig, results_all_Ghanalean_psig)
# write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_AGEGhanalean.txt",
	# sep="\t",col.names=TRUE,row.names=TRUE)
# if (!is.null(dim(tax_table_psig))){
	# most_sig_OTUs<-rbind(most_sig_OTUs, 
		# cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			# rep("AGE_Ghanalean",nrow(tax_table_psig))))
# }

sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"leptin"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+leptin )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.1)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs, 
		cbind(rownames(tax_table_psig), tax_table_psig$pvalue, tax_table_psig$padj,
			rep("AGE_Country_BMI_leptin",nrow(tax_table_psig))))
}
	
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"Percentbodyfatmass"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+Country+Percentbodyfatmass )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_Percentbodyfatmass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_Percentbodyfatmass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_Percentbodyfat",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"waist"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+Country+waist )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_waist.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_waist.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_waist",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"meanenergy"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+Country+BMI+meanenergy)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_meanenergy.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_meanenergy.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_meanenergy",nrow(tax_table_psig))))
}
	
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"meanfat"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+meanfat)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_meanfat.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_meanfat.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_meanfat",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"GlucoseResultGLUCRSLT"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+GlucoseResultGLUCRSLT)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_GlucoseResultGLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_GlucoseResultGLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_GlucoseResultGLUCRSLT",nrow(tax_table_psig))))
}


sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"totscfa"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+totscfa)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_totscfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_totscfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_totscfa",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"formicacid"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+formicacid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_formicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_formicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_formicacid",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"aceticacid"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+formicacid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_aceticacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_aceticacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_aceticacid",nrow(tax_table_psig))))
}

		
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"propionicacid"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+propionicacid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_propionicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_propionicacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_propionicacid",nrow(tax_table_psig))))
}


sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"butyricacid"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+butyricacid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_butyricacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_butyricacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_butyricacid",nrow(tax_table_psig))))
}


sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"isovalericacid"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+isovalericacid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_isovalericacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_isovalericacid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_isovalericacid",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa4"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa4)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa4.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa4.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa4",nrow(tax_table_psig))))
}

	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa6"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa6)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa6.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa6.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	

if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa6",nrow(tax_table_psig))))
}	

sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa8"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa8)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa8.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa8.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa8",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa10"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa10)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa10.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa10.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa10",nrow(tax_table_psig))))
}

sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa12"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa12)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa12.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa12.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa12",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa14"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa14)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa14.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa14.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa14",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa16"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa16)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa16.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa16.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa16",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa17"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa17)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa17.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa17.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa17",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa18"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa18)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa18.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa18.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa18",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa20"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa20)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa20.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa20.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa20",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"sfa22"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+sfa22)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_sfa22.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_sfa22.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_sfa22",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"glucrslt"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+glucrslt)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_glucrslt.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_glucrslt.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Countryobese_glucrslt",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"leptin"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+leptin)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_leptin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Countryobese_leptin",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"insulin"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+insulin)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_insulin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_insulin.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	

if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_insulin",nrow(tax_table_psig))))
}
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"glucosemmol"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+glucosemmol)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_glucosemmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_glucosemmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_glucosemmol",nrow(tax_table_psig))))
}	
	
sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"insulinmmol"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+insulinmmol)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_insulinmmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_insulinmmol.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	

if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_insulinmmol",nrow(tax_table_psig))))
}

sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"homa_ir"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+homa_ir)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_homa_ir.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_homa_ir.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_homa_ir",nrow(tax_table_psig))))
}

sampledataf<-sample_data(final_filtered)[!is.na(sample_data(final_filtered)[,"ratio"]),]
otudataf<-otu_table(final_filtered)
taxf<-tax_table(final_filtered)
finalv2 = merge_phyloseq(otudataf,sampledataf, taxf,tree)
final_filtered2DESeq = phyloseq_to_deseq2(finalv2, ~ AGE+country_obese+ratio)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.15)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(finalv2)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(finalv2)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/mapfile_pvalues_homa_ir.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(finalv2)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(finalv2)))),]
tax_table_psig<-cbind(tax_table_psig, results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/taxtable_pvalues_ratio.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)	
	
if (!is.null(dim(tax_table_psig))){
	most_sig_OTUs<-rbind(most_sig_OTUs,
				cbind( OTUs=rownames(tax_table_psig), padj=tax_table_psig$padj,
					variable=rep("AGE_Country_BMI_ratio",nrow(tax_table_psig))))
}

 most_sig_OTUs_merged<-merge(tax_table(final_filtered),most_sig_OTUs,by.x=0,by.y=1)
 most_sig_OTUs_filtered<-most_sig_OTUs_merged[as.numeric(as.character(most_sig_OTUs_merged$V3))<0.001,]
 most_sig_OTUs_filtered<-most_sig_OTUs_filtered[order(most_sig_OTUs_filtered$V4),]
 
 write.table(most_sig_OTUs_filtered, 
 	file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq_openref/most_sig_OTUs_filtered.txt", sep="\t",
 	row.names=FALSE)
 
# # Try unifrac with only those variables
# #=================================================================================================

# otu_table<-otu_table(final)

# data_psig_corr2<-data_psig_corr[,grep("OTU", colnames(data_psig_corr))]
# data_psig_corr2<-apply(data_psig_corr2,c(1,2),function(x){
	# rr<-ifelse(as.numeric(x)<=0, 0,round(as.numeric(x),0))
	# return(rr)
# })
# data_psig_corr3<-data_psig_corr2[,grep("_adj", colnames(data_psig_corr2))]
# colnames(data_psig_corr3)<-gsub("_adj", "",colnames(data_psig_corr3))
# data_psig_corr3<-t(data_psig_corr3)
# data_psig_corr4<-matrix(data_psig_corr3, ncol=ncol(data_psig_corr3))
# colnames(data_psig_corr4)<-colnames(data_psig_corr3)
# rownames(data_psig_corr4)<-rownames(data_psig_corr3)
# data_psig_corr4<-otu_table(data_psig_corr4, taxa_are_rows=TRUE)
# final_DESeq<-merge_phyloseq(data_psig_corr4, mapfile,tax_split,tree)

# unifrac_weighted_normalized=UniFrac(final_DESeq, weighted=TRUE, normalized=TRUE, parallel=FALSE, 
	# fast=TRUE)
	
# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/UniFrac_Weigthed_Normalized_DESeq.pdf")
# weighted_n.pco2 = ordinate(final_DESeq, method = "PCoA", distance = unifrac_weighted_normalized)
# p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
# p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/Bray_DESeq.pdf")
# weighted_n.pco2 = ordinate(final_DESeq, method = "PCoA", distance = "bray")
# p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
# p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/UniFrac_Weigthed_Normalized_NMDS_DESeq.pdf")
# weighted_n.pco2 = ordinate(final_DESeq, method = "NMDS", distance = unifrac_weighted_normalized)
# p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
# p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()

# pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/Bray_NMDS_DESeq.pdf")
# weighted_n.pco2 = ordinate(final_DESeq, method = "NMDS", distance = "bray")
# p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "obese", shape = "Country") 
# p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# dev.off()


# # otu_table<-otu_table(final)
# # otu_table_filter_BMI<-otu_table[match(BMI_corrected_OTUs,rownames(otu_table)),]
# # otu_table_filter_BMI<-otu_table_filter_BMI[,colSums(otu_table_filter_BMI)>0]
# # final_DESeq_BMI<-merge_phyloseq(otu_table_filter_BMI, mapfile,tax_split,tree)
# # unifrac_weighted_normalized=UniFrac(final_DESeq_BMI, weighted=TRUE, normalized=TRUE, parallel=FALSE, 
	# # fast=TRUE)
	
# # pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/UniFrac_Weigthed_Normalized_DESeqBMI.pdf")
# # weighted_n.pco2 = ordinate(final_DESeq_BMI, method = "PCoA", distance = unifrac_weighted_normalized)
# # p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "BMI", shape = "Country") 
# # p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# # dev.off()

# # pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/UniFrac_Weigthed_Normalized_NMDS_DESeqBMI.pdf")
# # weighted_n.pco2 = ordinate(final_DESeq_BMI, method = "NMDS", distance = unifrac_weighted_normalized)
# # p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "BMI", shape = "Country") 
# # p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# # dev.off()

# # pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/bray_DESeqBMI.pdf")
# # weighted_n.pco2 = ordinate(final_DESeq_BMI, method = "PCoA", distance = "bray")
# # p2.2 = plot_ordination(final_DESeq, weighted_n.pco2, type = "samples", color = "BMI", shape = "Country") 
# # p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# # dev.off()


# # pdf(file="/Users/beatrizp/Dropbox/SCFA_Obesity/BetaDiversity/bray_NMDS_DESeqBMI.pdf")
# # weighted_n.pco2 = ordinate(final_DESeq_BMI, method = "NMDS", distance = "bray")
# # p2.2 = plot_ordination(final_DESeq_BMI, weighted_n.pco2, type = "samples", color = "BMI", shape = "Country") 
# # p2.2 + theme_bw() + theme(text = element_text(size = 8)) + geom_point(size = 3) + scale_colour_gradient(low="red", high="blue") +
        # # guides(color = guide_legend(order=2, override.aes = list(shape = "__", size = 10))) 
# # dev.off()


# Random forest code
#=================================================================================================
featureData =data.frame(tax_table(final_mGnorm))
matrixData<-matrix(otu_table(final_mGnorm),ncol=ncol(otu_table(final_mGnorm)))
rownames(matrixData)<-rownames(otu_table(final_mGnorm))
colnames(matrixData)<-colnames(otu_table(final_mGnorm))
aaa<-data.frame(Country=sample_data(final_mGnorm)[,"Country"],obese=sample_data(final_mGnorm)[,"obese"])
country_obese=paste(aaa[,"Country"],aaa[,"obese"],sep="")
updated_sampledata<-data.frame(sample_data(final_mGnorm),country_obese,stringsAsFactors=FALSE)
metadata<-updated_sampledata[match(colnames(matrixData),rownames(updated_sampledata)), ]

library(randomForest)
library(mlbench)
library(caret)
 
# Load Dataset
data(Sonar)
dataset <- Sonar
x <- dataset[,1:60]
y <- dataset[,61]

# Create model with default paramters
control <- trainControl(method="repeatedcv", number=10, repeats=3)
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)

# Random Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
set.seed(seed)
mtry <- sqrt(ncol(x))
rf_random <- train(Class~., data=dataset, method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)



set.seed(647) 
mG_RF <- cbind(sample_data(final_mGnorm)[,"country_obese"], t(matrixData)) 
result <- rfcv(mG_RF, mG_RF$country_obese, cv.fold=3) 
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2)) 
## The following can take a while to run, so if you really want to try 
## it, copy and paste the code into R. 
## Not run: 
result <- replicate(5, rfcv(myiris, iris$Species), simplify=FALSE) 
error.cv <- sapply(result, "[[", "error.cv") 
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l", lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
	xlab="Number of variables", ylab="CV Error") #


# XGBoost code with only the variables from DESeq
#=================================================================================================
dir.create("/Users/beatrizp/Dropbox/SCFA_Obesity/XGBoost/")
otu_table_xgboost<-otu_table(final_DESeq)
sample_data_xgboost<-sample_data(final_DESeq)
data_xgboost<-data.frame(sample_data_xgboost,
	t(otu_table_xgboost[,match(rownames(sample_data_xgboost),
		colnames(otu_table_xgboost))]),
	stringsAsFactors=FALSE)

# one-hot-encoding categorical features
ohe_feats = c('us','obese')
dt_xgboost<-data.table(data_xgboost)
dummies <- dummyVars(~   us+obese, data = dt_xgboost)
df_all_ohe <- as.data.frame(predict(dummies, newdata = dt_xgboost))
data_xgboostv2<-as.data.frame(dt_xgboost)
binary_OTU<-apply(data_xgboostv2[,grep("OTU",colnames(data_xgboostv2))],c(1,2),function(x){ifelse (as.numeric(x)>5,1,0)})
data_xgboost_cat_fil<-cbind(AGE=data_xgboostv2[,"AGE"], BMI=data_xgboostv2[,"BMI"], df_all_ohe, binary_OTU)

#Convert numeric data into factors
dt_xgboostv2<-data.table(data_xgboost_cat_fil)
dt_xgboostv2[,AGE_discrete:=as.factor(round(AGE/10,0))]
output_vector=dt_xgboostv2[,AGE_discrete]
dt_xgboostv2[,AGE:=NULL]
dt_xgboostv2[,BMI_discrete:=as.factor(round(BMI/10,0))]
output_vector=dt_xgboostv2[,BMI_discrete]
dt_xgboostv2[,BMI:=NULL]

dt_xgboost_df<-data.frame(dt_xgboostv2)

bootstrapping_results<-c()
nboot=100
for (iboot in 1:nboot){
	dt_xgboost_iboot<-dt_xgboost_df[sample(seq(1,nrow(dt_xgboost_df),by=1),
		round(0.7*nrow(dt_xgboost_df)),replace=FALSE),]
	dt_xgboost_iboot<-data.table(dt_xgboost_iboot)
	sparse_matrix_iboot <- sparse.model.matrix(obese ~ .-1, data = dt_xgboost_iboot)
	output_vector_iboot=dt_xgboost_iboot[,obese]
	bst_iboot <- xgboost(data = sparse_matrix_iboot , label = output_vector_iboot, 
		max.depth = 4, eta = 1, nthread = 2, nround = 1000)
	model_iboot <- xgb.dump(bst_iboot, with.stats = T)
	names_iboot <- dimnames(dt_xgboost_iboot)[[2]]
	importance_matrix_iboot <- xgb.importance(names_iboot, model = bst_iboot)
	importance_matrix_top<-data.frame(importance_matrix_iboot)
	if (iboot==1){
		bootstrapping_results<-cbind(trial=rep(iboot,nrow(importance_matrix_top)),
			importance_matrix_top)
	} else {
		bootstrapping_results<-rbind(bootstrapping_results,
			cbind(trial=rep(iboot,nrow(importance_matrix_top)),importance_matrix_top))
	}	
}

rank_variables<-tapply(bootstrapping_results$Gain,bootstrapping_results$Feature,sum)/nboot
rank_variables<-rank_variables[order(rank_variables,decreasing=TRUE)]
pdf("/Users/beatrizp/Dropbox/SCFA_Obesity/XGBoost/barplots_boostrapping")
barplot(rank_variables[1:5],names.arg=names(rank_variables)[1:5])
dev.off()

#Make model with top more robust predictors
	dt_xgboost_df_top<-dt_xgboost_df[,match(names(rank_variables),colnames(dt_xgboost_df))]
	data_xgboost_df_top_train<-dt_xgboost_df_top[sample(seq(1,nrow(dt_xgboost_df_top),by=1),
		round(0.7*nrow(dt_xgboost_df_top)),replace=FALSE),]
	data_xgboost_df_top_test<-dt_xgboost_df_top[is.na(match(rownames(dt_xgboost_df_top),
		rownames(data_xgboost_df_top_train))),]
	dt_xgboost_df_top_train<-data.table(data_xgboost_df_top_train)
	dt_xgboost_df_top_test<-data.table(data_xgboost_df_top_test)
	sparse_matrix<- sparse.model.matrix(obese ~ .-1, data = data_xgboost_df_top_train)
	output_vector=dt_xgboost_df_top_train[,obese]
	bst_top <- xgboost(data = sparse_matrix , 
		label = output_vector, 
		max.depth =10, 
		eta = 0.1, 
		nthread = 4, 
		nround = 1000,
		subsample=2,
		nthread=4,
		gamma=0.01,
		min_child_weight=10,
		colsample_bytree=12,
		eval_metric = "auc",
		alpha=0.5)
	pred<-predict(bst_top, data.matrix(data_xgboost_df_top_test[,-1]))
	cor(round(pred,0),as.numeric(dt_xgboost_df_top_test[,obese]))
	prediction_table<-cbind(round(pred,0),dt_xgboost_df_top_test[,obese], dt_xgboost_df_top_test[,us])
	plot(pred,as.numeric(dt_xgboost_df_top_test[,obese]))
	TP=length(prediction_table[(prediction_table[,2]==1 & prediction_table [,1]==1),1])
	TN=length(prediction_table[(prediction_table[,2]==0 & prediction_table [,1]==0),1])
	FN=length(prediction_table[(prediction_table[,2]==1 & prediction_table [,1]==0),1])
	FP=length(prediction_table[(prediction_table[,2]==0 & prediction_table [,1]==1),1])
	Sensitivity=TP/(TP+FN)
	Specificity=TN/(TN+FP)
	Precision=TP/(TP+FP)
	model_top <- xgb.dump(bst_top, with.stats = T)
	names_top <- dimnames(dt_xgboost_df_top_train)[[2]][-1]
	importance_matrix_top <- xgb.importance(names_top, model = bst_top)
	importance_matrix_top<-data.frame(importance_matrix_top)
	pdf("/Users/beatrizp/Dropbox/SCFA_Obesity/XGBoost/barplots_topvariables")
	barplot(as.numeric(importance_matrix_top$Gain),names.arg=importance_matrix_top$Feature)
	dev.off()
    taxa_filtered<-tax_table(final_DESeq)
    taxa_sig<-taxa_filtered[na.omit(match(importance_matrix_top$Feature,rownames(taxa_filtered))),]
    write.table(taxa_sig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/XGBoost/top_train_taxonomy",
    	sep="\t",col.names=TRUE,row.names=TRUE)

# Predict other significant variables
#==========================================================================	
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"Glucose_Result_GLUCRSLT"]),]
final_filteredGLU= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)

final_filtered2DESeq = phyloseq_to_deseq2(final_filteredGLU, ~ us+ Glucose_Result_GLUCRSLT )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res,contrasts=c("Glucose_Result_GLUCRSLT","B","C"))
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.9)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_Glucose_Result_GLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_Glucose_Result_GLUCRSLT.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"acetic_acid"]),]
interval<-(max(as.numeric(unlist(sample_data_filter[,"acetic_acid"])))-min(as.numeric(unlist(sample_data_filter[,"acetic_acid"]))))/4
sample_data_filter[,"acetic_acid"]<-round(sample_data_filter[,"acetic_acid"]/interval,0)

final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us+ acetic_acid  )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_acetic_acid .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_acetic_acid .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us+ formic_acid  )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.9)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_formic_acid .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_formic_acid .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us+ propionic_acid  )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_propionic_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_propionic_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us+ butyric_acid  )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_butyric_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_butyric_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered, ~ us+ isovaleric_acid)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_isovaleric_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_isovaleric_acid.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa"]),]
final_filteredsfa= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filteredsfa, ~ us+ mean_sfa)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
		
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"Mean_duration_mins_of_moderate_activity_bouts_minimum_1_min_bouts_per_day_fr"]),]
final_filteredexer= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filteredexer, ~ us+ Mean_duration_mins_of_moderate_activity_bouts_minimum_1_min_bouts_per_day_fr)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.9999)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_Mean_duration_mins_of_moderate_activity_bouts_minimum_1_min_bouts_per_day_fr.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_Mean_duration_mins_of_moderate_activity_bouts_minimum_1_min_bouts_per_day_fr.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa4"]),]
final_filteredexer= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filteredexer, ~ us+ mean_sfa4)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa4.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa4.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
	
sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa6"]),]
final_filtered6= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered6, ~ us+ mean_sfa6)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa6.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa6.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa8"]),]
final_filtered8= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered8, ~ us+ mean_sfa8)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa8.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa8.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa10"]),]
final_filtered10= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered10, ~ us+ mean_sfa10)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa10.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa10.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa12"]),]
final_filtered12= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered12, ~ us+ mean_sfa12)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa12.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa12.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa14"]),]
final_filtered14= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered14, ~ us+ mean_sfa14)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa14.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa14.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa16"]),]
final_filtered16= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered16, ~ us+ mean_sfa16)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa16.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa16.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa17"]),]
final_filtered17= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered17, ~ us+ mean_sfa17)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa17.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa17.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa18"]),]
final_filtered18= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered18, ~ us+ mean_sfa18)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.5)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa18.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa18.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)


sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa20"]),]
final_filtered20= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered20, ~ us+ mean_sfa20)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.9)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa20.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa20.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"mean_sfa22"]),]
final_filtered22= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filtered22, ~ us+ mean_sfa22)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.9)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_mean_sfa22.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_mean_sfa22.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"tot_scfa"]),]
final_filteredtot_scfa= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filteredtot_scfa , ~ us+ tot_scfa )
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_tot_scfa .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_tot_scfa .txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

sample_data_filter<-sample_data(final_filtered)
sample_data_filter<-sample_data_filter[!is.na(sample_data_filter[,"Percent_body_fat_mass"]),]
final_filteredbodyfat= merge_phyloseq(otu_filtered, sample_data_filter,tax_split,tree)
	
final_filtered2DESeq = phyloseq_to_deseq2(final_filteredbodyfat, ~ us+ Percent_body_fat_mass)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(final_filtered2DESeq), 1, gm_mean)
final_filtered2DESeq = estimateSizeFactors(final_filtered2DESeq, geoMeans = geoMeans)
final_filtered2DESeq_res= DESeq(final_filtered2DESeq, test="Wald", fitType="parametric")
results_final<-results(final_filtered2DESeq_res)
results_all<-cbind(mcols(final_filtered2DESeq_res),results_final)
#rownames(results_all)<-rownames(results_final)

results_all_filtered<-results_all[which((results_final$padj<=0.01)==TRUE),]	
results_psig<-data.frame(results_all_filtered,stringsAsFactors=FALSE)
otu_table_psig<-otu_table(final_filtered)
otu_table_psig<-otu_table_psig[rownames(results_psig),]
sample_data_psig<-sample_data(final_filtered)
data_psig<-cbind(sample_data_psig,t(otu_table_psig[,match(rownames(sample_data_psig),colnames(otu_table_psig))]))
write.table(data_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/mapfile_pvalues_Percent_body_fat_mass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)
tax_table_psig<-tax_table(final_filtered)
tax_table_psig<-tax_table_psig[na.omit(match(rownames(otu_table_psig),
	rownames(tax_table(final_filtered)))),]
tax_table_psig<-cbind(tax_table_psig,results_psig)
write.table(tax_table_psig,file="/Users/beatrizp/Dropbox/SCFA_Obesity/DESeq/taxtable_pvalues_Percent_body_fat_mass.txt",
	sep="\t",col.names=TRUE,row.names=TRUE)

# Identify OTUs that are different just in Ghana lean people that will protect them against obesity





	