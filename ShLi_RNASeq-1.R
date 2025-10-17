###############################################
#    Shuang Liu Transcriptome Analysis        #
#           Date: July 17, 2021               #
###############################################

#Libraries
###############################################
library(DESeq2)
library(ggplot2)
library(stringr)
library(vegan)
library(RColorBrewer)

#Functions
###############################################
#create PCOA plot 
plot_pcoa <-function(pcoa, conditions, shapeby, title) {
  eig<- eigenvals(pcoa)
  prop<-eig/sum(eig)
  PCOA1 <- paste("PCoA1 ",100*round(prop[1],3),"%")
  PCOA2 <- paste("PCoA2 ",100*round(prop[2],3),"%")
  
  pcoa.sum <- summary(pcoa)
  pcoa.sum.sites  <- data.frame(pcoa.sum$sites[,1:2])       # PC1 and PC2
  pcoa.sum.species  <- data.frame(pcoa.sum$species[,1:2])     # loadings for PC1 and PC2
  
  if(!is.null(conditions)) {
    pcoa.sum.sites["colorby"] <- conditions
  }
  else { pcoa.sum.sites["colorby"] = rep("1", nrow(pcoa.sum.sites)) }
  
  if(!is.null(shapeby)) {
    pcoa.sum.sites["shapeby"] <- shapeby
  }
  else { pcoa.sum.sites["shapeby"] = rep("1", nrow(pcoa.sum.sites)) }
  
  
  pcoa_plot <- ggplot(pcoa.sum.sites, aes(x=MDS1, y=MDS2)) + 
    geom_jitter(aes(color = colorby, shape = shapeby), size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() +
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           text = element_text(size=14), 
           legend.background = element_blank(), 
           legend.title = element_blank(), 
           legend.key = element_blank()) +
    ggtitle(title) + 
    xlab (PCOA1) + 
    ylab (PCOA2) + 
    coord_fixed()
  
  return(pcoa_plot)
}

CompileResults <- function(res, padj=0.001, l2fc=5) {
  res.df <- data.frame(res)
  res.df <- res.df[order(res.df$log2FoldChange, decreasing = T), ]
  
  #Filter the significant DESeq2 Results
  res.df <- res.df[which(abs(res.df$log2FoldChange) > l2fc & res.df$padj < padj), ]
  res.df$Sequence <- factor(rownames(res.df), levels = rev(rownames(res.df)))
  return(res.df)
}

'%!in%' <- function(x,y)!('%in%'(x,y))


# Import Data
###############################################
setwd("T:/Bioinformatics/Analysis/ShLi")

#Import Metadata (prepared by Shuang Liu)
metadata <- read.delim("PS_gill_RNA_metafile.txt")
tail(metadata)

#import TPM data
TPM_quantmerge <- read.delim("TPM_quantmerge.filter.txt",
                             header=T)
rownames(TPM_quantmerge) <- TPM_quantmerge$Name; TPM_quantmerge$Name <- NULL

#Sort the dataset and keep the top 10,000 transcripts
TPM_quantmerge <- TPM_quantmerge[order(rowSums(TPM_quantmerge), decreasing = T), ]
TPM_quantmerge.top <- TPM_quantmerge[1:10000, ]

#Remove any transcripts with < 1 mean TPM across samples
TPM_quantmerge.keep <- TPM_quantmerge[which(rowMeans(TPM_quantmerge) > 1.0), ]


#Import Annotations
Annotations <- read.delim("annotation/est_transcripts95.annotation.tsv",
                          header = T)
Annotations$Transcript <- gsub("\\.p[1-9]+", "", Annotations$Peptide)
Annotations <- Annotations[!duplicated(Annotations$Transcript), ]
head(Annotations)


#Create Summary Plots
#################################################
#PCOA Plot
#################################################
pcoa <- capscale(t(TPM_quantmerge.top) ~ 1, dist="bray")
pcoa.plot <- plot_pcoa(pcoa, metadata$Region , metadata$treatment, "PCoA") + scale_color_manual(values = brewer.pal(9,"Set3"))
pcoa.plot


#HEATMAP
#################################################
TPM_quantmerge.top <- TPM_quantmerge.top[, match(metadata$Sample, colnames(TPM_quantmerge.top))]
colnames(TPM_quantmerge.top) == metadata$Sample

n <- nlevels(metadata$Region)
dat.col <- data.frame(Region=unique(metadata$Region),
                      RegionColors=brewer.pal(n,"Set3"))
SampleCol <- as.character(merge(metadata,dat.col)$RegionColors)

plot.new()
heatmap3::heatmap3(TPM_quantmerge.top[1:500,],
                   labRow=FALSE,
                   cexCol = 0.5,
                   showRowDendro = T,
                   showColDendro = T,
                   ColSideColors = SampleCol,
                   legendfun=function() heatmap3::showLegend(legend=as.character(dat.col$Region),
                                                   col=as.character(dat.col$RegionColors),
                                                   cex=1.5))


###############################################
#DeSeq2
###############################################

#import numReads data
quantmerge <- read.delim("est_numreads_quantmerge.txt",
                             header=T)
rownames(quantmerge) <- quantmerge$Name; quantmerge$Name <- NULL

#Match numeric data with filtered TPM dataset
quantmerge <- quantmerge[which(rownames(quantmerge) %in% rownames(TPM_quantmerge.keep)), match(metadata$Sample, colnames(quantmerge))]

#Convert to integer
quantmerge <- as.matrix(quantmerge); mode(quantmerge) <- "integer"
quantmerge[is.na(quantmerge)] <- 0


#Run DESeq2 Analysis
#############################################
quantmerge.filter <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.filter <- quantmerge.filter[which(rowSums(quantmerge.filter) != 0), ]
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]

dds <- DESeqDataSetFromMatrix(countData = quantmerge.filter, colData = metadata.filter, design = ~ treatment + Region) #
dds <- DESeq(dds)

#Extract Normalized count data
dds.normcounts <- counts(dds, normalized=T)
#sizeFactors(dds)

#Plot PCoA for deseq2 normalized count table (outliers removed)
pcoa <- capscale(t(dds.normcounts) ~ 1, dist="bray")
pcoa.plot <- plot_pcoa(pcoa, metadata.filter$Region , metadata.filter$treatment, "PCoA") + scale_color_manual(values = brewer.pal(9,"Set3"))
pcoa.plot + geom_text(label = metadata.filter$Sample, color = "#cecece", size = 2.4, nudge_x = 0.05, nudge_y = 0.05)

pcoa.model <-  adonis(t(dds.normcounts) ~ Region*treatment, strata=Population, data = metadata.filter)
pcoa.model

#Possible Outliers
#CR4SW
#LCR1S2

#Extract the DESeq2 Results and order
#############################################
#Extract Results
res.cR.cL <- results(dds, contrast=c("Region","coastR","coastL"))
res.iL.cL <- results(dds, contrast=c("Region","inlandL","coastL"))
res.cR.iL <- results(dds, contrast=c("Region","coastR","inlandL"))
res.SW.AW <- results(dds, contrast=c("treatment","SW","AW"))

#Compile Results
res.cR.cL.df <- CompileResults(res.cR.cL)
res.iL.cL.df <- CompileResults(res.iL.cL)
res.cR.iL.df <- CompileResults(res.cR.iL)
res.SW.AW.df <- CompileResults(res.SW.AW)

#Add annotations to dataframe cR.iL
Annotations.cR.iL <- Annotations[match(rownames(res.cR.iL.df), Annotations$Transcript), ]
res.cR.iL.df <- cbind(res.cR.iL.df, Annotations.cR.iL)
head(res.cR.iL.df, 50)

#Add annotations to dataframe SW.AW
Annotations.SW.AW <- Annotations[match(rownames(res.SW.AW.df), Annotations$Transcript), ]
res.SW.AW.df <- cbind(res.SW.AW.df, Annotations.SW.AW)
tail(res.SW.AW.df, 50)

#Isolate significant sequences
SigSequences <- c(rownames(res.cR.cL.df),
                  rownames(res.iL.cL.df),
                  rownames(res.cR.iL.df))
SigSequences <- SigSequences[!duplicated(SigSequences)]
length(SigSequences)

#Annotations for significant sequences (look for GO terms)
Sig.Annotations <- Annotations[which(Annotations$Transcript %in% SigSequences), ]
RefSeq <- gsub(" .*", "", Sig.Annotations$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations$PFAMID, 
                       PGAMGO = Sig.Annotations$PFAMGO), 
            "annotation/Sig.Annotations.tsv",
            sep="\t", row.names = F, quote = F)

#Heatmap of Significant Results for Region
#############################################
dds.normcounts.sig <- dds.normcounts[which(rownames(dds.normcounts) %in% SigSequences), ]
colnames(dds.normcounts.sig) == metadata.filter$Sample

n <- nlevels(metadata.filter$Region)
dat.col <- data.frame(Region=unique(metadata.filter$Region),
                      RegionColors=brewer.pal(n,"Set3"))
SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)

plot.new()
heatmap3::heatmap3(dds.normcounts.sig,
                   labRow=FALSE,
                   cexCol = 0.9,
                   showRowDendro = T,
                   showColDendro = T,
                   ColSideColors = SampleCol,
                   legendfun=function() heatmap3::showLegend(legend=as.character(dat.col$Region),
                                                             col=as.character(dat.col$RegionColors),
                                                             cex=1.5))
#Possible Outliers
#KL4AW
#FRP5SW


#Plot differential expression transcripts
ggplot(res.df, aes(y = Sequence, x = log2FoldChange)) +
  geom_bar(stat="identity", aes(fill = log2FoldChange)) +
  theme(axis.text.y = element_text(size = 5),
        legend.position = 'none')

#MA Plot
plotMA(res, ylim=c(-27,27))
