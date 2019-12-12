################This is an R Script to process RNAseq Data (quant files) from salmon##########################
#Vincent Perez
#01/23/2019
#===================================================================================================#
###SECTION 1 - CREATING THE TX2GENE FILE FOR TXIMPORT (PLEASE NOTE, THIS FIRST SCRIPT WILL NOT WORK UNLESS YOU RUN IT
### IN THE SAME DIRECTORY AS THE .GTF FILE MENTIONED BELOW)

library(tidyverse)
library(GenomicFeatures)
library(tximport)
library(readr)
library(edgeR)
library(DOSE)
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)
library(ReactomePA)
setwd("Z:/Vincent Perez/RNA Seq Data/Exp 134/")


gtf <- file.path(dir2,"Mus_musculus.GRCm38.95.gtf")                               ### Set up filenames. Only change these lines to point to the GTF file.                               
txdb.filename <- gsub("gtf", "sqlite", gtf)                                       ### Lets generate some file names. The "gsub" command replaces patterns with other patterns.
tx2gene.filename <- gsub("gtf", "tx2gene.csv", gtf)
if (!file.exists(txdb.filename)) {                                                ### Check to see if the sqlite file exists in dir2.
  message(paste(txdb.filename, "doesn't exist. Creating and saving to disk..."))  ### By default we assume the sqlite file doesn't exist.
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  message("Done.")
} else {
  message(paste(txdb.filename, "found. Loading..."))                              ### However if we actually do have that, then we'll load it.
  txdb <- loadDb(txdb.filename)
  message("Done.")
}

keytypes(txdb)
tx2gene <- mapIds(txdb,
                  keys=keys(txdb, "GENEID"),                                      ### Create the tx2gene
                  column="TXNAME",
                  keytype="GENEID",
                  multiVals="list") %>%
  enframe(name="ensgene", value="enstxp") %>%
  unnest() %>%
  dplyr::select(enstxp, ensgene) %>%
  distinct()
write_csv(tx2gene, tx2gene.filename)                                              ### Write out tx2gene with the header

#===================================================================================================#
###SECTION 2 - GENERATING COUNT FILES ASSOCIATED WITH GENE NAMES WITH TXIMPORT

dir<-("/Z:/Vincent Perez/RNA Seq Data/Exp 134")
samples <- read.csv(file.path(dir, "samples_info.txt"), header = TRUE)
samples
files <- file.path(dir, "quants", samples$Salmon_out, "quant.sf")                 ### Tell R where each and every 'quant.sf' file is located, then check and see if those files exist
names(files) <- paste0("sample", 1:16)
all(file.exists(files))
tx2gene<-read_tsv("./Mus_musculus.GRCm38.95.tx2gene.txt")                         ### Importing our annotation file, tx ID, then gene ID 
head(tx2gene)                                                                     ### Note: this file was generated in the above section
quant_ex<-read_tsv('./quants/black629_quant/quant.sf')                            ### (Note necessary step) Count the number of IDs that are matching between your txToogene file and your quant files
results<-intersect(quant_ex$Name, tx2gene$Name)
length(quant_ex$Name)
length(results)
txi = tximport(files, type = "salmon",                                            ### Importing our annotation file, tx ID, then gene ID....generating counts (All arugments shown below)
               txIn = T, txOut = F, countsFromAbundance = "scaledTPM", 
               tx2gene = tx2gene, varReduce = FALSE,
               dropInfReps = FALSE, ignoreTxVersion = TRUE, 
               geneIdCol = tx2gene[,2], txIdCol = tx2gene[,1],
               abundanceCol = "TPM", countsCol = "counts", 
               lengthCol = "Length", importer = NULL)
names(txi)                                                                        ### txi=tximport(files,type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE)
head(txi$counts)
counts = txi$counts
write.csv(counts, "./counts.csv")
gene_length = txi$length

#===================================================================================================#
### SECTION 3 - GETTING DIFFERENTIALLY EXPRESSED GENES WITH edgeR 

setwd("C:/Users/vincent/Desktop/jmp_genomics_workspace/RNA-seq 135 136 Workspace/counts files/")
counts = read.csv("counts.csv", header=TRUE)                                      ### Make the counts data.frame we'll be using
head(counts)
tail(counts)
group<-factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8))
dgList<-DGEList(counts=counts[,2:33], genes=counts$Gene.name, group=group)
dgList$samples
head(dgList$counts)
head(dgList$genes) 
countsPerMillion <- cpm(dgList)                                                  ### Here I'm filtering out transcripts based on their CPM (Counts per Million). I use CPM to account for..
summary(countsPerMillion)                                                        ### ..library size differences. I also am selecting for transcripts that are found in all samples, rather..
countCheck <- countsPerMillion > 1                                               ### ..than transcripts found in, say, 3/8 of the samples.
head(countCheck)
keep <- rowSums(cpm(dgList)>1)>=3
dgList <- dgList[keep,,keep.lib.sizes=FALSE]
summary(cpm(dgList)) 
dgList <- calcNormFactors(dgList, method="TMM")                                 ### NORMALIZATION
dgList$samples                                                                  ### "calcNormFactors" function normalizes for RNA composition by finding a set of scaling factors 
points<-c(0,1,15,19,2,6,17,25)                                                  ### Trimmed mean of M-values (TMM) is the default method for calcNormFactors, and it is the recommended
colors<-rep(c('black', 'red'), 4)
plotMDS(dgList, top=20000, col=colors[group], pch=points[group], cex = 1.0)     ### Create a PCA plot
legend('topright', legend=c( 'Group 1',
                             'Group 2', 
                             "Group 3",
                             "Group 4",
                             'Group 5',
                             'Group 6', 
                             "Group 7",
                             "Group 8"), 
       pch=points, col=colors,ncol=1)
designMat <- model.matrix(0~group)                                             ### Lets setup the model, we're using GLM for MULTIPLE FACTORS, if single factor remove GLM in lines
designMat                                                                       ### GLM dispersion estimate
y<-estimateDisp(dgList, designMat)
y<-estimateTagwiseDisp(y)
y$common.dispersion
plotBCV(y)
fit<-glmQLFit(y, designMat, robust=TRUE)                                        ### Perform QL or F-test. F-test is robust and useful when replicate numbers are small. QL test is
plotQLDisp(fit)                                                                 ### These tests are only useful in experiments with MULTIPLE FACTORS! "coef" tells the model
qlf.2vs1<-glmQLFTest(fit, coef=c(2))                                            ### What to compare. "coef=2" will be 2vs1, whereas, "contrast=c(0,-1,1) would ve 3 vs 2. 
topTags(qlf.2vs1)                                                               ### If you're looking for genes that are different between any of the 4 groups it'd be "coef=3:2"
edgeR_results <- topTags(qlf.2vs1, n=Inf)
deGenes <- decideTestsDGE(qlf.2vs1, p=0.005)
deGenes <- rownames(qlf.2vs1)[as.logical(deGenes)]
plotSmear(qlf.2vs1, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
write.csv(edgeR_results, "./edgeR_results Chow Females.csv")                    ### Now write the results of the DEGs to a file for downstream analysis outside of R
DEGs<-edgeR_results$table

#===================================================================================================#
### SECTION 4 - PERFORMING PATHWAY ENRICHMENT ANALYSIS WITH clusterProfiler

rm(list = ls())                                                                 ### Clear your workspace
setwd("C:/pathway/to/the/edgeR/results/csv/")
df1<-read.csv("clusterProfiler_input.csv", header = T)
head(df1)
formula_res.twogroups <- compareCluster(Entrez~Treatment+Group,                 ### Formula interface of compareCluster
                                        data=df1, fun="enrichPathway", 
                                        organism="mouse")
xx=summary(as.data.frame(formula_res.twogroups))                                ### Summarize the data
xx
p=dotplot(formula_res.twogroups, x=~Treatment, showCategory = 15) +             ### Make a plot
  ggplot2::facet_grid(~Group)
p+theme_pubr(base_size = 10, x.text.angle = 45, legend = "right")
p+labs_pubr()

#===================================================================================================#
save.image("C:/Users/vincent/Google Drive/yehlab/MockAnalysis-master/MockAnalysis-master/Vincent Perez Mock Analysis/Q5_DEG_from_raw_files.RData")
