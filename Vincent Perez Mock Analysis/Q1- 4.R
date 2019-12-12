################This is an R Script for a mock analysis in R using RNA-seq data##########################
###By Vincent Perez
###12/06/2019
#===================================================================================================#
### In this R script, I will answer the 4 questions found in the README file (Q1- Q4). All of my markup will be indicated by 3 hashtags ("###") 
### and 1 hashtag ("#") indicates code that I want to keep but does not need to be ran for the analysis to work. There are 4 sections of code for each question. Each section is 
### And each section is broken up by "==". I hope it is clear, and thanks!

#BiocManager::install("ConsensusClusterPlus")
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(tidyr)
library(dplyr)
library(DESeq2)
library(stringr)
library(stats)
library(grid)
library(gridExtra)
setwd("C:/Users/vincent/Google Drive/yehlab/MockAnalysis-master/MockAnalysis-master/Vincent Perez Mock Analysis/")
load(file = "Q1-4.RData")

#===================================================================================================#
### Q1 - Perform a consensus clustering for 150 samples using genes in Moffitt-basal25.txt (orange) 
### and Moffitt-classical25.txt (blue). Use K=2 for both column and row. Plot a heatmap to show the 
### result. The resultant two clusters of samples are basal-like (orange) and classical (blue) subtypes respectively. 

data_basal<-read.csv("data_Moffitt_basal25.csv")                                    ### I have cleaned up the ex.csv and featinfo.csv in excel to now be 2 separate files corresponding to Moffitt_basal25.txt and Moffitt_classical25.txt genes
data_classical<-read.csv("data_Moffitt_classical25.csv")                            ### These 2 new datasets will only have expression values for the 25 genes in each file of the 150 paties. Here I will load them in and view them.
row.names(data_basal)<-data_basal[,1]                                               ### Note that there were 2 genes in the classical dataset not found in the RNA-sequencing datasets, so there are only 23 genes in "data_classical"
data_basal<-data_basal[,2:151]                                                      ### Remove the redundant first column
data_basal<-na.omit(data_basal)
data_basal<-data.matrix(data_basal)                                                 ### Convert to a matrix
row.names(data_classical)<-data_classical[,1]                                       ### Do the same for the Moffitt_classical25 file
data_classical<-data_classical[,2:151]
data_classical<-na.omit(data_classical)
data_classical<-data.matrix(data_classical)
d<-rbind(data_basal,data_classical)                                                 ### For an additional downstream analysis, I am combining these to with rbind
db<-data_basal                                                                     
dc<-data_classical                                                                  ### Now I have 3 datasets to work with (d, db, and dc)
mads<-apply(d,1,mad)                                                                ### Begin consensus cluster analysis on the 3 datasets (d = data_basal and data_classical combing, db = data_basal, dc = data_classical)
d<-d[rev(order(mads))[1:5000],]                                                     ### Here we are selecting the top 5000 most variable genes measured by median absolute deviation (mads)
d<-na.omit(d)
mads<-apply(db,1,mad)
db<-db[rev(order(mads))[1:5000],]                                                   ### Do for all 3 datasets
db<-na.omit(db)
mads<-apply(dc,1,mad)
dc<-dc[rev(order(mads))[1:5000],]
dc<-na.omit(dc)
d <- sweep(d,1, apply(d,1,median,na.rm=T))                                          ### clean the dataset by removing NAs
db <- sweep(db,1, apply(db,1,median,na.rm=T))
dc <- sweep(dc,1, apply(dc,1,median,na.rm=T))

subfolder_names <- c("Combined","Basal","Classical")                                ### Write a loop to create 3 new directories for the results
for (j in 1:length(subfolder_names)){
  folder<-dir.create(paste0("Consensus Results ",subfolder_names[j]))
}

d_results <- ConsensusClusterPlus(d=d,maxK=6,reps=150,pItem=0.8,pFeature=1,         ### Create plots for each comparison and write to their directories
                                  title="Consensus Results Combined", 
                                  clusterAlg="hc",
                                  distance="pearson",
                                  seed=1262118388.71279,                            
                                  tmyPal=brewer.pal(11, "RdYlBu"),
                                  plot="png")
db_results <- ConsensusClusterPlus(d=db,maxK=6,reps=150,pItem=0.8,pFeature=1,
                                  title="Consensus Results Basal",
                                  clusterAlg="hc",
                                  distance="pearson",
                                  seed=1262118388.71279,
                                  tmyPal=brewer.pal(11, "RdYlBu"),
                                  plot="png")
dc_results <- ConsensusClusterPlus(d=dc,maxK=6,reps=150,pItem=0.8,pFeature=1,
                                  title="Consensus Results Classical",
                                  clusterAlg="hc",
                                  distance="pearson",
                                  seed=1262118388.71279,
                                  tmyPal=brewer.pal(11, "RdYlBu"),
                                  plot="png")

#d_results[[2]][["consensusMatrix"]]                                                ### Additional results can be viewed with the code on the left
#d_results[[2]][["consensusTree"]]
#icl <- calcICL(d_results,title="C:/Users/vincent/Desktop/results/dbdc/",plot="png")
#icl[["clusterConsensus"]]
#icl[["itemConsensus"]][1:5,]

#===================================================================================================#
### Q3 - Identify all the differentially expressed (DE) genes between basal-like and classical samples 
### using DESeq2 (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). 
### Plot MA plot and volcano plot to show the results. Cut-offs to use: p=0.05, log2(fold-change)=1.

                                                                                  ### Here we are prepping 2 objects for DESeq2 Input. 1) the coldata data.frame and 2) the counts (cts) matrix
coldata<-TCGA_PAAD_plus$sampInfo                                                  ### load the sample info
coldata<-coldata%>%drop_na(Grade)                                                 ### Clean up the experimental information by removing "NAs"
coldata<-coldata[order(coldata$Tumor.Sample.ID),]                                 ### now order them alphabetically
row.names(coldata)<-coldata[,1]                             
coldata<-coldata[,c(-3:-48,-52:-62,-64:-73)]                                      ### subset the data
coldata$SSC.subtype<-as.factor(coldata$SSC.subtype)                               ### convert subtype to factor
cts<-read.csv("raw_counts.csv", check.names = F ,header=T)                        ### create the counts file
cts<-cts[!duplicated(cts[,1]),]                                                   ### Remove duplicates
row.names(cts)<-cts[,1]
cts<-cts[,colnames(cts) %in% coldata[,1]]                                         ### clean up cts by indexing with coldata
cts<-cts[,order(colnames(cts))]  
cts<-round(cts,0)                                                                 ### The counts data needs to be in integer form for DESeq2 Analysis. This is a minor alteration that should not affect the final analysis.
dds <- DESeqDataSetFromMatrix(countData = cts,                                    ### Create the Deseq2 model
                              colData = coldata,
                              design= ~ Batch + SSC.subtype)
dds$treatment <- relevel(dds$SSC.subtype, "Classical" )                           ### build the statistical model to make sure that control (classical) is the first level in the treatment factor
                                                                                  ### this will ensure log2fc is made relative to control not other way around
dds$Gender <- droplevels( dds$Gender )                                            ### drop the additional uneccesary factors 
dds$Ethnicity <- droplevels( dds$Ethnicity )
dds$Race <- droplevels( dds$Race )
as.data.frame(colData(dds))                                                       ### Double check we have the right samples
dds <- DESeq(dds)                                                                 ### Run the deg pipeline
res<-results(dds)                                                                 ### Check the results
#res2<-results(dds, contrast = c("treatment", "Basal", "Classical"))              ### Here we can compare any 2 levels of a variable with contrast
mcols(res,use.names=TRUE)                                                         ### Take a look at the analysis summary
#sum( res$pvalue < 0.005, na.rm=TRUE )                                            ### Check how many genes have p<0.01
#sum( res$padj < 0.05, na.rm=TRUE )                                               ### Check the BH-adjusted pvalues
#resSig <- res[ which(res$padj < 0.1 ), ]  
#head( resSig[ order( resSig$log2FoldChange ), ] )                                ### subset with log2 fc
#tail( resSig[ order( resSig$log2FoldChange ), ] )
png("MA Plot.png", height = 5, width = 5, units = 'in',  res=300)                 ### Here we will make the MA plot and write it to a file (p-value =<0.05 and |logFc2|>= 1)
plotMA(res, alpha = 0.05, 
       main = "MA Plot of Classical vs Basal Subtypes", 
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change",
       ylim=c(-5,5))
dev.off()
png("Volcano Plot.png", height = 5, width = 5, units = 'in',  res=300)
alpha <- 0.05                                                                     ### Here we will make the volcano plot (p-value =<0.05 and |logFc2|>= 1)
cols <- densCols(res$log2FoldChange, -log10(res$pvalue), nbin=1000)                
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="DEGs in Basal Cells Vs Classical", 
     xlab="Log2Fold Change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.5)
abline(v=c(-1,1), col="brown")                                                    ### Add p-value and log2 FC cut-offs
abline(h=-log10(alpha), col="brown")
gn.selected <- abs(res$log2FoldChange) > 4 & res$padj < 0.0005                    ### Label genes that hit a desired threshold
text(res$log2FoldChange[gn.selected],                                             ### Add the text to the graph
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.5)
dev.off()
#===================================================================================================#
### Q2 - Test if each of the genes in AURKA_sig_genes.txt is differentially expressed between basal-like 
### and classical subtypes. Use boxplot to illustrate the results with p-values labeled on the figures.

#a_genes<-read.csv("AURKA_sig_genes.csv",header=T)                                ### Define the genes of interest
a_genes=c("AURKA",
          "MYCN",
          "NFKBIA",
          "RALA",
          "SRC",
          "CHUK",
          "STAT3",
          "JAK2",
          "BRCA1",
          "BRCA2",
          "RELA",
          "TP53",
          "TP73",
          "GSK3B")                                                         
a_genes %in% names(dds)                                                           ### Check and make sure they are all in your dds object
tcounts <- t(log2((counts(dds[a_genes, ], normalized=TRUE, 
                          replaced=FALSE)+.5))) %>%
merge(colData(dds), ., by="row.names") %>%                                        ### Build a new matrix that can be graphed into a boxplot
gather(gene, expression, (ncol(.)-length(a_genes)+1):ncol(.))

tcounts %>%                                                                       ### Check your new matrix
  select(Row.names, gene, expression, SSC.subtype) %>% 
  head %>% 
  knitr::kable()

png("Boxplots.png", height = 6, width = 10, units = 'in', res=300)
ggplot(tcounts, aes(SSC.subtype, expression, fill=SSC.subtype)) +                 ### Make and save the boxplots for each gene
  geom_boxplot(outlier.size = 0.5) + 
  stat_compare_means(size = 2, label="p.format", label.x.npc = "center")  +
  facet_wrap(~gene, scales="free_y",shrink=TRUE) + 
  geom_point(position=position_dodge(width=0.1),size = 0.5) +
  labs(x="Tissue Subtype", 
       y="Expression (Log2 Normalized Counts)", 
       fill="Tissue Subtype", 
       title="Boxplots of AURKA_sig_genes.txt") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold")) +
  scale_fill_brewer(palette="Set1") +
  theme_classic()
dev.off()

#===================================================================================================#
### Q4 - Is subtype assocaited with the clinical variables (sex, age, gender and race)? Use a table with
### p-values to show the results.

summary(coldata)                                                                  ### First, we take a look at the categorical variables
grouped_gender<-group_by(coldata, Gender, SSC.subtype)                            ### Then we need to subset our categorical variables
sum_gender<-summarise(grouped_gender, count=n())                                  ### Count how many females and males have the two subtypes
tbl <- matrix(data=c(14, 55, 19, 62), nrow=2, ncol=2, byrow=T)                    ### Create a 2X2 table of these sub-totals shown in "summary"
dimnames(tbl) = list(Gender=c('F', 'M'), Subtype=c('B', 'C'))                     ### Add labels
chi2_gender <- chisq.test(tbl, correct=F)                                         ### Because these are categorical, we can use a chi-squared test to test if there is correlation between the categories
chi2_gender                                                                       ### Our test shows a x-squared of 0.2177 and p-value of 0.6407
V<-sqrt(chi2$statistic/sum(tbl))                                                  ### Our V is 0.0381 (smaller the v means low correlation)
V
coldata2<-coldata%>%drop_na(Ethnicity)
grouped_ethnicity<-group_by(coldata2, Ethnicity, SSC.subtype)                     ### We will do the same for ethnicity, but first must drop NAs
sum_ethnicity<-summarise(grouped_ethnicity, count=n())                            ### Count how many latinos/non latinos have the two subtypes
tbl2 <- matrix(data=c(1, 4, 25, 82), nrow=2, ncol=2, byrow=T)                     ### Create a 2X2 table of these sub-totals shown in "summary"
dimnames(tbl2) = list(Ethnicity=c('hispanic or latino',
                              'not hispanic or latino'), Subtype=c('B', 'C'))     ### Add labels
chi2_ethnicity <- chisq.test(tbl2, correct=F)                                     ### Again chi-squared test to test if there is correlation between the categories
chi2_ethnicity                                                                    ### Our test shows a x-squared of 0.30335 and p-value of 0.8617
V<-sqrt(chi2_ethnicity$statistic/sum(tbl2))                                       ### Our V is 0.04409536 (smaller the v means low correlation)
V
coldata3<-coldata%>%drop_na(Race)
grouped_race<-group_by(coldata3, Race, SSC.subtype)                               ### We will do the same for ethnicity, but first must drop NAs
sum_race<-summarise(grouped_race, count=n())                                      ### Count how many latinos/non latinos have the two subtypes
tbl3 <- matrix(data=c(6, 8, 31,101), nrow=2, ncol=2, byrow=T)                     ### Create a 2X2 table of these sub-totals shown in "summary"
dimnames(tbl3) = list(Race=c('Other (Asian & Black)', 'White'),                   ### Add labels 
                      Subtype=c('B', 'C'))     
chi2_race <- chisq.test(tbl3, correct=F)                                          ### Again chi-squared test to test if there is correlation between the categories
chi2_race                                                                         ### Our test shows a x-squared of 2.5107 and p-value of 0.1131
V<-sqrt(chi2_race$statistic/sum(tbl2))                                            ### Our V is 0.1497219 (smaller the v means low correlation)
V
corr_tbl=matrix(data=c(0.2177,0.3033,2.511,                                       ### Now here we will build a table to summarize our results
                       0.6407,0.8617,0.1131,
                       0.0381,0.04410,0.1497), nrow=3, ncol=3, byrow=T)           
dimnames(corr_tbl)=list(Chi_Squared_Values=c("Correlation (X^2)","P-Value","V"),  ### Label the x and y
                        Demographic=c("Gender","Ethnicity","Race"))

tt <- ttheme_minimal()                                                            ### Choose a theme and plot a table
png("Correlation Table.png", height = 3, width = 5, units = 'in', res=300)
grid.table(corr_tbl, theme=tt)
dev.off()

#===================================================================================================#
### Lastly, save the work to an .RData file.

save.image("C:/Users/vincent/Google Drive/yehlab/MockAnalysis-master/MockAnalysis-master/Vincent Perez Mock Analysis/Q1-4.RData")




