#This is a function designed to generate plots facilitated the visualisation of large datasets.
#Copyright (C) 2017  Aurelien Dugourd

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(gskb)
library(piano)
library(BioNet)
library(igraph)
library(reshape)
library(data.table)
library(pheatmap)
library(limma)
library(grid)
library(gridExtra)
library(GSEABase)

data("mm_metabolic")

unfactor <- function(df){ 
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df
}

TPO_vs_control_diffExpression <- read.csv2("~/Documents/Kramann/M2/TPO_vs_control_diffExpression.csv")

gene_to_term <- data.frame(NA,NA)
names(gene_to_term) <- c("gene","term")

i <- 0
for (pathway in mm_metabolic)
{
  pathway <- unlist(pathway)
  pathway <- as.data.frame(pathway)
  pathway$term <-  pathway[1,]
  pathway <- pathway[c(-1,-2),]
  names(pathway) <- c("gene","term")
  gene_to_term <- rbind(gene_to_term,pathway)
}

gene_to_term_no_reac <- gene_to_term[grep("REACTOME",gene_to_term[,"term"], invert = TRUE),]
geneSet <- loadGSC(gene_to_term_no_reac)

myFC <- TPO_vs_control_diffExpression$log2FoldChange
names(myFC) <- toupper(TPO_vs_control_diffExpression$X)

myPval <- TPO_vs_control_diffExpression$padj
names(myPval) <- toupper(TPO_vs_control_diffExpression$X)

myTval <- TPO_vs_control_diffExpression$stat
names(myTval) <- toupper(TPO_vs_control_diffExpression$X)

###Run the GSA
gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean")
gsaRes2 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "median")
gsaRes3 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "sum")
gsaRes4 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "maxmean")
gsaRes5 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fisher")
gsaRes6 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "stouffer")
gsaRes7 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "tailStrength")
gsaRes8 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "wilcoxon")
gsaRes9 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "page")


resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes8,gsaRes9)
names(resList) <- c("mean","median","sum","maxmean","fisher", "stouffer","tailStrength","wilcoxon","page")

ch <- consensusHeatmap(resList,cutoff=50,method="median", ncharLabel = 100, cellnote = "medianPvalue", cex = 0.2) ##The results are strange

consensus <- ch$pMat

write.csv(consensus,"~/Documents/Kramann/M2/consensus_mm_metabolic.csv")
