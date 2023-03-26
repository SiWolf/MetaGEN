# --------------------------------------------------------------------------------------------------------
# Title: pathway_analysis.R
# Author: Silver A. Wolf
# Last Modified: Sun, 26.03.2023
# Version: 0.0.1
# --------------------------------------------------------------------------------------------------------

# Libraries
library("gtools")
library("openxlsx")
library("splitstackshape")
library("tidyverse")

# Read metadata
meta.raw <- read.xlsx("../metadata/23_03_Horses_Overview.xlsx", sheet = 1)
meta.raw <- meta.raw[,-c(1)]
colnames(meta.raw)[1] <- "SampleID"

# Read pathway abundancies
c1 <- read.csv("../output/03_functional_analysis/humann3/D61001/humann_D61001_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61001 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61002/humann_D61002_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61002 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61003/humann_D61003_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61003 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61004/humann_D61004_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61004 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61005/humann_D61005_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61005 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61006/humann_D61006_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61006 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61007/humann_D61007_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61007 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61008/humann_D61008_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61008 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61009/humann_D61009_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61009 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61010/humann_D61010_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61010 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61011/humann_D61011_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61011 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61012/humann_D61012_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61012 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61013/humann_D61013_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61013 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61014/humann_D61014_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61014 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61015/humann_D61015_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
D61015 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25286/humann_H25286_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25286 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25287/humann_H25287_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25287 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25288/humann_H25288_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25288 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25289/humann_H25289_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25289 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25290/humann_H25290_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25290 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25294/humann_H25294_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25294 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25296/humann_H25296_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25296 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25297/humann_H25297_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
H25297 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12271/humann_I12271_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12271 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12272/humann_I12272_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12272 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12273/humann_I12273_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12273 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12276/humann_I12276_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12276 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12277/humann_I12277_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12277 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12278/humann_I12278_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12278 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12279/humann_I12279_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12279 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12280/humann_I12280_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12280 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12281/humann_I12281_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12281 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12282/humann_I12282_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
I12282 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10519/humann_J10519_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10519 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10520/humann_J10520_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10520 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10521/humann_J10521_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10521 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10522/humann_J10522_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10522 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10523/humann_J10523_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10523 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10524/humann_J10524_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10524 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10525/humann_J10525_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10525 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10526/humann_J10526_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10526 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10527/humann_J10527_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10527 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10528/humann_J10528_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10528 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10529/humann_J10529_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10529 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10530/humann_J10530_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
J10530 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505693/humann_SRR10505693_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
SRR10505693 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505694/humann_SRR10505694_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
SRR10505694 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505695/humann_SRR10505695_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
SRR10505695 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505696/humann_SRR10505696_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
SRR10505696 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505697/humann_SRR10505697_pathabundance.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Pathway", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Pathway_2),]
SRR10505697 <- data.frame(PW = c3[,2], RPK = c3[,1])

pathways.raw <- list(D61001, D61002, D61003, D61004, D61005, D61006, D61007, D61008, D61009, D61010,
                     D61011, D61012, D61013, D61014, D61015, H25286, H25287, H25288, H25289, H25290,
                     H25294, H25296, H25297, I12271, I12272, I12273, I12276, I12277, I12278, I12279,
                     I12280, I12281, I12282, J10519, J10520, J10521, J10522, J10523, J10524, J10525,
                     J10526, J10527, J10528, J10529, J10530, SRR10505693, SRR10505694, SRR10505695,
                     SRR10505696, SRR10505697) %>% reduce(full_join, by = "X..Pathway_1")

pathways.raw[is.na(pathways.raw)] <- 0
write_tsv(pathways.raw, file = "../output/08_visualization/tab_pathway_abundances.csv")

# Significantly different pathways between SSG and 5DG at t0
pw = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t0", ]$SampleID, "_Abundance", sep = ""))
fdg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t0", ]$SampleID, "_Abundance", sep = ""))

ssg.pathways <- pathways.raw[, colnames(pathways.raw) %in% ssg.samples]
fdg.pathways <- pathways.raw[, colnames(pathways.raw) %in% fdg.samples]

for (row in pathways.raw$X..Pathway_1){
  ssg.abundance = ssg.pathways[ssg.pathways$X..Pathway_1 == row, ]
  fdg.abundance = fdg.pathways[fdg.pathways$X..Pathway_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  pw[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Pathway = pw, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_pathways_diff_t0.tsv")

# Significantly different pathways between SSG and 5DG at t1
pw = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t1", ]$SampleID, "_Abundance", sep = ""))
fdg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t1", ]$SampleID, "_Abundance", sep = ""))

ssg.pathways <- pathways.raw[, colnames(pathways.raw) %in% ssg.samples]
fdg.pathways <- pathways.raw[, colnames(pathways.raw) %in% fdg.samples]

for (row in pathways.raw$X..Pathway_1){
  ssg.abundance = ssg.pathways[ssg.pathways$X..Pathway_1 == row, ]
  fdg.abundance = fdg.pathways[fdg.pathways$X..Pathway_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  pw[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Pathway = pw, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_pathways_diff_t1.tsv")

# Significantly different pathways between SSG and 5DG at t2
pw = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t2", ]$SampleID, "_Abundance", sep = ""))
fdg.samples <- c("X..Pathway_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t2", ]$SampleID, "_Abundance", sep = ""))

ssg.pathways <- pathways.raw[, colnames(pathways.raw) %in% ssg.samples]
fdg.pathways <- pathways.raw[, colnames(pathways.raw) %in% fdg.samples]

for (row in pathways.raw$X..Pathway_1){
  ssg.abundance = ssg.pathways[ssg.pathways$X..Pathway_1 == row, ]
  fdg.abundance = fdg.pathways[fdg.pathways$X..Pathway_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  pw[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Pathway = pw, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_pathways_diff_t2.tsv")

# Read gene abundancies
c1 <- read.csv("../output/03_functional_analysis/humann3/D61001/humann_D61001_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61001 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61002/humann_D61002_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61002 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61003/humann_D61003_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61003 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61004/humann_D61004_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61004 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61005/humann_D61005_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61005 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61006/humann_D61006_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61006 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61007/humann_D61007_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61007 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61008/humann_D61008_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61008 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61009/humann_D61009_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61009 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61010/humann_D61010_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61010 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61011/humann_D61011_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61011 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61012/humann_D61012_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61012 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61013/humann_D61013_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61013 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61014/humann_D61014_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61014 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/D61015/humann_D61015_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
D61015 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25286/humann_H25286_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25286 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25287/humann_H25287_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25287 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25288/humann_H25288_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25288 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25289/humann_H25289_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25289 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25290/humann_H25290_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25290 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25294/humann_H25294_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25294 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25296/humann_H25296_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25296 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/H25297/humann_H25297_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
H25297 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12271/humann_I12271_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12271 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12272/humann_I12272_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12272 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12273/humann_I12273_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12273 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12276/humann_I12276_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12276 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12277/humann_I12277_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12277 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12278/humann_I12278_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12278 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12279/humann_I12279_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12279 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12280/humann_I12280_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12280 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12281/humann_I12281_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12281 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/I12282/humann_I12282_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
I12282 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10519/humann_J10519_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10519 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10520/humann_J10520_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10520 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10521/humann_J10521_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10521 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10522/humann_J10522_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10522 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10523/humann_J10523_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10523 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10524/humann_J10524_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10524 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10525/humann_J10525_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10525 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10526/humann_J10526_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10526 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10527/humann_J10527_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10527 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10528/humann_J10528_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10528 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10529/humann_J10529_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10529 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/J10530/humann_J10530_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
J10530 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505693/humann_SRR10505693_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
SRR10505693 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505694/humann_SRR10505694_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
SRR10505694 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505695/humann_SRR10505695_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
SRR10505695 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505696/humann_SRR10505696_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
SRR10505696 <- data.frame(PW = c3[,2], RPK = c3[,1])

c1 <- read.csv("../output/03_functional_analysis/humann3/SRR10505697/humann_SRR10505697_genefamilies.tsv", sep = "\t")
c2 <- concat.split(data = c1, split.col = "X..Gene.Family", sep = "\\|", drop = TRUE)
c3 <- c2[is.na(c2$X..Gene.Family_2),]
SRR10505697 <- data.frame(PW = c3[,2], RPK = c3[,1])

genes.raw <- list(D61001, D61002, D61003, D61004, D61005, D61006, D61007, D61008, D61009, D61010,
                     D61011, D61012, D61013, D61014, D61015, H25286, H25287, H25288, H25289, H25290,
                     H25294, H25296, H25297, I12271, I12272, I12273, I12276, I12277, I12278, I12279,
                     I12280, I12281, I12282, J10519, J10520, J10521, J10522, J10523, J10524, J10525,
                     J10526, J10527, J10528, J10529, J10530, SRR10505693, SRR10505694, SRR10505695,
                     SRR10505696, SRR10505697) %>% reduce(full_join, by = "X..Gene.Family_1")

genes.raw[is.na(genes.raw)] <- 0
write_tsv(genes.raw, file = "../output/08_visualization/tab_gene_abundances.csv")

# Significantly different gene families between SSG and 5DG at t0
gf = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t0", ]$SampleID, "_Abundance.RPKs", sep = ""))
fdg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t0", ]$SampleID, "_Abundance.RPKs", sep = ""))

ssg.genes <- genes.raw[, colnames(genes.raw) %in% ssg.samples]
fdg.genes <- genes.raw[, colnames(genes.raw) %in% fdg.samples]

for (row in genes.raw$X..Gene.Family_1){
  ssg.abundance = ssg.genes[ssg.genes$X..Gene.Family_1 == row, ]
  fdg.abundance = fdg.genes[fdg.genes$X..Gene.Family_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  gf[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Gene_Family = gf, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_genes_diff_t0.tsv")

# Significantly different gene families between SSG and 5DG at t1
gf = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t1", ]$SampleID, "_Abundance.RPKs", sep = ""))
fdg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t1", ]$SampleID, "_Abundance.RPKs", sep = ""))

ssg.genes <- genes.raw[, colnames(genes.raw) %in% ssg.samples]
fdg.genes <- genes.raw[, colnames(genes.raw) %in% fdg.samples]

for (row in genes.raw$X..Gene.Family_1){
  ssg.abundance = ssg.genes[ssg.genes$X..Gene.Family_1 == row, ]
  fdg.abundance = fdg.genes[fdg.genes$X..Gene.Family_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  pw[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Gene_Family = gf, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_genes_diff_t1.tsv")

# Significantly different gene families between SSG and 5DG at t2
gf = c()
sr = c()
fr = c()
fc = c()
pv = c()
i = 1

ssg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "SSG" & meta.raw$Timepoint == "t2", ]$SampleID, "_Abundance.RPKs", sep = ""))
fdg.samples <- c("X..Gene.Family_1", paste("humann_", meta.raw[meta.raw$AB_Group == "5DG" & meta.raw$Timepoint == "t2", ]$SampleID, "_Abundance.RPKs", sep = ""))

ssg.genes <- genes.raw[, colnames(genes.raw) %in% ssg.samples]
fdg.genes <- genes.raw[, colnames(genes.raw) %in% fdg.samples]

for (row in genes.raw$X..Gene.Family_1){
  ssg.abundance = ssg.genes[ssg.genes$X..Gene.Family_1 == row, ]
  fdg.abundance = fdg.genes[fdg.genes$X..Gene.Family_1 == row, ]
  
  stat.df <- data.frame(GROUP = c(rep("SSG", length(ssg.samples)-1), rep("5DG", length(fdg.samples)-1)), ABN = c(as.character(ssg.abundance[-1]), as.character(fdg.abundance[-1])))
  stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABN), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE, exact = FALSE)
  
  pw[i] <- row
  sr[i] <- sum(ssg.abundance[-1])/(length(ssg.samples)-1)
  fr[i] <- sum(fdg.abundance[-1])/(length(fdg.samples)-1)
  fc[i] <- foldchange(fr[i], sr[i])
  pv[i] <- as.numeric(stat.res$p.value)
  
  i = i + 1
}

pathways.comparison = data.frame(Gene_Family = gf, SSG_Avg_RPK = sr, FDG_Avg_RPK = fr, FoldChange = fc, pval = pv)
write_tsv(pathways.comparison, file = "../output/08_visualization/tab_genes_diff_t2.tsv")