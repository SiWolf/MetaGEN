# --------------------------------------------------------------------------------------------------------
# Title: MetaGEN.R
# Author: Silver A. Wolf
# Last Modified: Wed, 26.01.2022
# Version: 0.4.7
# --------------------------------------------------------------------------------------------------------

# Libraries

library("circlize")
library("ComplexHeatmap")
library("dplyr")
library("edgeR")
library("ggpubr")
library("metagMisc")
library("microbiome")
library("microbiomeExplorer")
library("randomcoloR")
library("PathoStat")
library("phyloseq")
library("stringr")
library("taxonomizr")
library("tidyr")

# --------------------------------------------------------------------------------------------------------

# [1] Import and preprocess analysis data

# Abricate Results
abricate <- read.csv("output/06_amr/abricate/amr/abricate.summary", sep = "\t")
abricate$X.FILE <- gsub(".tab", "", abricate$X.FILE)

# MegaRes database
megares.db <- read.csv("db/megares_drugs_annotations_v2.00.csv")
megares.db.filtered <- megares.db %>% separate(header, into = c("ID", "C1", "C2", "C3", "SYMBOL"), sep = "\\|", extra = "drop")
megares.db.filtered$group <- gsub("-", ".", megares.db.filtered$group)
megares.db.filtered <- megares.db.filtered[megares.db.filtered$group %in% colnames(abricate)[-c(1,2)],][,-c(1)]
megares.db.filtered <- unique(megares.db.filtered)
megares.db.clean <- data.frame(V1 = megares.db.filtered$C1, V2 = megares.db.filtered$C2, V3 = megares.db.filtered$C3, V4 = megares.db.filtered$group)
megares.db.clean <- megares.db.clean[order(megares.db.clean$V4),]
rownames(megares.db.clean) <- megares.db.clean$V4
megares.db.clean$V2 <- tolower(megares.db.clean$V2)

# Sequence Stats
seqreport = read.csv("output/01_preprocessing/seqfu/stats.tsv", header = TRUE, sep = "\t", row.names = 1)
tmp1 <- do.call(rbind, strsplit(rownames(seqreport), "/"))
tmp2 <- do.call(rbind, strsplit(data.frame(tmp1)$X4, "_"))
tmp3 <- data.frame(tmp2)$X1
tmp4 <- do.call(rbind, strsplit(data.frame(tmp2)$X2, ".fastq"))
tmp5 <- do.call(rbind, strsplit(data.frame(tmp4)$X1, "_"))
tmp6 <- data.frame(tmp5)$tmp5
seqreport$SAMPLE <- tmp3
seqreport$READ <- tmp6

# CoverM Results
rawCountTable <- read.table("output/06_amr/coverm/coverm.summary", header = TRUE, sep = "\t", row.names = 1)
colnames(rawCountTable) <- sub("output.06_amr.coverm.", "", colnames(rawCountTable))

# BIOM File
data.biom <- import_biom("output/02_taxonomic_profiling/kraken_biom/kraken2.biom", parseFunction = parse_taxonomy_default)

# Metadata
meta <- read.csv("metadata/16s_Horses_Overview_Reordered.csv", sep = "\t")
meta.filtered <- meta[meta$Microbiome == "Gut" & meta$Type == "WGS",]

# Group order
groups.order <- c("SSG", "5DG", "SWITCHED")

# --------------------------------------------------------------------------------------------------------

# [2] Diversity Estimations

# Sample Depth
sample_depth <- sort(sample_sums(data.biom))

# Alpha Diversity (Raw)
data.alpha <- microbiome::alpha(data.biom)

# Alpha diversity (Rarefy)
data.rarefy <- rarefy_even_depth(data.biom, rngseed = 1, sample.size = 0.9*min(sample_sums(data.biom)), replace = FALSE, trimOTUs = FALSE)
data.alpha.rarefy <- microbiome::alpha(data.rarefy)
#data.rarefy <- aggregate_top_taxa(data.rarefy, 22, "Rank2")

# Beta Diversity (Bray-Curtis distance)
braycurtis <- phyloseq::distance(data.rarefy, method = "bray")
data.bray <- as.matrix(braycurtis)

# Beta Diversity (PCoA)
braycurtis.pcoa <- ordinate(physeq = data.rarefy, method = "PCoA", distance = "bray")
data.pcoa <- as.data.frame(braycurtis.pcoa$vectors, row.names = NULL, optional = FALSE, cut.names = FALSE, col.names = names(braycurtis.pcoa$vectors), fix.empty.names = TRUE, stringsAsFactors = default.stringsAsFactors())

# Add Metadata
meta.sorted = meta.filtered[match(rownames(data.alpha), meta.filtered$SampleID),]
data.alpha$HORSE = meta.sorted$HorseID
data.alpha$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.alpha$TIMEPOINT = meta.sorted$Timepoint
data.alpha.rarefy$HORSE = meta.sorted$HorseID
data.alpha.rarefy$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.alpha.rarefy$TIMEPOINT = meta.sorted$Timepoint
data.pcoa$HORSE = meta.sorted$HorseID
data.pcoa$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.pcoa$TIMEPOINT = meta.sorted$Timepoint
data.pcoa$TIME_GROUP = paste(meta.sorted$Timepoint, meta.sorted$AB_Group, sep = )
data.pcoa.filtered <- data.pcoa[data.pcoa$AB_GROUP != "SWITCHED", ]

# Export Diversities
write.csv(data.alpha, file = "output/07_visualization/tab_div_alpha_raw.csv", quote = FALSE)
write.csv(data.alpha.rarefy, file = "output/07_visualization/tab_div_alpha_rarefy.csv", quote = FALSE)
write.csv(data.bray, file = "output/07_visualization/tab_div_beta_distance.csv", quote = FALSE)
write.csv(data.pcoa, file = "output/07_visualization/tab_div_beta_pcoa.csv", quote = FALSE)

# Export human-readable OTU table
data.otu <- phyloseq_to_df(data.rarefy)
write.csv(data.otu, file = "output/07_visualization/tab_otu.csv", row.names = FALSE, quote = FALSE)

# PCA
colours.days = c("t0" = "#00BA38",
                 "t1" = "#F8766D",
                 "t2" = "#619CFF"
                 )

colours.groups = c("SSG" = "#00ff7f",
                   "5DG" = "#ffa500",
                   "SWITCHED" = "#00bfff"
                   )

eigenvalue_pc1 = round(braycurtis.pcoa$values$Relative_eig[1]*100, 1)
eigenvalue_pc2 = round(braycurtis.pcoa$values$Relative_eig[2]*100, 1)

# Timepoints - All Samples
png("output/07_visualization/div_pca_time_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "All Samples - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - SSG
png("output/07_visualization/div_pca_time_ssg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered[data.pcoa.filtered$AB_GROUP == "SSG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "SSG - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - 5DG
png("output/07_visualization/div_pca_time_5dg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered[data.pcoa.filtered$AB_GROUP == "5DG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "5DG - Beta Diversity PCA (Bray–Curtis Dissimilarity)",
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - All Samples
png("output/07_visualization/div_pca_group_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups,
          title = "All Samples - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t0
png("output/07_visualization/div_pca_group_t0.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered[data.pcoa.filtered$TIMEPOINT == "t0", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups[1:2],
          title = "t0 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t1
png("output/07_visualization/div_pca_group_t1.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered[data.pcoa.filtered$TIMEPOINT == "t1", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups[1:2],
          title = "t1 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t2
png("output/07_visualization/div_pca_group_t2.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered[data.pcoa.filtered$TIMEPOINT == "t2", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups[1:2],
          title = "t2 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints-Groups - All Samples
png("output/07_visualization/div_pca_group_time_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa.filtered,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.15, 0.7),
          ylim = c(-0.35, 0.6),
          color = "TIME_GROUP",
          #shape = "TIME_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          title = "All Samples - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Boxplots
boxplot.comparisions <- list(c("t0", "t1"), c("t1", "t2"), c("t0","t2"))

png("output/07_visualization/div_box.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy, aes(x = TIMEPOINT, y = diversity_shannon, fill = TIMEPOINT)) +
        geom_boxplot() +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(2, 9)) +
        scale_y_continuous(breaks = c(3, 5, 7)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(7.5, 8.1, 8.7),
                           size = 3,
                           paired = TRUE)
dev.off()

png("output/07_visualization/div_box_even.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy, aes(x = TIMEPOINT, y = evenness_simpson, fill = TIMEPOINT)) +
        geom_boxplot() +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(0.0007, 0.03)) +
        scale_y_continuous(breaks = c(0.005, 0.015, 0.025)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(0.023, 0.0255, 0.028),
                           size = 3,
                           paired = TRUE)
dev.off()

# Barplots
bar_data_aggregated <- aggregate_top_taxa(data.rarefy, 9, "Rank2")
bar_data_melted <- psmelt(bar_data_aggregated)
num_taxa <- length(unique(bar_data_melted$OTU))
palette <- distinctColorPalette(num_taxa)

bar_ext_horse = c()
bar_ext_group = c()
bar_ext_time = c()
i = 1

for(e in bar_data_melted$Sample){
        current_line = data.alpha.rarefy[rownames(data.alpha.rarefy) == e,]
        bar_ext_horse[i] <- current_line$HORSE
        bar_ext_group[i] <- as.character(current_line$AB_GROUP)
        bar_ext_time[i] <- current_line$TIMEPOINT
        i = i + 1
        }

bar_data_melted$HORSE <- bar_ext_horse
bar_data_melted$AB_GROUP <- factor(bar_ext_group, levels = groups.order)
bar_data_melted$TIMEPOINT <- bar_ext_time

# Rename OTUs accordingly
bar_data_melted$OTU[bar_data_melted$OTU == "p__"] <- "Unclassified"
bar_data_melted$OTU <- gsub("p__", "", bar_data_melted$OTU)
#bar_data_melted$OTU <- gsub(" 1", "", bar_data_melted$OTU)
bar_data_melted$OTU <- bar_data_melted$OTU

# Prepare order based on taxa abundance
# Deleting this section will result in alphabetical order within the barcharts
#bar_order <- unique(bar_data_melted$OTU)
#bar_order <- bar_order[bar_order != "Other" & bar_order != "Unclassified"]
#bar_order[length(bar_order) + 1] <- "Other"
#bar_order[length(bar_order) + 1] <- "Unclassified"
#bar_data_melted$OTU <- factor(bar_data_melted$OTU, levels = bar_order)

# Individual Horses - 5DG
bar_data_5dg <- bar_data_melted[bar_data_melted$AB_GROUP == "5DG",]

png("output/07_visualization/tax_bar_horses_5dg.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_5dg, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Relative Abundance (5DG)",
             x = "Horses",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
             ) +
        facet_grid(~ HORSE) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              )
dev.off()

# Individual Horses - SSG
bar_data_ssg <- bar_data_melted[bar_data_melted$AB_GROUP == "SSG",]

png("output/07_visualization/tax_bar_horses_ssg.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_ssg, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Relative Abundance (SSG)",
             x = "Horses",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
             ) +
        facet_grid(~ HORSE) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              )
dev.off()

# Individual Horses - Switched
bar_data_switched <- bar_data_melted[bar_data_melted$AB_GROUP == "SWITCHED",]

png("output/07_visualization/tax_bar_horses_switched.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_switched, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Relative Abundance (Switched)",
             x = "Horses",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
             ) +
        facet_grid(~ HORSE) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              )
dev.off()

# Summarized groups
png("output/07_visualization/tax_bar_sum.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Mean Abundance (Groups)",
             x = "AB Groups",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
             ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_grid(~ AB_GROUP) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              )
dev.off()

# Percentages of individual taxa
per_taxa_total = sum(bar_data_melted$Abundance)/36

per_bact_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacteroidetes", ]$Abundance)/36
per_bact_norm = (per_bact_total/per_taxa_total) * 100

per_firm_total = sum(bar_data_melted[bar_data_melted$OTU == "Firmicutes", ]$Abundance)/36
per_firm_norm = (per_firm_total/per_taxa_total) * 100

per_prot_total = sum(bar_data_melted[bar_data_melted$OTU == "Proteobacteria", ]$Abundance)/36
per_prot_norm = (per_prot_total/per_taxa_total) * 100

per_spir_total = sum(bar_data_melted[bar_data_melted$OTU == "Spirochaetes", ]$Abundance)/36
per_spir_norm = (per_spir_total/per_taxa_total) * 100

per_verr_total = sum(bar_data_melted[bar_data_melted$OTU == "Verrucomicrobia", ]$Abundance)/36
per_verr_norm = (per_verr_total/per_taxa_total) * 100

per_bact_norm
per_firm_norm
per_prot_norm
per_spir_norm
per_verr_norm

# Abundance heatmap
abundance_matrix <- matrix(0, 9, length(unique(bar_data_melted$OTU)))
abundance_columns <- unique(bar_data_melted$OTU)
abundance_rows <- c("SSG_t0", "5DG_t0", "SWITCHED_t0", "SSG_t1", "5DG_t1", "SWITCHED_t1", "SSG_t2", "5DG_t2", "SWITCHED_t2")

j = 1

for (c in abundance_columns){
        i = 1
        for (r in abundance_rows){
                group <- strsplit(r, "_")[[1]][1]
                timepoint <- strsplit(r, "_")[[1]][2]
                abundance_matrix[i, j] <- mean(bar_data_melted[bar_data_melted$TIMEPOINT == timepoint & bar_data_melted$AB_GROUP == group & bar_data_melted$OTU == c,]$Abundance)
                i = i + 1
                }
        j = j + 1
        }

abundance_matrix <- t(abundance_matrix)

rownames(abundance_matrix) <- abundance_columns
colnames(abundance_matrix) <- c("SSG", "5DG", "SWITCHED", "SSG", "5DG", "SWITCHED", "SSG", "5DG", "SWITCHED")

png("output/07_visualization/tax_heat.png", width = 20, height = 15, units = "cm", res = 500)
Heatmap(log2(abundance_matrix),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        col = c("black", "darkred", "red", "orange", "yellow"),
        column_split = c(rep("t0", 3), rep("t1", 3), rep("t2", 3)),
        row_title = "Top 10 Phyla",
        name = "log2(Abundance)",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
        )
dev.off()

# --------------------------------------------------------------------------------------------------------

# [3] AMR Heatmap (Assembly)

# Count individual resistance genes
abricate.count <- abricate
for (i in 1:36){
        for (j in 3:ncol(abricate)){
                d = abricate[i,j]
                if (d == "."){
                        k = 0
                } else{
                        k = str_count(d,"\\.")
                }
                abricate.count[i,j] <- k
        }
}

abricate.filtered <- abricate.count[,-c(1,2)]
rownames(abricate.filtered) <- abricate[,1]
abricate.matrix <- data.matrix(abricate.filtered) - 1
amr.counts <- rowSums(abricate.matrix)
abricate.meta <- meta.filtered[match(rownames(abricate.matrix), meta.filtered$SampleID),]
abricate.meta$AMR_FOUND <- amr.counts
abricate.meta$DIV = data.alpha.rarefy$diversity_shannon

# Update AMR counts
amr.df.counts <- data.frame(ids = rownames(abricate), counts = amr.counts)
seqreport.filtered <- seqreport[seqreport$READ == "R1",]
rownames(seqreport.filtered) <- seqreport.filtered$SAMPLE
amr.df.stats <- merge(x = seqreport.filtered, y = amr.df.counts, by = "row.names")
amr.df.stats <- subset(amr.df.stats, select = c("SAMPLE", "X.Seq", "counts"))

# Read AMR Normalization by #Reads
amr.df.stats$X.Seq <- amr.df.stats$X.Seq*2
amr.df.stats$Norm <- amr.df.stats$counts/amr.df.stats$X.Seq
amr.df.stats$CPM <- amr.df.stats$counts/(amr.df.stats$X.Seq/1000000)
amr.df.stats$CP60M <- amr.df.stats$CPM*60
colnames(amr.df.stats) <- c("SAMPLE", "#PE_READS", "AMR_COUNTS", "NORM_COUNT", "CPM", "CP60M")
write.csv(amr.df.stats, file = "output/07_visualization/tab_amr_counts.csv", quote = FALSE)

# Prepare dataframe for heatmap
amr.norm.reads <- subset(amr.df.stats, select = c("SAMPLE", "CPM", "CP60M"))
amr.norm.reads$DIV = data.alpha.rarefy$diversity_shannon
amr.norm.reads$TIMEPOINT = data.alpha.rarefy$TIMEPOINT
amr.norm.reads$AB_GROUP = data.alpha.rarefy$AB_GROUP
amr.norm.reads.filtered <- amr.norm.reads[amr.norm.reads$AB_GROUP != "SWITCHED",]

# Define colours for heatmap
colours.amr = c("aminoglycosides" = "#e6194B",
                "betalactams" = "#3cb44b",
                "drug_and_biocide_resistance" = "#ffe119",
                "drug_and_biocide_and_metal_resistance" = "#4363d8",
                "rifampin" = "#f58231",
                "multi-drug_resistance" = "#911eb4",
                "bacitracin" = "#42d4f4",
                "phenicol" = "#f032e6",
                "trimethoprim" = "#bfef45",
                "cationic_antimicrobial_peptides" = "#fabed4",
                "mls" = "#469990",
                "fosfomycin" = "#dcbeff",
                "lipopeptides" = "#9A6324",
                "metronidazole" = "#fffac8",
                "fluoroquinolones" = "#800000",
                "nucleosides" = "#aaffc3",
                "sulfonamides" = "#808000",
                "tetracyclines" = "#ffd8b1",
                "glycopeptides" = "#000075"
                )

colours.genes = colorRamp2(c(1, 50, 100, 180), c("grey", "yellow", "orange", "red"))

# Set Annotations for Heatmap
annot.column = HeatmapAnnotation(class = megares.db.clean$V2,
                                 col = list(class = colours.amr)
                                 )

annot.row.left = rowAnnotation(timepoint = abricate.meta$Timepoint,
                               ab_group = abricate.meta$AB_Group,
                               "#amr_genes" = amr.counts,
                               col = list(ab_group = colours.groups,
                                          timepoint = colours.days,
                                          "#amr_genes" = colours.genes
                                          )
                               )

annot.row.right = rowAnnotation(horse = anno_text(abricate.meta$HorseID,
                                                   gp = gpar(fontsize = 8)
                                                  )
                                )

re.order.cols = megares.db.clean[order(megares.db.clean$V2),]
re.order.rows <- abricate.meta[with(abricate.meta, order(Day, AB_Group, HorseID)), ]

# Group all counts > 1 into a single class
abricate.matrix[abricate.matrix >= 2] = 2

# Plot Heatmap
png("output/07_visualization/amr_heat_abricate.png", width = 40, height = 20, units = "cm", res = 500)
Heatmap(abricate.matrix,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        #col = c("black", "red4", "firebrick3", "firebrick2", "red1"),
        col = c("black", "gold", "red"),
        
        column_names_gp = gpar(fontsize = 5),
        column_order = re.order.cols$V4,
        column_title = "AMR Genes",
        
        row_order = re.order.rows$SampleID,
        row_split = abricate.meta$Timepoint,
        row_title = "Metagenome Samples",
        show_row_names = FALSE,
        
        heatmap_legend_param = list(
                title = "#copies", at = c(0, 1, 2),
                labels = c("0", "1", "\u2265 2")),
        
        left_annotation = annot.row.left,
        right_annotation = annot.row.right,
        top_annotation = annot.column
        )
dev.off()

# --------------------------------------------------------------------------------------------------------

# [4] Statistical Correlation (Assembly)

# Correlation - AMR and Diversity
cor.test(amr.norm.reads$CP60M, amr.norm.reads$DIV, method = "spearman")

png("output/07_visualization/amr_div_cor.png", width = 17, height = 16, units = "cm", res = 500)
ggscatter(amr.norm.reads,
          x = "CP60M",
          y = "DIV",
          color = "TIMEPOINT",
          shape = "AB_GROUP",
          xlab = "AMR PER 60M READS",
          ylab = "SHANNON INDEX",
          add = "reg.line",
          add.params = list(color = "black",
                            fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman",
                                label.x = 100,
                                label.sep = "\n"),
          title = "CORRELATION - AMR / DIVERSITY - ASSEMBLY",
          palette = colours.days,
          size = 2.5
          ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(shape = "AB GROUP")
dev.off()

# Wilcoxon Rank Sum Tests

# Diversity - SSG

# Significant Differences between t0 and t1 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t2" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t0" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t1" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Diversity - 5DG

# Significant Differences between t0 and t1 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t2" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t0" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t1" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, p.adjust.method = "BH", paired = TRUE)
stat.res

# Diversity - SSG vs. 5DG

# Significant Differences at t0-t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t0 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t1 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t2" & data.alpha.rarefy$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# AMR - SSG

# Significant Differences between t0 and t1 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t2" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t0" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t1" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# AMR - 5DG

# Significant Differences between t0 and t1 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t2" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t0" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t1" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH", paired = TRUE)
stat.res

# AMR - SSG vs. 5DG

# Significant Differences at t0-t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t0 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t0" & amr.norm.reads$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t1 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t1" & amr.norm.reads$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences at t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t2" & amr.norm.reads$AB_GROUP != "SWITCHED",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# --------------------------------------------------------------------------------------------------------

# [5] AMR Normalization (Reads)

# [5.1] AMR Gene Heatmap

# Read raw data
Gene_sumCountTable <- aggregate(. ~ Resistance_Symbol, rawCountTable[,-c(1:3)], sum)
rownames(Gene_sumCountTable) <- Gene_sumCountTable$Resistance_Symbol
Gene_filterCountTable <- Gene_sumCountTable[, -c(1)]

# Perform edgeR TMM normalization
Gene_dgeFullRaw <- DGEList(Gene_filterCountTable, group = colnames(Gene_sumCountTable)[-c(1)])
Gene_dgeFullRaw <- edgeR::calcNormFactors(Gene_dgeFullRaw, method = "TMM")

# Filter 0 values
Gene_Keep <- filterByExpr(Gene_dgeFullRaw)
Gene_dgeFullFilter <- Gene_dgeFullRaw[Gene_Keep, , keep.lib.sizes = FALSE]
Gene_dgeFullFilter <- edgeR::calcNormFactors(Gene_dgeFullFilter, method = "TMM")

# Output TMM values
Gene_normCounts <- cpm(Gene_dgeFullFilter)
write.csv(Gene_normCounts, file = "output/07_visualization/tab_amr_tmm.csv", quote = FALSE)

# Bin TMM values by counts
Gene_groupedCounts = Gene_normCounts
Gene_groupedCounts[Gene_groupedCounts < 500] = 0
Gene_groupedCounts[Gene_groupedCounts >= 500 & Gene_groupedCounts < 1000] = 1
Gene_groupedCounts[Gene_groupedCounts >= 1000 & Gene_groupedCounts < 5000] = 2
Gene_groupedCounts[Gene_groupedCounts >= 5000] = 3
# Filter AMR genes to make heatmap less dense
# Might want to improve this with a better filtering mechanism!
Gene_groupedCounts <- subset(Gene_groupedCounts, rowSums(Gene_groupedCounts) > 10)
Gene_groupedCounts <- t(Gene_groupedCounts)

# Prepare metadata for AMR Read Heatmap
g = c()
h = c()
j = 1

for (i in colnames(Gene_groupedCounts)) {
        g[j] <- subset(rawCountTable, Resistance_Symbol == i)$Resistance_Class[1]
        h[j] <- i
        j = j + 1
        }

Gene_AMR_DF = data.frame(AMR_Gene = h, AMR_Class = tolower(g))
re.order.cols.reads = Gene_AMR_DF[order(Gene_AMR_DF$AMR_Class),]

# Plot AMR Read Heatmap
png("output/07_visualization/amr_heat_coverm_genes.png", width = 40, height = 20, units = "cm", res = 500)
Heatmap(Gene_groupedCounts,
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        col = c("black", "red", "orange", "yellow"),
        
        column_names_gp = gpar(fontsize = 4),
        #column_order = re.order.cols.reads$AMR_Gene,
        column_title = "AMR Genes",
        
        row_order = re.order.rows$SampleID,
        row_split = abricate.meta$Timepoint,
        row_title = "Metagenome Samples",
        show_row_names = FALSE,
        
        heatmap_legend_param = list(
                title = "abundance", at = c(0, 1, 2, 3), 
                labels = c("<500", "500-1000", "1000-5000", ">5000")),
        
        left_annotation = rowAnnotation(timepoint = abricate.meta$Timepoint,
                                        ab_group = abricate.meta$AB_Group,
                                        col = list(ab_group = colours.groups,
                                                   timepoint = colours.days)
                                        ),
        
        right_annotation = rowAnnotation(horse = anno_text(abricate.meta$HorseID,
                                                           gp = gpar(fontsize = 8))),
        
        top_annotation = columnAnnotation(amr_class = Gene_AMR_DF$AMR_Class)
        )
dev.off()

# [5.2] AMR Class Heatmap

# Read raw data
Class_sumCountTable <- aggregate(. ~ Resistance_Class, rawCountTable[,-c(1,3,4)], sum)
rownames(Class_sumCountTable) <- Class_sumCountTable$Resistance_Class
Class_filterCountTable <- Class_sumCountTable[, -c(1)]

# Perform edgeR TMM normalization
Class_dgeFullRaw <- DGEList(Class_filterCountTable, group = colnames(Class_sumCountTable)[-c(1)])
Class_dgeFullRaw <- edgeR::calcNormFactors(Class_dgeFullRaw, method = "TMM")

# Filter 0 values
Class_Keep <- filterByExpr(Class_dgeFullRaw)
Class_dgeFullFilter <- Class_dgeFullRaw[Class_Keep, , keep.lib.sizes = FALSE]
Class_dgeFullFilter <- edgeR::calcNormFactors(Class_dgeFullFilter, method = "TMM")

# Output TMM values
Class_normCounts <- cpm(Class_dgeFullFilter)
write.csv(Class_normCounts, file = "output/07_visualization/tab_amr_tmm_classes.csv", quote = FALSE)

# Bin TMM values by counts
Class_groupedCounts = Class_normCounts
Class_groupedCounts[Class_groupedCounts < 500] = 0
Class_groupedCounts[Class_groupedCounts >= 500 & Class_groupedCounts < 1000] = 1
Class_groupedCounts[Class_groupedCounts >= 1000 & Class_groupedCounts < 5000] = 2
Class_groupedCounts[Class_groupedCounts >= 5000] = 3
Class_groupedCounts <- t(Class_groupedCounts)
colnames(Class_groupedCounts) <- tolower(colnames(Class_groupedCounts))

# Plot AMR Class Heatmap
png("output/07_visualization/amr_heat_coverm_classes.png", width = 30, height = 20, units = "cm", res = 500)
Heatmap(Class_groupedCounts,
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        col = c("black", "red", "orange", "yellow"),
        
        column_title = "AMR Classes",
        column_names_max_height = unit(10, "cm"),
        column_names_gp = gpar(fontsize = 8),
        
        row_title = "Metagenome Samples",
        row_order = re.order.rows$SampleID,
        row_split = abricate.meta$Timepoint,
        show_row_names = FALSE,
        
        heatmap_legend_param = list(
                title = "abundance", at = c(0, 1, 2, 3), 
                labels = c("<500", "500-1000", "1000-5000", ">5000")),
        
        left_annotation = rowAnnotation(timepoint = abricate.meta$Timepoint,
                                        ab_group = abricate.meta$AB_Group,
                                        col = list(ab_group = colours.groups,
                                                   timepoint = colours.days)
                                        ),
        
        right_annotation = rowAnnotation(horse = anno_text(abricate.meta$HorseID,
                                         gp = gpar(fontsize = 8)))
        )
dev.off()

# [5.3] Scatterplots / Boxplots of AMR Classes
e = c()
f = c()
j = 1

for (i in rownames(data.alpha.rarefy)) {
        e[j] <- sum(Class_normCounts[,i])
        f[j] <- i
        j = j + 1
        }

Class_AMR_SUM <- data.frame(SAMPLE = f, AMR = e, TIMEPOINT = data.alpha.rarefy$TIMEPOINT, AB_GROUP = data.alpha.rarefy$AB_GROUP, HORSE = data.alpha.rarefy$HORSE)

png("output/07_visualization/amr_sum_scatter.png", width = 30, height = 20, units = "cm", res = 500)
ggplot(Class_AMR_SUM, aes(x = TIMEPOINT, y = log2(AMR))) +
        coord_cartesian(ylim = c(17, 24)) +
        geom_line(aes(group = HORSE)) +
        geom_point() +
        geom_point(size = 2) +
        facet_wrap(~AB_GROUP + HORSE, scale = "free")
dev.off()

# Total AMR Sum
png("output/07_visualization/amr_sum_box.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(Class_AMR_SUM, aes(x = TIMEPOINT, y = log2(AMR), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(17, 25)) +
        geom_boxplot() +
        facet_wrap(~AB_GROUP, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(22.5, 23.5, 24.5),
                           size = 3,
                           paired = TRUE)
dev.off()

amr.class.df <- data.frame(t(Class_normCounts))
amr.class.df$AB_GROUP <- data.pcoa$AB_GROUP
amr.class.df$HORSE <- data.pcoa$HORSE
amr.class.df$TIMEPOINT <- data.pcoa$TIMEPOINT
amr.class.df$TIME_GROUP <- data.pcoa$TIME_GROUP

# SSG
amr.class.df.filtered <- amr.class.df[amr.class.df$AB_GROUP == "SSG", ]
amr.class.length <- length(colnames(amr.class.df.filtered)) - 4
amr.class.df.ssg <- data.frame(AMR_Class = rep(colnames(amr.class.df.filtered[1:amr.class.length]), each = 15),
                               AMR_TMM = unlist(amr.class.df.filtered[1:amr.class.length]),
                               AB_GROUP = rep(amr.class.df.filtered$AB_GROUP, amr.class.length),
                               HORSE = rep(amr.class.df.filtered$HORSE, amr.class.length),
                               TIMEPOINT = rep(amr.class.df.filtered$TIMEPOINT, amr.class.length),
                               TIME_GROUP = rep(amr.class.df.filtered$TIME_GROUP, amr.class.length)
                               )

png("output/07_visualization/amr_box_group_ssg.png", width = 40, height = 30, units = "cm", res = 500)
ggplot(amr.class.df.ssg, aes(x = TIMEPOINT, y = log2(AMR_TMM + 1), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(0, 25)) +
        geom_boxplot() +
        facet_wrap(~AMR_Class, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(21, 22.5, 24),
                           size = 2,
                           paired = TRUE,
                           method.args = list(exact = FALSE))
dev.off()

# 5DG
amr.class.df.filtered <- amr.class.df[amr.class.df$AB_GROUP == "5DG", ]
amr.class.length <- length(colnames(amr.class.df.filtered)) - 4
amr.class.df.5dg <- data.frame(AMR_Class = rep(colnames(amr.class.df.filtered[1:amr.class.length]), each = 15),
                               AMR_TMM = unlist(amr.class.df.filtered[1:amr.class.length]),
                               AB_GROUP = rep(amr.class.df.filtered$AB_GROUP, amr.class.length),
                               HORSE = rep(amr.class.df.filtered$HORSE, amr.class.length),
                               TIMEPOINT = rep(amr.class.df.filtered$TIMEPOINT, amr.class.length),
                               TIME_GROUP = rep(amr.class.df.filtered$TIME_GROUP, amr.class.length)
                               )

png("output/07_visualization/amr_box_group_5dg.png", width = 40, height = 30, units = "cm", res = 500)
ggplot(amr.class.df.5dg, aes(x = TIMEPOINT, y = log2(AMR_TMM + 1), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(0, 25)) +
        geom_boxplot() +
        facet_wrap(~AMR_Class, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(21, 22.5, 24),
                           size = 2,
                           paired = TRUE,
                           method.args = list(exact = FALSE))
dev.off()

# SWITCHED
amr.class.df.filtered <- amr.class.df[amr.class.df$AB_GROUP == "SWITCHED", ]
amr.class.length <- length(colnames(amr.class.df.filtered)) - 4
amr.class.df.switched <- data.frame(AMR_Class = rep(colnames(amr.class.df.filtered[1:amr.class.length]), each = 6),
                                    AMR_TMM = unlist(amr.class.df.filtered[1:amr.class.length]),
                                    AB_GROUP = rep(amr.class.df.filtered$AB_GROUP, amr.class.length),
                                    HORSE = rep(amr.class.df.filtered$HORSE, amr.class.length),
                                    TIMEPOINT = rep(amr.class.df.filtered$TIMEPOINT, amr.class.length),
                                    TIME_GROUP = rep(amr.class.df.filtered$TIME_GROUP, amr.class.length)
                                    )

png("output/07_visualization/amr_box_group_switched.png", width = 40, height = 30, units = "cm", res = 500)
ggplot(amr.class.df.switched, aes(x = TIMEPOINT, y = log2(AMR_TMM + 1), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(0, 25)) +
        geom_boxplot() +
        facet_wrap(~AMR_Class, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(21, 22.5, 24),
                           size = 2,
                           paired = TRUE,
                           method.args = list(exact = FALSE))
dev.off()

# [5.4] Statistical Correlation (Reads)
Class_AMR_SUM$AMR_LOG2 <- log2(Class_AMR_SUM$AMR)
Class_AMR_SUM$DIV = amr.norm.reads$DIV

cor.test(Class_AMR_SUM$AMR, Class_AMR_SUM$DIV, method = "spearman")

png("output/07_visualization/amr_div_cor_reads.png", width = 17, height = 16, units = "cm", res = 500)
ggscatter(Class_AMR_SUM,
          x = "AMR_LOG2",
          y = "DIV",
          color = "TIMEPOINT",
          shape = "AB_GROUP",
          xlab = "LOG2(AMR SUM TMM)",
          ylab = "SHANNON INDEX",
          add = "reg.line",
          add.params = list(color = "black",
                            fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman",
                                label.sep = "\n"),
          title = "CORRELATION - AMR / DIVERSITY - READS",
          palette = colours.days,
          size = 2.5
          ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(shape = "AB GROUP")
dev.off()

# --------------------------------------------------------------------------------------------------------

# [6] Differential Taxa Analysis

# Family Level
edger.aggretated <- aggregate_top_taxa(data.biom, 10^5, "Rank5")

# Prepare data for edgeR analysis
edger.groups <- data.pcoa[with(data.pcoa, order(rownames(data.pcoa))), ]$TIME_GROUP
edger.dge <- phyloseq_to_edgeR(edger.aggretated, group = edger.groups)

# C1 - SSG: t0 vs. t1
edger.res.c1 <- exactTest(edger.dge, pair = c("t0 SSG", "t1 SSG"))
edger.top.c1 <- topTags(edger.res.c1, n = nrow(edger.res.c1$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c1 <- edger.top.c1$table[edger.top.c1$table$FDR < 0.05 & abs(edger.top.c1$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_ssg_t0_t1.csv", quote = FALSE)

# C2 - SSG: t1 vs. t2
edger.res.c2 <- exactTest(edger.dge, pair = c("t1 SSG", "t2 SSG"))
edger.top.c2 <- topTags(edger.res.c2, n = nrow(edger.res.c2$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c2 <- edger.top.c2$table[edger.top.c2$table$FDR < 0.05 & abs(edger.top.c2$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_ssg_t1_t2.csv", quote = FALSE)

# C3 - SSG: t0 vs. t2
edger.res.c3 <- exactTest(edger.dge, pair = c("t0 SSG", "t2 SSG"))
edger.top.c3 <- topTags(edger.res.c3, n = nrow(edger.res.c3$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c3 <- edger.top.c3$table[edger.top.c3$table$FDR < 0.05 & abs(edger.top.c3$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_ssg_t0_t2.csv", quote = FALSE)

# C4 - 5DG: t0 vs. t1
edger.res.c4 <- exactTest(edger.dge, pair = c("t0 5DG", "t1 5DG"))
edger.top.c4 <- topTags(edger.res.c4, n = nrow(edger.res.c4$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c4 <- edger.top.c4$table[edger.top.c4$table$FDR < 0.05 & abs(edger.top.c4$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_5dg_t0_t1.csv", quote = FALSE)

# C5 - 5DG: t1 vs. t2
edger.res.c5 <- exactTest(edger.dge, pair = c("t1 5DG", "t2 5DG"))
edger.top.c5 <- topTags(edger.res.c5, n = nrow(edger.res.c5$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c5 <- edger.top.c5$table[edger.top.c5$table$FDR < 0.05 & abs(edger.top.c5$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_5dg_t1_t2.csv", quote = FALSE)

# C6 - 5DG: t0 vs. t2
edger.res.c6 <- exactTest(edger.dge, pair = c("t0 5DG", "t2 5DG"))
edger.top.c6 <- topTags(edger.res.c6, n = nrow(edger.res.c6$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c6 <- edger.top.c6$table[edger.top.c6$table$FDR < 0.05 & abs(edger.top.c6$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_5dg_t0_t2.csv", quote = FALSE)

# C7 - t0: SSG vs. 5DG
edger.res.c7 <- exactTest(edger.dge, pair = c("t0 SSG", "t0 5DG"))
edger.top.c7 <- topTags(edger.res.c7, n = nrow(edger.res.c7$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c7 <- edger.top.c7$table[edger.top.c7$table$FDR < 0.05 & abs(edger.top.c7$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_t0_ssg_5dg.csv", quote = FALSE)

# C8 - t1: SSG vs. 5DG
edger.res.c8 <- exactTest(edger.dge, pair = c("t1 SSG", "t1 5DG"))
edger.top.c8 <- topTags(edger.res.c8, n = nrow(edger.res.c8$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c8 <- edger.top.c8$table[edger.top.c8$table$FDR < 0.05 & abs(edger.top.c8$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_t1_ssg_5dg.csv", quote = FALSE)

# C9 - t2: SSG vs. 5DG
edger.res.c9 <- exactTest(edger.dge, pair = c("t2 SSG", "t2 5DG"))
edger.top.c9 <- topTags(edger.res.c9, n = nrow(edger.res.c9$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c9 <- edger.top.c9$table[edger.top.c9$table$FDR < 0.05 & abs(edger.top.c9$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/07_visualization/tab_edger_t2_ssg_5dg.csv", quote = FALSE)

# --------------------------------------------------------------------------------------------------------

# [7] Boxplots for abundance of specific families

boxplot.families.biom <- aggregate_top_taxa(data.rarefy, 22, "Rank5")
boxplot.families.melted <- psmelt(boxplot.families.biom)

box.ext.horse = c()
box.ext.group = c()
box.ext.time = c()
j = 1

for(f in boxplot.families.melted$Sample){
        current.line = data.alpha.rarefy[rownames(data.alpha.rarefy) == f,]
        box.ext.horse[j] <- current.line$HORSE
        box.ext.group[j] <- as.character(current.line$AB_GROUP)
        box.ext.time[j] <- current.line$TIMEPOINT
        j = j + 1
}

boxplot.families.melted$HORSE <- box.ext.horse
boxplot.families.melted$AB_GROUP <- factor(box.ext.group, levels = groups.order)
boxplot.families.melted$TIMEPOINT <- box.ext.time

G2 <- boxplot.families.melted[boxplot.families.melted$Rank5 == "f__Enterobacteriaceae", ]
png("output/07_visualization/tax_box_f_enterobacteriaceae.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(G2, aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
        geom_boxplot() +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(9, 19.5)) +
        scale_y_continuous(breaks = c(9, 12, 15, 18)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.comparisions,
                           method = "wilcox.test",
                           label.y = c(17, 18, 19),
                           size = 3,
                           paired = TRUE,
                           method.args = list(exact = FALSE)) +
        ggtitle("Abundance - Enterobacteriaceae (Family)") +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# --------------------------------------------------------------------------------------------------------

# [8] Calculate differences between timepoints

# Alpha diversity
alpha.5dg.t0 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0", ]$diversity_shannon)
alpha.5dg.t1 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t1", ]$diversity_shannon)
alpha.5dg.t2 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t2", ]$diversity_shannon)

alpha.5dg.diff.t0.t1 = round(alpha.5dg.t1 - alpha.5dg.t0, 2)
alpha.5dg.diff.t1.t2 = round(alpha.5dg.t2 - alpha.5dg.t1, 2)

alpha.ssg.t0 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0", ]$diversity_shannon)
alpha.ssg.t1 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t1", ]$diversity_shannon)
alpha.ssg.t2 = mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t2", ]$diversity_shannon)

alpha.ssg.diff.t0.t1 = round(alpha.ssg.t1 - alpha.ssg.t0, 2)
alpha.ssg.diff.t1.t2 = round(alpha.ssg.t2 - alpha.ssg.t1, 2)

# ARG abundance (assembly)
arg.5dg.t0 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "5DG" & amr.norm.reads.filtered$TIMEPOINT == "t0", ]$CP60M)
arg.5dg.t1 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "5DG" & amr.norm.reads.filtered$TIMEPOINT == "t1", ]$CP60M)
arg.5dg.t2 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "5DG" & amr.norm.reads.filtered$TIMEPOINT == "t2", ]$CP60M)

arg.5dg.diff.t0.t1 = round(arg.5dg.t1 - arg.5dg.t0, 2)
arg.5dg.diff.t0.t2 = round(arg.5dg.t2 - arg.5dg.t1, 2)

arg.ssg.t0 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "SSG" & amr.norm.reads.filtered$TIMEPOINT == "t0", ]$CP60M)
arg.ssg.t1 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "SSG" & amr.norm.reads.filtered$TIMEPOINT == "t1", ]$CP60M)
arg.ssg.t2 = mean(amr.norm.reads.filtered[amr.norm.reads.filtered$AB_GROUP == "SSG" & amr.norm.reads.filtered$TIMEPOINT == "t2", ]$CP60M)

arg.ssg.diff.t0.t1 = round(arg.ssg.t1 - arg.ssg.t0, 2)
arg.ssg.diff.t0.t2 = round(arg.ssg.t2 - arg.ssg.t1, 2)

# Abundance of Enterobacteriaceae
entero.5dg.t0 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t0", ]$Abundance))
entero.5dg.t1 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t1", ]$Abundance))
entero.5dg.t2 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t2", ]$Abundance))

entero.5dg.diff.t0.t1 = round(entero.5dg.t1 - entero.5dg.t0, 2)
entero.5dg.diff.t1.t2 = round(entero.5dg.t2 - entero.5dg.t1, 2)

entero.ssg.t0 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t0", ]$Abundance))
entero.ssg.t1 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t1", ]$Abundance))
entero.ssg.t2 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t2", ]$Abundance))

entero.ssg.diff.t0.t1 = round(entero.ssg.t1 - entero.ssg.t0, 2)
entero.ssg.diff.t1.t2 = round(entero.ssg.t2 - entero.ssg.t1, 2)

# --------------------------------------------------------------------------------------------------------

# [9] Visualizing ARG abundance per family

# Prepare input data
kreport = read.csv("output/06_amr/abricate/amr/kraken2.summary", header = TRUE, sep = "\t", row.names = 1)
taxids <- sub("X", "", colnames(kreport))
colnames(kreport) <- taxids

# Generate database
prepareDatabase(sqlFile = "tmp/accessionTaxa.sql", tmpDir = "tmp/")

# Annotate taxids
taxnames <- as.data.frame(getTaxonomy(taxids, sqlFile = "tmp/accessionTaxa.sql", desiredTaxa = "family"))
taxnames$taxid <- taxids

# Remove ARGs of higher classification level
# And merge identical taxas
taxa.matrix <- matrix(, nrow = nrow(kreport), ncol = length(unique(taxnames$family)) - 1)
names = c()
c = 1
for (t in unique(taxnames$family)){
        if (!is.na(t)) {
                names[c] <- t
                taxa.list <- taxnames[taxnames$family == t & !is.na(taxnames$family),]$taxid
                taxa.sums <- rowSums(kreport[taxa.list])
                taxa.matrix[,c] <- taxa.sums
                c = c + 1
        }
}

colnames(taxa.matrix) <- names
rownames(taxa.matrix) <- names(taxa.sums)

# Select top taxas and merge lower taxas
taxa.filter.top <- names(sort(colSums(taxa.matrix), T)[1:9])
taxa.filter.other <- names(sort(colSums(taxa.matrix), T)[10:ncol(taxa.matrix)])
taxa.matrix.filtered <- as.data.frame(taxa.matrix[, colnames(taxa.matrix) %in% taxa.filter.top])
taxa.matrix.filtered$Other <- rowSums(taxa.matrix[, colnames(taxa.matrix) %in% taxa.filter.other])

# Add sample metadata
taxa.matrix.filtered$AB_GROUP <- data.alpha.rarefy$AB_GROUP
taxa.matrix.filtered$HORSE <- data.alpha.rarefy$HORSE
taxa.matrix.filtered$TIMEPOINT <- data.alpha.rarefy$TIMEPOINT

# Create dataframe for visualization
taxa.g1.t0 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SSG" & taxa.matrix.filtered$TIMEPOINT == "t0",]
taxa.g1.t1 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SSG" & taxa.matrix.filtered$TIMEPOINT == "t1",]
taxa.g1.t2 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SSG" & taxa.matrix.filtered$TIMEPOINT == "t2",]

taxa.g2.t0 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "5DG" & taxa.matrix.filtered$TIMEPOINT == "t0",]
taxa.g2.t1 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "5DG" & taxa.matrix.filtered$TIMEPOINT == "t1",]
taxa.g2.t2 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "5DG" & taxa.matrix.filtered$TIMEPOINT == "t2",]

taxa.g3.t0 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SWITCHED" & taxa.matrix.filtered$TIMEPOINT == "t0",]
taxa.g3.t1 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SWITCHED" & taxa.matrix.filtered$TIMEPOINT == "t1",]
taxa.g3.t2 <- taxa.matrix.filtered[taxa.matrix.filtered$AB_GROUP == "SWITCHED" & taxa.matrix.filtered$TIMEPOINT == "t2",]

# Bubbleplots (Absolute)
taxa.g1.t0.nom <- taxa.g1.t0
taxa.g1.t1.nom <- taxa.g1.t1
taxa.g1.t2.nom <- taxa.g1.t2

taxa.g2.t0.nom <- taxa.g2.t0
taxa.g2.t1.nom <- taxa.g2.t1
taxa.g2.t2.nom <- taxa.g2.t2

taxa.g3.t0.nom <- taxa.g3.t0
taxa.g3.t1.nom <- taxa.g3.t1
taxa.g3.t2.nom <- taxa.g3.t2

for (s in unique(seqreport$SAMPLE)) {
        r = sum(seqreport[seqreport$SAMPLE == s & seqreport$READ != "R3", ]$X.Seq)
        if (s %in% rownames(taxa.g1.t0.nom)){
                taxa.g1.t0.nom[rownames(taxa.g1.t0.nom) == s, ][1:10] <- ((taxa.g1.t0.nom[rownames(taxa.g1.t0.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g1.t1.nom)){
                taxa.g1.t1.nom[rownames(taxa.g1.t1.nom) == s, ][1:10] <- ((taxa.g1.t1.nom[rownames(taxa.g1.t1.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g1.t2.nom)){
                taxa.g1.t2.nom[rownames(taxa.g1.t2.nom) == s, ][1:10] <- ((taxa.g1.t2.nom[rownames(taxa.g1.t2.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g2.t0.nom)){
                taxa.g2.t0.nom[rownames(taxa.g2.t0.nom) == s, ][1:10] <- ((taxa.g2.t0.nom[rownames(taxa.g2.t0.nom) == s, ][1:10]/r)*60000000) 
        }
        if (s %in% rownames(taxa.g2.t1.nom)){
                taxa.g2.t1.nom[rownames(taxa.g2.t1.nom) == s, ][1:10] <- ((taxa.g2.t1.nom[rownames(taxa.g2.t1.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g2.t2.nom)){
                taxa.g2.t2.nom[rownames(taxa.g2.t2.nom) == s, ][1:10] <- ((taxa.g2.t2.nom[rownames(taxa.g2.t2.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g3.t0.nom)){
                taxa.g3.t0.nom[rownames(taxa.g3.t0.nom) == s, ][1:10] <- ((taxa.g3.t0.nom[rownames(taxa.g3.t0.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g3.t1.nom)){
                taxa.g3.t1.nom[rownames(taxa.g3.t1.nom) == s, ][1:10] <- ((taxa.g3.t1.nom[rownames(taxa.g3.t1.nom) == s, ][1:10]/r)*60000000)
        }
        if (s %in% rownames(taxa.g3.t2.nom)){
                taxa.g3.t2.nom[rownames(taxa.g3.t2.nom) == s, ][1:10] <- ((taxa.g3.t2.nom[rownames(taxa.g3.t2.nom) == s, ][1:10]/r)*60000000)
        }
}

taxa.g1.abs.df <- data.frame(SIZE = c(colSums(taxa.g1.t0.nom[1:10]),
                                      colSums(taxa.g1.t1.nom[1:10]),
                                      colSums(taxa.g1.t2.nom[1:10])),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g1.t0.nom)[1:10], 3), levels = sort(colnames(taxa.g1.t0.nom), decreasing = TRUE))
                             )

taxa.g2.abs.df <- data.frame(SIZE = c(colSums(taxa.g2.t0.nom[1:10]),
                                      colSums(taxa.g2.t1.nom[1:10]),
                                      colSums(taxa.g2.t2.nom[1:10])),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g2.t0.nom)[1:10], 3), levels = sort(colnames(taxa.g2.t0.nom), decreasing = TRUE))
                             )

taxa.g3.abs.df <- data.frame(SIZE = c(colSums(taxa.g3.t0.nom[1:10]),
                                      colSums(taxa.g3.t1.nom[1:10]),
                                      colSums(taxa.g3.t2.nom[1:10])),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g3.t0.nom)[1:10], 3), levels = sort(colnames(taxa.g3.t0.nom), decreasing = TRUE))
                             )

png("output/07_visualization/amr_bubble_families_abs_ssg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g1.abs.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("SSG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 385), range = c(1, 12), breaks = c(10, 100, 200, 300), name = "ARG Count (CP60M)")
dev.off()

png("output/07_visualization/amr_bubble_families_abs_5dg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g2.abs.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("5DG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 385), range = c(1, 12), breaks = c(10, 100, 200, 300), name = "ARG Count (CP60M)")
dev.off()

png("output/07_visualization/amr_bubble_families_abs_switched.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g3.abs.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("SWITCHED - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 385), range = c(1, 12), breaks = c(10, 100, 200, 300), name = "ARG Count (CP60M)")
dev.off()

# Bubbleplots (Percentage)
taxa.g1.per.df <- data.frame(SIZE = c((colSums(taxa.g1.t0[1:10])/sum(colSums(taxa.g1.t0[1:10])))*100,
                                      (colSums(taxa.g1.t1[1:10])/sum(colSums(taxa.g1.t1[1:10])))*100,
                                      (colSums(taxa.g1.t2[1:10])/sum(colSums(taxa.g1.t2[1:10])))*100),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g1.t0)[1:10], 3), levels = sort(colnames(taxa.g1.t0), decreasing = TRUE))
                             )

taxa.g2.per.df <- data.frame(SIZE = c((colSums(taxa.g2.t0[1:10])/sum(colSums(taxa.g2.t0[1:10])))*100,
                                      (colSums(taxa.g2.t1[1:10])/sum(colSums(taxa.g2.t1[1:10])))*100,
                                      (colSums(taxa.g2.t2[1:10])/sum(colSums(taxa.g2.t2[1:10])))*100),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g2.t0)[1:10], 3), levels = sort(colnames(taxa.g2.t0), decreasing = TRUE))
                             )

taxa.g3.per.df <- data.frame(SIZE = c((colSums(taxa.g3.t0[1:10])/sum(colSums(taxa.g3.t0[1:10])))*100,
                                      (colSums(taxa.g3.t1[1:10])/sum(colSums(taxa.g3.t1[1:10])))*100,
                                      (colSums(taxa.g3.t2[1:10])/sum(colSums(taxa.g3.t2[1:10])))*100),
                             TIME = c(rep("t0", 10), rep("t1", 10), rep("t2", 10)),
                             SPECIES = factor(rep(colnames(taxa.g3.t0)[1:10], 3), levels = sort(colnames(taxa.g3.t0), decreasing = TRUE))
                             )

png("output/07_visualization/amr_bubble_families_rel_ssg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g1.per.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("SSG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 100), range = c(1, 12), breaks = c(10, 40, 80), name = "ARG Abundance (%)")
dev.off()

png("output/07_visualization/amr_bubble_families_rel_5dg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g2.per.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("5DG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 100), range = c(1, 12), breaks = c(10, 40, 80), name = "ARG Abundance (%)")
dev.off()

png("output/07_visualization/amr_bubble_families_rel_switched.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g3.per.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("SWITCHED - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 100), range = c(1, 12), breaks = c(10, 40, 80), name = "ARG Abundance (%)")
dev.off()

# --------------------------------------------------------------------------------------------------------

# [10] Run MicrobiomeExplorer (16s/WGS)

converted_biom <- readData(filepath = "output/02_taxonomic_profiling/kraken_biom/kraken2.biom", type = "BIOM")
saveRDS(converted_biom, "output/07_visualization/kraken2.rds")

#runMicrobiomeExplorer()