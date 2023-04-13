# --------------------------------------------------------------------------------------------------------
# Title: MetaGEN.R
# Author: Silver A. Wolf
# Last Modified: Thu, 13.04.2023
# Version: 0.6.7
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
library("microbiomeutilities")
library("randomcoloR")
library("openxlsx")
library("PathoStat")
library("phyloseq")
library("stringr")
library("taxonomizr")
library("tidyr")
library("tidyverse")
library("vegan")

# --------------------------------------------------------------------------------------------------------

# [01] Import and preprocess analysis data

# AMR Results
abricate <- read.csv("output/07_amr/abricate/amr/abricate.summary", sep = "\t")
abricate$X.FILE <- gsub(".tab", "", abricate$X.FILE)

# Virulence Results
vir.table <- read.csv("output/07_amr/abricate/vir/abricate.summary", sep = "\t")
vir.table$X.FILE <- gsub(".tab", "", vir.table$X.FILE)

# MegaRes database
megares.db <- read.csv("db/megares_filtered.csv")
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
rawCountTable <- read.table("output/07_amr/coverm/coverm.summary", header = TRUE, sep = "\t", row.names = 1)
colnames(rawCountTable) <- sub("output.07_amr.coverm.", "", colnames(rawCountTable))

# BIOM File
data.biom <- import_biom("output/02_taxonomic_profiling/kraken_biom/bracken_update.biom", parseFunction = parse_taxonomy_default)

# Metadata
meta.raw <- read.xlsx("metadata/23_03_Horses_Overview.xlsx", sheet = 1)
meta.raw <- meta.raw[,-c(1)]
colnames(meta.raw)[1] <- "SampleID"

# Group order
groups.order <- c("SSG", "5DG", "REF")
timepoints.order <- c("t0", "t1", "t2", "REF")

# --------------------------------------------------------------------------------------------------------

# [02] Diversity Estimations

# Sample Depth
sample_depth <- sort(sample_sums(data.biom))
min(sample_depth)
max(sample_depth)
median(sample_depth)
sum(sample_depth)

# Alpha Diversity (Raw)
data.alpha <- microbiome::alpha(data.biom)

# Plot Rarefaction Curve
#png("output/08_visualization/div_rarefaction_curve.png", width = 16, height = 16, units = "cm", res = 500)
#rarecurve(as(t(otu_table(data.biom)), "matrix"), step = 50, cex = 0.5, xlab = "Sequencing Depth (#Reads)", ylab = "Species (#Taxa)")
#dev.off()

# Alpha diversity (Rarefy)
data.rarefy <- rarefy_even_depth(data.biom, rngseed = 1, sample.size = min(sample_depth), replace = FALSE, trimOTUs = FALSE)
data.alpha.rarefy <- microbiome::alpha(data.rarefy)
#data.rarefy <- aggregate_top_taxa2(data.rarefy, 22, "Rank2")

# Beta Diversity (Bray-Curtis distance)
braycurtis <- phyloseq::distance(data.rarefy, method = "bray")
data.bray <- as.matrix(braycurtis)

# Beta Diversity (PCoA)
braycurtis.pcoa <- ordinate(physeq = data.rarefy, method = "PCoA", distance = "bray")
data.pcoa <- as.data.frame(braycurtis.pcoa$vectors, row.names = NULL, optional = FALSE, cut.names = FALSE, col.names = names(braycurtis.pcoa$vectors), fix.empty.names = TRUE, stringsAsFactors = default.stringsAsFactors())

# Add Metadata
meta.sorted = meta.raw[match(rownames(data.alpha), meta.raw$SampleID),]
data.alpha$HORSE = meta.sorted$HorseID
data.alpha$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.alpha$TIMEPOINT = factor(meta.sorted$Timepoint, levels = timepoints.order)
data.alpha.rarefy$HORSE = meta.sorted$HorseID
data.alpha.rarefy$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.alpha.rarefy$TIMEPOINT = factor(meta.sorted$Timepoint, levels = timepoints.order)
data.alpha.filtered <- data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "REF", ]
data.pcoa$HORSE = meta.sorted$HorseID
data.pcoa$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.pcoa$TIMEPOINT = factor(meta.sorted$Timepoint, levels = timepoints.order)
data.pcoa$TIME_GROUP = paste(meta.sorted$Timepoint, meta.sorted$AB_Group, sep = " ")

# Experimental
rownames(meta.sorted) <- meta.sorted$SampleID
data.biom@sam_data <- sample_data(meta.sorted)
data.rarefy@sam_data <- sample_data(meta.sorted)

# Export Diversities
write.csv(data.alpha, file = "output/08_visualization/tab_div_alpha_raw.csv", quote = FALSE)
write.csv(data.alpha.rarefy, file = "output/08_visualization/tab_div_alpha_rarefy.csv", quote = FALSE)
write.csv(data.bray, file = "output/08_visualization/tab_div_beta_distance.csv", quote = FALSE)
write.csv(data.pcoa, file = "output/08_visualization/tab_div_beta_pcoa.csv", quote = FALSE)

# Export human-readable OTU table
data.otu <- phyloseq_to_df(data.biom)
write.csv(data.otu, file = "output/08_visualization/tab_otu.csv", row.names = FALSE, quote = FALSE)
data.otu.rarefy <- phyloseq_to_df(data.rarefy)
write.csv(data.otu.rarefy, file = "output/08_visualization/tab_otu_rarefy.csv", row.names = FALSE, quote = FALSE)

# Specific taxa level counts
length(unique(data.otu$Rank6))

# PCA
colours.days = c("t0" = "#00BA38",
                 "t1" = "#F8766D",
                 "t2" = "#619CFF",
                 "REF" = "#606060"
                 )

colours.groups = c("SSG" = "#00ff7f",
                   "5DG" = "#ffa500",
                   "REF" = "#606060"
                   )

eigenvalue_pc1 = round(braycurtis.pcoa$values$Relative_eig[1]*100, 1)
eigenvalue_pc2 = round(braycurtis.pcoa$values$Relative_eig[2]*100, 1)

# Timepoints - All Samples
png("output/08_visualization/div_pca_time_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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
png("output/08_visualization/div_pca_time_ssg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "SSG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days[1:3],
          title = "SSG - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - 5DG
png("output/08_visualization/div_pca_time_5dg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "5DG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days[1:3],
          title = "5DG - Beta Diversity PCA (Bray–Curtis Dissimilarity)",
          ) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - All Samples
png("output/08_visualization/div_pca_group_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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
png("output/08_visualization/div_pca_group_t0.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t0", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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
png("output/08_visualization/div_pca_group_t1.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t1", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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
png("output/08_visualization/div_pca_group_t2.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t2", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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
png("output/08_visualization/div_pca_group_time_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.57, 0.17),
          ylim = c(-0.43, 0.46),
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

# Boxplots (Timewise)
boxplot.timepoints <- list(c("t0", "t1"), c("t1", "t2"), c("t0","t2"))

png("output/08_visualization/div_box_shan_time.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.filtered, aes(x = TIMEPOINT, y = diversity_shannon, fill = TIMEPOINT)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(1.9, 9)) +
        scale_y_continuous(breaks = c(3, 5, 7)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days[1:3]) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
                           method = "wilcox.test",
                           label.y = c(7.5, 8.1, 8.7),
                           size = 3,
                           paired = TRUE)
dev.off()

png("output/08_visualization/div_box_even_time.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.filtered, aes(x = TIMEPOINT, y = evenness_simpson, fill = TIMEPOINT)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(0.0007, 0.055)) +
        scale_y_continuous(breaks = c(0.01, 0.03, 0.05)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days[1:3]) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
                           method = "wilcox.test",
                           label.y = c(0.042, 0.047, 0.052),
                           size = 3,
                           paired = TRUE)
dev.off()

# Barplots
bar_data_aggregated <- aggregate_top_taxa2(data.rarefy, 9, "Rank2")
bar_data_melted <- psmelt(bar_data_aggregated)
num.taxa <- length(unique(bar_data_melted$OTU))
palette <- distinctColorPalette(num.taxa)

bar_ext_horse = c()
bar_ext_group = c()
bar_ext_time = c()
i = 1

for(e in bar_data_melted$Sample){
        current_line = data.alpha.rarefy[rownames(data.alpha.rarefy) == e,]
        bar_ext_horse[i] <- current_line$HORSE
        bar_ext_group[i] <- as.character(current_line$AB_GROUP)
        bar_ext_time[i] <- as.character(current_line$TIMEPOINT)
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

png("output/08_visualization/tax_bar_horses_5dg.png", width = 20, height = 15, units = "cm", res = 500)
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

png("output/08_visualization/tax_bar_horses_ssg.png", width = 15, height = 15, units = "cm", res = 500)
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

# Individual Horses - Reference
bar_data_ref <- bar_data_melted[bar_data_melted$AB_GROUP == "REF",]

png("output/08_visualization/tax_bar_horses_ref.png", width = 10, height = 15, units = "cm", res = 500)
ggplot(bar_data_ref, aes(fill = OTU, y = Abundance, x = HORSE)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Relative Abundance (REF)",
             x = "Horses",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
        ) +
        facet_grid(~ TIMEPOINT) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
dev.off()

# Summarized groups
png("output/08_visualization/tax_bar_sum_groups.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted[bar_data_melted$AB_GROUP != "REF",], aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
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

png("output/08_visualization/tax_bar_sum_timepoints.png", width = 20, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted[bar_data_melted$AB_GROUP != "REF",], aes(fill = OTU, y = Abundance, x = AB_GROUP)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Mean Abundance (Groups)",
             x = "AB Groups",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
        ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_grid(~ TIMEPOINT) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
dev.off()

png("output/08_visualization/tax_bar_sum_refs.png", width = 10, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted[bar_data_melted$AB_GROUP == "REF",], aes(fill = OTU, y = Abundance, x = AB_GROUP)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = palette) +
        labs(title = "Mean Abundance (Groups)",
             x = "AB Groups",
             y = "Relative Abundance (%)",
             fill = "Top 10 Phyla"
        ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_grid(~ TIMEPOINT) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
dev.off()

# Percentages of individual taxa
per_taxa_total = sum(bar_data_melted$Abundance)/length(sample_depth)

per_bact_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacteroidota", ]$Abundance)/length(sample_depth)
per_bact_norm = (per_bact_total/per_taxa_total) * 100

per_firm_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacillota", ]$Abundance)/length(sample_depth)
per_firm_norm = (per_firm_total/per_taxa_total) * 100

per_prot_total = sum(bar_data_melted[bar_data_melted$OTU == "Pseudomonadota", ]$Abundance)/length(sample_depth)
per_prot_norm = (per_prot_total/per_taxa_total) * 100

per_spir_total = sum(bar_data_melted[bar_data_melted$OTU == "Spirochaetota", ]$Abundance)/length(sample_depth)
per_spir_norm = (per_spir_total/per_taxa_total) * 100

per_verr_total = sum(bar_data_melted[bar_data_melted$OTU == "Verrucomicrobiota", ]$Abundance)/length(sample_depth)
per_verr_norm = (per_verr_total/per_taxa_total) * 100

per_bact_norm
per_firm_norm
per_prot_norm
per_spir_norm
per_verr_norm

# Abundance heatmap
abundance_columns <- unique(bar_data_melted$OTU)
abundance_rows <- c("SSG_t0", "5DG_t0", "SSG_t1", "5DG_t1", "SSG_t2", "5DG_t2", "REF")
abundance_matrix <- matrix(0, length(abundance_rows), length(unique(bar_data_melted$OTU)))

j = 1

for (c in abundance_columns){
        i = 1
        for (r in abundance_rows){
                if (r != "REF"){
                        group <- strsplit(r, "_")[[1]][1]
                        timepoint <- strsplit(r, "_")[[1]][2]
                        abundance_matrix[i, j] <- mean(bar_data_melted[bar_data_melted$TIMEPOINT == timepoint & bar_data_melted$AB_GROUP == group & bar_data_melted$OTU == c,]$Abundance)
                        i = i + 1 
                } else {
                        abundance_matrix[i, j] <- mean(bar_data_melted[bar_data_melted$AB_GROUP == r & bar_data_melted$OTU == c,]$Abundance)
                        i = i + 1 
                        }
                }
        j = j + 1
        }

abundance_matrix <- t(abundance_matrix)
rownames(abundance_matrix) <- abundance_columns

colnames(abundance_matrix) <- c("SSG (t0)", "5DG (t0)", "SSG (t1)", "5DG (t1)", "SSG (t2)", "5DG (t2)", "REF")

png("output/08_visualization/tax_heat_sim.png", width = 20, height = 15, units = "cm", res = 500)
Heatmap(log2(abundance_matrix),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        col = c("black", "darkred", "red", "orange", "yellow"),
        row_title = "Top 10 Phyla",
        name = "log2(Abundance)",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
        )
dev.off()

colnames(abundance_matrix) <- c("SSG", "5DG", "SSG", "5DG", "SSG", "5DG", "REF")

png("output/08_visualization/tax_heat_time.png", width = 20, height = 15, units = "cm", res = 500)
Heatmap(log2(abundance_matrix),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        col = c("black", "darkred", "red", "orange", "yellow"),
        column_split = c(rep("t0", 2), rep("t1", 2), rep("t2", 2), "REF"),
        row_title = "Top 10 Phyla",
        name = "log2(Abundance)",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
        )
dev.off()

# --------------------------------------------------------------------------------------------------------

# [03] AMR Heatmap (Assembly)

# Count individual resistance genes
abricate.count <- abricate
for (i in 1:length(sample_depth)){
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
abricate.meta <- meta.raw[match(rownames(abricate.matrix), meta.raw$SampleID),]
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
write.csv(amr.df.stats, file = "output/08_visualization/tab_amr_counts.csv", quote = FALSE)

# Prepare dataframe for heatmap
amr.norm.reads <- subset(amr.df.stats, select = c("SAMPLE", "CPM", "CP60M"))
amr.norm.reads$DIV = data.alpha.rarefy$diversity_shannon
amr.norm.reads$TIMEPOINT = data.alpha.rarefy$TIMEPOINT
amr.norm.reads$AB_GROUP = data.alpha.rarefy$AB_GROUP

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
colours.div = colorRamp2(c(1, 3, 6, 8), c("red", "orange", "yellow", "green"))

# Set Annotations for Heatmap
annot.column = HeatmapAnnotation(class = megares.db.clean$V2,
                                 col = list(class = colours.amr)
                                 )

annot.row.left = rowAnnotation(timepoint = abricate.meta$Timepoint,
                               ab_group = abricate.meta$AB_Group,
                               alpha_div = abricate.meta$DIV,
                               "#amr_genes" = amr.counts,
                               col = list(ab_group = colours.groups,
                                          timepoint = colours.days,
                                          alpha_div = colours.div,
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
png("output/08_visualization/amr_heat_abricate.png", width = 40, height = 20, units = "cm", res = 500)
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

# [04] Statistical Correlation (Assembly)

# Correlation - AMR and Diversity
cor.test(amr.norm.reads$CP60M, amr.norm.reads$DIV, method = "spearman")

png("output/08_visualization/amr_div_cor.png", width = 17, height = 16, units = "cm", res = 500)
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

# Significant Differences between t0 and t1 -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t2" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t0" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t1" & data.alpha.rarefy$AB_GROUP == "SSG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Diversity - 5DG

# Significant Differences between t0 and t1 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t2" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t0" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT != "t1" & data.alpha.rarefy$AB_GROUP == "5DG",]
stat.df <- data.frame(TIME = stat.data$TIMEPOINT, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$TIME, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Diversity - SSG vs. 5DG

# Significant Differences at t0-t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t0 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t1 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t2 -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t2" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, DIV = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$DIV), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# AMR - SSG

# Significant Differences between t0 and t1 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t2" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t0" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t1" & amr.norm.reads$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# AMR - 5DG

# Significant Differences between t0 and t1 -> yes
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t2" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t1 and t2 -> yes
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t0" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# Significant Differences between t0 and t2 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT != "t1" & amr.norm.reads$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH", paired = TRUE)
stat.res

# AMR - SSG vs. 5DG

# Significant Differences at t0-t2 -> no
stat.data <- amr.norm.reads[data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t0 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t1 -> no
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences at t2 -> yes
stat.data <- amr.norm.reads[amr.norm.reads$TIMEPOINT == "t2" & data.alpha.rarefy$AB_GROUP != "REF",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$CP60M)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# --------------------------------------------------------------------------------------------------------

# [05] AMR Normalization (Reads)

# [05.01] AMR Gene Heatmap

# Read raw data
Gene_sumCountTable <- aggregate(. ~ Resistance_Symbol, rawCountTable[,-c(1:3)], sum)
rownames(Gene_sumCountTable) <- Gene_sumCountTable$Resistance_Symbol
colnames(Gene_sumCountTable) <- gsub("output.06_amr.coverm.", "", colnames(Gene_sumCountTable))
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
write.csv(Gene_normCounts, file = "output/08_visualization/tab_amr_tmm.csv", quote = FALSE)

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
png("output/08_visualization/amr_heat_coverm_genes.png", width = 40, height = 20, units = "cm", res = 500)
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
                                        alpha_div = abricate.meta$DIV,
                                        col = list(ab_group = colours.groups,
                                                   timepoint = colours.days,
                                                   alpha_div = colours.div)
                                        ),
        
        right_annotation = rowAnnotation(horse = anno_text(abricate.meta$HorseID,
                                                           gp = gpar(fontsize = 8))),
        
        top_annotation = columnAnnotation(amr_class = Gene_AMR_DF$AMR_Class)
        )
dev.off()

# [05.02] AMR Class Heatmap

# Read raw data
Class_sumCountTable <- aggregate(. ~ Resistance_Class, rawCountTable[,-c(1,3,4)], sum)
rownames(Class_sumCountTable) <- Class_sumCountTable$Resistance_Class
colnames(Class_sumCountTable) <- gsub("output.06_amr.coverm.", "", colnames(Class_sumCountTable))
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
write.csv(Class_normCounts, file = "output/08_visualization/tab_amr_tmm_classes.csv", quote = FALSE)

# Bin TMM values by counts
Class_groupedCounts = Class_normCounts
Class_groupedCounts[Class_groupedCounts < 500] = 0
Class_groupedCounts[Class_groupedCounts >= 500 & Class_groupedCounts < 1000] = 1
Class_groupedCounts[Class_groupedCounts >= 1000 & Class_groupedCounts < 5000] = 2
Class_groupedCounts[Class_groupedCounts >= 5000] = 3
Class_groupedCounts <- t(Class_groupedCounts)
colnames(Class_groupedCounts) <- tolower(colnames(Class_groupedCounts))

# Plot AMR Class Heatmap
png("output/08_visualization/amr_heat_coverm_classes.png", width = 20, height = 25, units = "cm", res = 500)
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
                                        alpha_div = abricate.meta$DIV,
                                        col = list(ab_group = colours.groups,
                                                   timepoint = colours.days,
                                                   alpha_div = colours.div)
                                        ),
        
        right_annotation = rowAnnotation(horse = anno_text(abricate.meta$HorseID,
                                         gp = gpar(fontsize = 8)))
        )
dev.off()

# [05.03] Scatterplots / Boxplots of AMR Classes
e = c()
f = c()
j = 1

for (i in rownames(data.alpha.rarefy)) {
        e[j] <- sum(Class_normCounts[,i])
        f[j] <- i
        j = j + 1
        }

Class_AMR_SUM <- data.frame(SAMPLE = f, AMR = e, TIMEPOINT = data.alpha.rarefy$TIMEPOINT, AB_GROUP = data.alpha.rarefy$AB_GROUP, HORSE = data.alpha.rarefy$HORSE)

png("output/08_visualization/amr_sum_scatter.png", width = 30, height = 20, units = "cm", res = 500)
ggplot(Class_AMR_SUM[Class_AMR_SUM$AB_GROUP != "REF",], aes(x = TIMEPOINT, y = log2(AMR))) +
        coord_cartesian(ylim = c(17, 24)) +
        geom_line(aes(group = HORSE)) +
        geom_point() +
        geom_point(size = 2) +
        facet_wrap(~AB_GROUP + HORSE, scale = "free")
dev.off()

# Total AMR Sum
png("output/08_visualization/amr_sum_box.png", width = 15, height = 10, units = "cm", res = 500)
ggplot(Class_AMR_SUM[Class_AMR_SUM$AB_GROUP != "REF",], aes(x = TIMEPOINT, y = log2(AMR), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(18, 25)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AB_GROUP, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
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
num.ssg <- length(meta.raw[meta.raw$AB_Group == "SSG",]$SampleID)
amr.class.df.ssg <- data.frame(AMR_Class = rep(colnames(amr.class.df.filtered[1:amr.class.length]), each = num.ssg),
                               AMR_TMM = unlist(amr.class.df.filtered[1:amr.class.length]),
                               AB_GROUP = rep(amr.class.df.filtered$AB_GROUP, amr.class.length),
                               HORSE = rep(amr.class.df.filtered$HORSE, amr.class.length),
                               TIMEPOINT = rep(amr.class.df.filtered$TIMEPOINT, amr.class.length),
                               TIME_GROUP = rep(amr.class.df.filtered$TIME_GROUP, amr.class.length)
                               )

png("output/08_visualization/amr_box_group_ssg.png", width = 40, height = 30, units = "cm", res = 500)
ggplot(amr.class.df.ssg, aes(x = TIMEPOINT, y = log2(AMR_TMM + 1), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(0, 30)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AMR_Class, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days[1:3]) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
                           method = "wilcox.test",
                           label.y = c(23, 25.5, 28),
                           size = 2,
                           paired = TRUE,
                           method.args = list(exact = FALSE))
dev.off()

# 5DG
amr.class.df.filtered <- amr.class.df[amr.class.df$AB_GROUP == "5DG", ]
amr.class.length <- length(colnames(amr.class.df.filtered)) - 4
num.5dg <- length(meta.raw[meta.raw$AB_Group == "5DG",]$SampleID)
amr.class.df.5dg <- data.frame(AMR_Class = rep(colnames(amr.class.df.filtered[1:amr.class.length]), each = num.5dg),
                               AMR_TMM = unlist(amr.class.df.filtered[1:amr.class.length]),
                               AB_GROUP = rep(amr.class.df.filtered$AB_GROUP, amr.class.length),
                               HORSE = rep(amr.class.df.filtered$HORSE, amr.class.length),
                               TIMEPOINT = rep(amr.class.df.filtered$TIMEPOINT, amr.class.length),
                               TIME_GROUP = rep(amr.class.df.filtered$TIME_GROUP, amr.class.length)
                               )

png("output/08_visualization/amr_box_group_5dg.png", width = 40, height = 30, units = "cm", res = 500)
ggplot(amr.class.df.5dg, aes(x = TIMEPOINT, y = log2(AMR_TMM + 1), fill = TIMEPOINT)) +
        coord_cartesian(ylim = c(0, 30)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AMR_Class, scale = "free") +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days[1:3]) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
                           method = "wilcox.test",
                           label.y = c(23, 25.5, 28),
                           size = 2,
                           paired = TRUE,
                           method.args = list(exact = FALSE))
dev.off()

# [05.04] Statistical Correlation (Reads)
Class_AMR_SUM$AMR_LOG2 <- log2(Class_AMR_SUM$AMR)
Class_AMR_SUM$DIV = amr.norm.reads$DIV

cor.test(Class_AMR_SUM$AMR, Class_AMR_SUM$DIV, method = "spearman")

png("output/08_visualization/amr_div_cor_reads.png", width = 17, height = 16, units = "cm", res = 500)
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
          size = 2.5,
          cor.coef.coord = c(22.5, 4)
          ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(shape = "AB GROUP")
dev.off()

# --------------------------------------------------------------------------------------------------------

# [06] Differential Taxa Analysis

# Family Level
edger.aggretated <- aggregate_top_taxa2(data.biom, 10^5, "Rank5")

# Prepare data for edgeR analysis
edger.groups <- data.pcoa[with(data.pcoa, order(rownames(data.pcoa))), ]$TIME_GROUP
edger.dge <- phyloseq_to_edgeR(edger.aggretated, group = edger.groups)

# C1 - SSG: t0 vs. t1
edger.res.c1 <- exactTest(edger.dge, pair = c("t0 SSG", "t1 SSG"))
edger.top.c1 <- topTags(edger.res.c1, n = nrow(edger.res.c1$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c1 <- edger.top.c1$table[edger.top.c1$table$FDR < 0.05 & abs(edger.top.c1$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_ssg_t0_t1.csv", quote = FALSE)

# C2 - SSG: t1 vs. t2
edger.res.c2 <- exactTest(edger.dge, pair = c("t1 SSG", "t2 SSG"))
edger.top.c2 <- topTags(edger.res.c2, n = nrow(edger.res.c2$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c2 <- edger.top.c2$table[edger.top.c2$table$FDR < 0.05 & abs(edger.top.c2$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_ssg_t1_t2.csv", quote = FALSE)

# C3 - SSG: t0 vs. t2
edger.res.c3 <- exactTest(edger.dge, pair = c("t0 SSG", "t2 SSG"))
edger.top.c3 <- topTags(edger.res.c3, n = nrow(edger.res.c3$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c3 <- edger.top.c3$table[edger.top.c3$table$FDR < 0.05 & abs(edger.top.c3$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_ssg_t0_t2.csv", quote = FALSE)

# C4 - 5DG: t0 vs. t1
edger.res.c4 <- exactTest(edger.dge, pair = c("t0 5DG", "t1 5DG"))
edger.top.c4 <- topTags(edger.res.c4, n = nrow(edger.res.c4$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c4 <- edger.top.c4$table[edger.top.c4$table$FDR < 0.05 & abs(edger.top.c4$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_5dg_t0_t1.csv", quote = FALSE)

# C5 - 5DG: t1 vs. t2
edger.res.c5 <- exactTest(edger.dge, pair = c("t1 5DG", "t2 5DG"))
edger.top.c5 <- topTags(edger.res.c5, n = nrow(edger.res.c5$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c5 <- edger.top.c5$table[edger.top.c5$table$FDR < 0.05 & abs(edger.top.c5$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_5dg_t1_t2.csv", quote = FALSE)

# C6 - 5DG: t0 vs. t2
edger.res.c6 <- exactTest(edger.dge, pair = c("t0 5DG", "t2 5DG"))
edger.top.c6 <- topTags(edger.res.c6, n = nrow(edger.res.c6$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c6 <- edger.top.c6$table[edger.top.c6$table$FDR < 0.05 & abs(edger.top.c6$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_5dg_t0_t2.csv", quote = FALSE)

# C7 - t0: SSG vs. 5DG
edger.res.c7 <- exactTest(edger.dge, pair = c("t0 SSG", "t0 5DG"))
edger.top.c7 <- topTags(edger.res.c7, n = nrow(edger.res.c7$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c7 <- edger.top.c7$table[edger.top.c7$table$FDR < 0.05 & abs(edger.top.c7$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_t0_ssg_5dg.csv", quote = FALSE)

# C8 - t1: SSG vs. 5DG
edger.res.c8 <- exactTest(edger.dge, pair = c("t1 SSG", "t1 5DG"))
edger.top.c8 <- topTags(edger.res.c8, n = nrow(edger.res.c8$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c8 <- edger.top.c8$table[edger.top.c8$table$FDR < 0.05 & abs(edger.top.c8$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_t1_ssg_5dg.csv", quote = FALSE)

# C9 - t2: SSG vs. 5DG
edger.res.c9 <- exactTest(edger.dge, pair = c("t2 SSG", "t2 5DG"))
edger.top.c9 <- topTags(edger.res.c9, n = nrow(edger.res.c9$table), adjust.method = "BH", sort.by = "logFC")
edger.sel.c9 <- edger.top.c9$table[edger.top.c9$table$FDR < 0.05 & abs(edger.top.c9$table$logFC) > 1, ]
write.csv(edger.top.c1$table, file = "output/08_visualization/tab_edger_t2_ssg_5dg.csv", quote = FALSE)

# --------------------------------------------------------------------------------------------------------

# [07] Boxplots for abundance of specific families

boxplot.families.biom <- aggregate_top_taxa2(data.rarefy, 22, "Rank5")
boxplot.families.melted <- psmelt(boxplot.families.biom)

box.ext.horse = c()
box.ext.group = c()
box.ext.time = c()
j = 1

for(f in boxplot.families.melted$Sample){
        current.line = data.alpha.rarefy[rownames(data.alpha.rarefy) == f,]
        box.ext.horse[j] <- current.line$HORSE
        box.ext.group[j] <- as.character(current.line$AB_GROUP)
        box.ext.time[j] <- as.character(current.line$TIMEPOINT)
        j = j + 1
}

boxplot.families.melted$HORSE <- box.ext.horse
boxplot.families.melted$AB_GROUP <- factor(box.ext.group, levels = groups.order)
boxplot.families.melted$TIMEPOINT <- box.ext.time

G2 <- boxplot.families.melted[boxplot.families.melted$Rank5 == "f__Enterobacteriaceae", ]
G2 <- G2[G2$AB_GROUP != "REF",]

png("output/08_visualization/tax_box_f_enterobacteriaceae_time.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(G2, aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
        geom_boxplot(alpha = 0.9) +
        geom_jitter(alpha = 0.5) +
        facet_wrap(~AB_GROUP, scale = "free") +
        coord_cartesian(ylim = c(11, 19.5)) +
        scale_y_continuous(breaks = c(9, 12, 15, 18)) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        scale_fill_manual(values = colours.days[1:3]) +
        stat_compare_means(comparisons = boxplot.timepoints,
                           alternative = "two.sided",
                           method = "wilcox.test",
                           label.y = c(17, 18, 19),
                           size = 3,
                           paired = TRUE,
                           method.args = list(exact = FALSE)) +
        ggtitle("Abundance - Enterobacteriaceae (Family)") +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

# --------------------------------------------------------------------------------------------------------

# [08] Calculate differences between timepoints

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
arg.5dg.t0 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "5DG" & amr.norm.reads$TIMEPOINT == "t0", ]$CP60M)
arg.5dg.t1 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "5DG" & amr.norm.reads$TIMEPOINT == "t1", ]$CP60M)
arg.5dg.t2 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "5DG" & amr.norm.reads$TIMEPOINT == "t2", ]$CP60M)

arg.5dg.diff.t0.t1 = round(arg.5dg.t1 - arg.5dg.t0, 2)
arg.5dg.diff.t0.t2 = round(arg.5dg.t2 - arg.5dg.t1, 2)

arg.ssg.t0 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "SSG" & amr.norm.reads$TIMEPOINT == "t0", ]$CP60M)
arg.ssg.t1 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "SSG" & amr.norm.reads$TIMEPOINT == "t1", ]$CP60M)
arg.ssg.t2 = mean(amr.norm.reads[amr.norm.reads$AB_GROUP == "SSG" & amr.norm.reads$TIMEPOINT == "t2", ]$CP60M)

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

# [09] Visualizing ARG abundance per family

# Prepare input data
kreport = read.csv("output/07_amr/abricate/amr/kraken2.summary", header = TRUE, sep = "\t", row.names = 1)
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

# Failsafe to secure correct sample list (can be removed)
taxa.matrix <- taxa.matrix[rownames(taxa.matrix) %in% rownames(data.alpha.rarefy),]

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

# Bubbleplots (Absolute)
taxa.g1.t0.nom <- taxa.g1.t0
taxa.g1.t1.nom <- taxa.g1.t1
taxa.g1.t2.nom <- taxa.g1.t2

taxa.g2.t0.nom <- taxa.g2.t0
taxa.g2.t1.nom <- taxa.g2.t1
taxa.g2.t2.nom <- taxa.g2.t2

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

png("output/08_visualization/amr_bubble_families_abs_ssg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g1.abs.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("SSG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 465), range = c(1, 12), breaks = c(10, 100, 200, 300), name = "ARG Count (CP60M)")
dev.off()

png("output/08_visualization/amr_bubble_families_abs_5dg.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(taxa.g2.abs.df, aes(x = TIME, y = SPECIES, size = SIZE)) +
        geom_point(alpha = 0.9, color = "red") +
        xlab("TIMEPOINT") +
        ylab("FAMILY") +
        ggtitle("5DG - FAMILIES") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 10)) +
        scale_size_continuous(limits = c(0, 465), range = c(1, 12), breaks = c(10, 100, 200, 300), name = "ARG Count (CP60M)")
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

png("output/08_visualization/amr_bubble_families_rel_ssg.png", width = 20, height = 10, units = "cm", res = 500)
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

png("output/08_visualization/amr_bubble_families_rel_5dg.png", width = 20, height = 10, units = "cm", res = 500)
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

# AMR table per taxa
taxnames.all <- as.data.frame(getTaxonomy(taxids, sqlFile = "tmp/accessionTaxa.sql", desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")))
taxnames.amr <- data.frame(taxnames.all, t(kreport))
write.csv(taxnames.amr, file = "output/08_visualization/tab_amr_taxa.csv", quote = FALSE, row.names = TRUE)

# --------------------------------------------------------------------------------------------------------

# [10] Run MicrobiomeExplorer (16s/WGS)

converted_biom <- readData(filepath = "output/02_taxonomic_profiling/kraken_biom/bracken_update.biom", type = "BIOM")
saveRDS(converted_biom, "output/08_visualization/bracken.rds")

#runMicrobiomeExplorer()

# --------------------------------------------------------------------------------------------------------

# [11] Calculation of specific taxa abundances

samples.ssg.t0 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t0", ]$SampleID
samples.ssg.t1 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t1", ]$SampleID
samples.ssg.t2 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t2", ]$SampleID

samples.5dg.t0 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t0", ]$SampleID
samples.5dg.t1 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t1", ]$SampleID
samples.5dg.t2 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t2", ]$SampleID

samples.ref <- meta.sorted[meta.sorted$AB_Group == "REF", ]$SampleID

samples.taxa <- c("c__Subdivision5",
                  "g__Allisonella",
                  "g__Escherichia/Shigella",
                  "g__Fusobacterium",
                  "g__Lactobacillus",
                  "g__Phascolarctobacterium",
                  "g__Rhodococcus",
                  "g__Roseburia",
                  "g__Ruminococcus",
                  "g__Salmonella",
                  "g__Staphylococcus",
                  "g__Streptococcus",
                  "f__Acidaminococcaceae",
                  "f__Bacteroidaceae",
                  "f__Enterobacteriaceae",
                  "f__Lachnospiraceae",
                  "f__Moraxellaceae",
                  "f__Planococcaceae",
                  "f__Ruminococcaceae",
                  "f__Veillonellaceae"
                  )

c1 = c()
c2 = c()
c3 = c()
c4 = c()
c5 = c()
c6 = c()
c7 = c()
x = 1

for (t in samples.taxa){
        t1 <- data.otu.rarefy[data.otu.rarefy$Rank1 == t |
                              data.otu.rarefy$Rank2 == t |
                              data.otu.rarefy$Rank3 == t |
                              data.otu.rarefy$Rank4 == t |
                              data.otu.rarefy$Rank5 == t |
                              data.otu.rarefy$Rank6 == t |
                              data.otu.rarefy$Rank7 == t, ]
        
        t2 <- t1[samples.ssg.t0]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.ssg.t0)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c1[x] <- t6
        
        t2 <- t1[samples.ssg.t1]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.ssg.t1)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c2[x] <- t6
        
        t2 <- t1[samples.ssg.t2]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.ssg.t2)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c3[x] <- t6
        
        t2 <- t1[samples.5dg.t0]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.5dg.t0)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c4[x] <- t6
        
        t2 <- t1[samples.5dg.t1]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.5dg.t1)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c5[x] <- t6
        
        t2 <- t1[samples.5dg.t2]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.5dg.t2)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c6[x] <- t6
        
        t2 <- t1[samples.ref]
        t3 <- colSums(t2)
        t4 <- sum(t3)/length(samples.ref)
        t5 <- (t4/min(sample_depth))*100
        t6 <- round(t5, 2)
        c7[x] <- t6
        
        x = x + 1
}

taxa.abundancies <- data.frame(TAXA = samples.taxa,
                               SSG_t0 = c1, SSG_t1 = c2, SSG_t2 = c3,
                               FDG_t0 = c4, FDG_t1 = c5, FDG_t2 = c6,
                               REF = c7
                               )

write.csv(taxa.abundancies, file = "output/08_visualization/tab_otu_species_mean.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------------------------------------------------

# [12] Assessing differences in beta diversities

beta.horses <- unique(meta.sorted[!(meta.sorted$AB_Group %in% c("REF")),]$HorseID)
beta.excluded <- unique(meta.sorted[meta.sorted$AB_Group %in% c("REF"),]$SampleID)
beta.intra <- c()
beta.inter <- c()

for (h in beta.horses){
        t <- meta.sorted[meta.sorted$HorseID == h, ]$SampleID
        r <- data.bray[t[1], ]
        # We select t0 as a starting point and query its ID to receive t1 and t2 distances
        beta.intra <- c(beta.intra, r[(names(r) %in% t)])
        # Afterwards, we assess the distances of this t0 sample to all other samples
        beta.inter <- c(beta.inter, r[!(names(r) %in% t) & !(names(r) %in% beta.excluded)])
        }

beta.intra <- beta.intra[!(beta.intra %in% c(0))]
beta.df <- data.frame(GROUP = c(rep("Intra", length(beta.intra)), rep("Inter", length(beta.inter))), VAL = c(beta.intra, beta.inter))
beta.res <- pairwise.wilcox.test(beta.df$VAL, beta.df$GROUP, paired = FALSE, alternative = "greater", p.adjust.method = "BH")

# --------------------------------------------------------------------------------------------------------

# [13] Virulence Heatmap (Assembly)

# Count individual virulence genes
vir.count <- vir.table
for (i in 1:length(sample_depth)){
        for (j in 3:ncol(vir.table)){
                d = vir.table[i,j]
                if (d == "."){
                        k = 0
                } else{
                        k = str_count(d,"\\.")
                }
                vir.count[i,j] <- k
        }
}

vir.filtered <- vir.count[,-c(1,2)]
rownames(vir.filtered) <- vir.table[,1]
vir.matrix <- data.matrix(vir.filtered) - 1
vir.counts <- rowSums(vir.matrix)
vir.meta <- meta.raw[match(rownames(vir.matrix), meta.raw$SampleID),]
vir.meta$VIR_FOUND <- vir.counts
vir.meta$DIV = data.alpha.rarefy$diversity_shannon

# Update virulence counts
vir.df.counts <- data.frame(ids = rownames(vir.table), counts = vir.counts)
vir.df.stats <- merge(x = seqreport.filtered, y = vir.df.counts, by = "row.names")
vir.df.stats <- subset(vir.df.stats, select = c("SAMPLE", "X.Seq", "counts"))

# Read virulence Normalization by #Reads
vir.df.stats$X.Seq <- vir.df.stats$X.Seq*2
vir.df.stats$Norm <- vir.df.stats$counts/vir.df.stats$X.Seq
vir.df.stats$CPM <- vir.df.stats$counts/(vir.df.stats$X.Seq/1000000)
vir.df.stats$CP60M <- vir.df.stats$CPM*60
colnames(vir.df.stats) <- c("SAMPLE", "#PE_READS", "VIR_COUNTS", "NORM_COUNT", "CPM", "CP60M")
write.csv(vir.df.stats, file = "output/08_visualization/tab_vir_counts.csv", quote = FALSE)

# Prepare dataframe for heatmap
vir.norm.reads <- subset(vir.df.stats, select = c("SAMPLE", "CPM", "CP60M"))
vir.norm.reads$DIV = data.alpha.rarefy$diversity_shannon
vir.norm.reads$TIMEPOINT = data.alpha.rarefy$TIMEPOINT
vir.norm.reads$AB_GROUP = data.alpha.rarefy$AB_GROUP

# Set Annotations for Heatmap
annot.row.left = rowAnnotation(timepoint = vir.meta$Timepoint,
                               ab_group = vir.meta$AB_Group,
                               alpha_div = vir.meta$DIV,
                               "#vir_genes" = vir.meta$VIR_FOUND,
                               col = list(timepoint = colours.days,
                                          ab_group = colours.groups,
                                          alpha_div = colours.div,
                                          "#vir_genes" = colours.genes
                                          )
                               )

annot.row.right = rowAnnotation(horse = anno_text(vir.meta$HorseID,
                                                  gp = gpar(fontsize = 8)
                                                  )
                                )

re.order.rows <- abricate.meta[with(abricate.meta, order(Day, AB_Group, HorseID)), ]

# Group all counts > 1 into a single class
vir.matrix[vir.matrix >= 2] = 2

# Plot Heatmap
png("output/08_visualization/vir_heat_abricate.png", width = 40, height = 20, units = "cm", res = 500)
Heatmap(vir.matrix,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        #col = c("black", "red4", "firebrick3", "firebrick2", "red1"),
        col = c("black", "gold", "red"),
        
        column_names_gp = gpar(fontsize = 5),
        column_title = "Virulence Genes",
        
        row_order = re.order.rows$SampleID,
        row_split = abricate.meta$Timepoint,
        row_title = "Metagenome Samples",
        show_row_names = FALSE,
        
        heatmap_legend_param = list(
                title = "#copies", at = c(0, 1, 2),
                labels = c("0", "1", "\u2265 2")),
        
        left_annotation = annot.row.left,
        right_annotation = annot.row.right
        )
dev.off()

# Correlation - Virulence and Diversity
cor.test(vir.norm.reads$CP60M, vir.norm.reads$DIV, method = "spearman", exact = FALSE)

png("output/08_visualization/vir_div_cor.png", width = 17, height = 16, units = "cm", res = 500)
ggscatter(vir.norm.reads,
          x = "CP60M",
          y = "DIV",
          color = "TIMEPOINT",
          shape = "AB_GROUP",
          xlab = "VIR PER 60M READS",
          ylab = "SHANNON INDEX",
          add = "reg.line",
          add.params = list(color = "black",
                            fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman",
                                label.x = 100,
                                label.sep = "\n"),
          title = "CORRELATION - VIR / DIVERSITY - ASSEMBLY",
          palette = colours.days,
          size = 2.5
          ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(shape = "AB GROUP")
dev.off()

# --------------------------------------------------------------------------------------------------------

# [14] Plasmid Analysis

plas.counts.args <- read.csv("output/08_visualization/tab_mobile_arg_counts.tsv", sep = "\t")
plas.counts.vir <- read.csv("output/08_visualization/tab_mobile_vir_counts.tsv", sep = "\t")
colnames(plas.counts.args) <- c("SampleID", "mobile_args")
colnames(plas.counts.vir) <- c("SampleID", "mobile_virs")
plas.df <- list(meta.sorted, plas.counts.args, plas.counts.vir) %>% reduce(full_join, by = "SampleID")

# Barchart of args per timepoint
png("output/08_visualization/plas_box_time_args.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(plas.df[plas.df$AB_Group != "REF",], aes(x = Timepoint, y = log(mobile_args + 1), fill = Timepoint)) +
  geom_boxplot(alpha = 0.9) +
  geom_jitter(alpha = 0.5) +
  facet_wrap(~AB_Group, scale = "free") +
  coord_cartesian(ylim = c(0, 4.9)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days[1:3]) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(4.1, 4.3, 4.6),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE))
dev.off()

# Barchart of virs per timepoint
png("output/08_visualization/plas_box_time_vir.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(plas.df[plas.df$AB_Group != "REF",], aes(x = Timepoint, y = log(mobile_virs + 1), fill = Timepoint)) +
  geom_boxplot(alpha = 0.9) +
  geom_jitter(alpha = 0.5) +
  facet_wrap(~AB_Group, scale = "free") +
  coord_cartesian(ylim = c(0, 4.9)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days[1:3]) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(4.1, 4.3, 4.6),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE))
dev.off()

# Piecharts for ARGs per family
plas.species.args <- read.csv("output/08_visualization/tab_mobile_arg_list.tsv", sep = "\t")
plas.species.args.taxa <- as.data.frame(getTaxonomy(plas.species.args$taxid, sqlFile = "tmp/accessionTaxa.sql", desiredTaxa = "family"))

# Excludes taxa at higher ranks
plas.species.args.sum <- table(plas.species.args.taxa)

plas.species.args.sum <- as.data.frame(plas.species.args.sum)
plas.species.args.sum <- plas.species.args.sum[order(plas.species.args.sum$Freq, decreasing = TRUE), ]
plas.species.args.top <- plas.species.args.sum[1:4, ]
plas.species.args.other <- c("Other", sum(plas.species.args.sum[5:nrow(plas.species.args.sum), ]$Freq))
levels(plas.species.args.top$family) <- c(levels(plas.species.args.top$family), "Other")
plas.species.args.top <- rbind(plas.species.args.top, plas.species.args.other)
plas.species.args.top$Freq <- as.numeric(plas.species.args.top$Freq)
colnames(plas.species.args.top) <- c("Family", "Freq")
plas.species.args.top[is.na(plas.species.args.top)] <- 0
plas.species.args.percentages <- round(((plas.species.args.top$Freq/sum(plas.species.args.top$Freq))*100), 0)
plas.species.args.percentages <- paste(plas.species.args.percentages, "%", sep = "")

png("output/08_visualization/plas_piechart_args.png", width = 15, height = 12, units = "cm", res = 500)
ggplot(plas.species.args.top, aes(x = "", y = Freq, fill = Family)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = plas.species.args.percentages),
            position = position_stack(vjust = 0.53)) +
  theme_void() +
  ggtitle("\nMobile ARGs per Family") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Piecharts for vir per family
plas.species.vir <- read.csv("output/08_visualization/tab_mobile_vir_list.tsv", sep = "\t")
plas.species.vir.taxa <- as.data.frame(getTaxonomy(plas.species.vir$taxid, sqlFile = "tmp/accessionTaxa.sql", desiredTaxa = "family"))

# Excludes taxa at higher ranks
plas.species.vir.sum <- table(plas.species.vir.taxa)

plas.species.vir.sum <- as.data.frame(plas.species.vir.sum)
plas.species.vir.sum <- plas.species.vir.sum[order(plas.species.vir.sum$Freq, decreasing = TRUE), ]
plas.species.vir.top <- plas.species.vir.sum[1, ]
levels(plas.species.vir.top$family) <- c(levels(plas.species.vir.top$family), "Other")
plas.species.vir.top$Freq <- as.numeric(plas.species.vir.top$Freq)
colnames(plas.species.vir.top) <- c("Family", "Freq")
plas.species.vir.top[is.na(plas.species.vir.top)] <- 0
plas.species.vir.percentages <- round(((plas.species.vir.top$Freq/sum(plas.species.vir.top$Freq))*100), 0)
plas.species.vir.percentages <- paste(plas.species.vir.percentages, "%", sep = "")

png("output/08_visualization/plas_piechart_vir.png", width = 15, height = 12, units = "cm", res = 500)
ggplot(plas.species.vir.top[1,], aes(x = "", y = Freq, fill = Family)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = plas.species.vir.percentages),
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  ggtitle("\nMobile VFs per Family") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# --------------------------------------------------------------------------------------------------------

# [15] Analysis of MDR MAGs

amr.mags.raw <- read.csv("output/07_amr/abricate/mags/amr.tab", sep = "\t")
amr.mags.split <- data.frame(do.call("rbind", strsplit(as.character(amr.mags.raw$PRODUCT), ":", fixed = TRUE)))
amr.mags.merged <- data.frame(amr.mags.raw, amr.mags.split)

s = c()
c = c()
i = 1

for (sample in unique(amr.mags.merged$X.FILE)){
  s[i] <- sample
  c[i] <- length(unique(amr.mags.merged[amr.mags.merged$X.FILE == sample & amr.mags.merged$X1 == "Drugs", ]$X2))
  i = i + 1
}

amr.mags.sum <- data.frame(SAMPLE = s, DRUG_RESISTANCE_CLASSES = c)
amr.mags.count = 8265

amr.mags.overview.rownames <- c("Total MAGs", "MDR (> 2 antibiotic classes)", "Resistant (1-2 antibiotic classes)", "Susceptible (0 antibiotic classes)")
amr.mags.overview.colnames <- c("Count", "Percentage of MAGs")
amr.mags.overview.counts <- c(amr.mags.count,
                              nrow(amr.mags.sum[amr.mags.sum$DRUG_RESISTANCE_CLASSES > 2, ]),
                              nrow(amr.mags.sum[amr.mags.sum$DRUG_RESISTANCE_CLASSES == 2 | amr.mags.sum$DRUG_RESISTANCE_CLASSES == 1, ]),
                              amr.mags.count - nrow(amr.mags.sum[amr.mags.sum$DRUG_RESISTANCE_CLASSES > 0, ])
                              )
amr.mags.overview.percentages <- round((amr.mags.overview.counts/amr.mags.count) * 100, 1)

amr.mags.overview <- data.frame(Counts = amr.mags.overview.counts, Percentages = amr.mags.overview.percentages)
colnames(amr.mags.overview) <- amr.mags.overview.colnames
rownames(amr.mags.overview) <- amr.mags.overview.rownames

# Filter any hits below 1 before exporting
amr.mags.sum.filtered <- amr.mags.sum[amr.mags.sum$DRUG_RESISTANCE_CLASSES > 0, ]

# Add taxonomic assignment from gtdbtk
gtdbtk.ar <- read.csv("output/05_genomic_bins/gtdbtk_full/classify/gtdbtk.ar53.summary.tsv", sep = "\t")
gtdbtk.bac <- read.csv("output/05_genomic_bins/gtdbtk_full/classify/gtdbtk.bac120.summary.tsv", sep = "\t")
gtdbtk.all <- rbind(gtdbtk.ar, gtdbtk.bac)
gtdbtk.all.split <- separate(data = gtdbtk.all, col = classification, into = c("DOMAIN", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"), sep = ";")

gtdbtk.all.filter <- gtdbtk.all.split[, 1:8]
gtdbtk.all.filter$user_genome <- paste(gtdbtk.all.filter$user_genome, ".fa", sep = "")

gtdbtk.all.merged <- merge(amr.mags.sum.filtered, gtdbtk.all.filter, by.x = "SAMPLE", by.y = "user_genome")
gtdbtk.all.merged <- gtdbtk.all.merged[, c(1, 3, 4, 5, 6, 7, 8, 9, 2)]

write.csv(gtdbtk.all.merged, file = "output/08_visualization/tab_mags_sum.csv", quote = FALSE, row.names = FALSE)
write.csv(amr.mags.overview, file = "output/08_visualization/tab_mags_mdr.csv", quote = FALSE)