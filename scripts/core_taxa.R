# --------------------------------------------------------------------------------------------------------
# Title: core_taxa.R
# Author: Silver A. Wolf
# Last Modified: Wed, 05.04.2023
# Version: 0.0.5
# --------------------------------------------------------------------------------------------------------

# Libraries
library("eulerr")
library("microbiome")
library("openxlsx")
library("taxonomizr")

# Import data
data.biom <- import_biom("../output/02_taxonomic_profiling/kraken_biom/bracken_update.biom", parseFunction = parse_taxonomy_default)
data.alpha <- microbiome::alpha(data.biom)
meta.raw <- read.xlsx("../metadata/23_03_Horses_Overview.xlsx", sheet = 1)
meta.raw <- meta.raw[,-c(1)]
colnames(meta.raw)[1] <- "SampleID"
meta.sorted = meta.raw[match(rownames(data.alpha), meta.raw$SampleID),]
rownames(meta.sorted) <- meta.sorted$SampleID
data.biom@sam_data <- sample_data(meta.sorted)

colours.days = c("t0" = "#00BA38",
                 "t1" = "#F8766D",
                 "t2" = "#619CFF"
                 )

colours.groups = c("SSG" = "#00ff7f",
                   "5DG" = "#ffa500",
                   "REF" = "#adadad"
                   )

# Taxa Tables

# Core Taxa Table
data.rel <- microbiome::transform(data.biom, "compositional")
micro.core <- core(data.rel, detection = 0, prevalence = 0.95)
core.taxa.ids <- taxa(micro.core)
core.taxa.names <- as.data.frame(getTaxonomy(core.taxa.ids, sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")))
write.csv(core.taxa.names, file = "../output/08_visualization/tab_taxa_core.csv", quote = FALSE)

# Accessory Taxa Table
data.taxa <- as.data.frame(tax_table(data.biom))
micro.accessory <- data.taxa[!(rownames(data.taxa) %in% core.taxa.ids), ]
accesssory.taxa.ids <- rownames(micro.accessory)
accesssory.taxa.names <- as.data.frame(getTaxonomy(accesssory.taxa.ids, sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")))
write.csv(accesssory.taxa.names, file = "../output/08_visualization/tab_taxa_accesssory.csv", quote = FALSE)

# 5DG Core Taxa Table
samples.5dg <- meta.sorted[meta.sorted$AB_Group == "5DG", ]$SampleID
data.5dg <- prune_samples(samples.5dg, data.biom)
data.5dg.rel <- microbiome::transform(data.5dg, "compositional")
micro.core.5dg <- core(data.5dg.rel, detection = 0, prevalence = 0.95)
core.taxa.ids.5dg <- taxa(micro.core.5dg)
core.taxa.names.5dg <- as.data.frame(getTaxonomy(core.taxa.ids.5dg, sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")))
write.csv(core.taxa.names.5dg, file = "../output/08_visualization/tab_taxa_core_5dg.csv", quote = FALSE)

# SSG Core Taxa Table
samples.ssg <- meta.sorted[meta.sorted$AB_Group == "SSG", ]$SampleID
data.ssg <- prune_samples(samples.ssg, data.biom)
data.ssg.rel <- microbiome::transform(data.ssg, "compositional")
micro.core.ssg <- core(data.ssg.rel, detection = 0, prevalence = 0.95)
core.taxa.ids.ssg <- taxa(micro.core.ssg)
core.taxa.names.ssg <- as.data.frame(getTaxonomy(core.taxa.ids.ssg, sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")))
write.csv(core.taxa.names.ssg, file = "../output/08_visualization/tab_taxa_core_ssg.csv", quote = FALSE)

# Venn Diagram Visualizations

# Venn Diagram (Groups)
pseq.rel <- microbiome::transform(data.biom, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$AB_Group))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, AB_Group == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_groups_all.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "Core Taxa (Groups)")
dev.off()

# Venn Diagram (SSG)
samples.ssg <- meta.sorted[meta.sorted$AB_Group == "SSG", ]$SampleID
data.ssg <- prune_samples(samples.ssg, data.biom)

pseq.rel <- microbiome::transform(data.ssg, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$Timepoint))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, Timepoint == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_groups_ssg.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.days, main = "Core Taxa (SSG)")
dev.off()

# Venn Diagram (5DG)
samples.5dg <- meta.sorted[meta.sorted$AB_Group == "5DG", ]$SampleID
data.5dg <- prune_samples(samples.5dg, data.biom)

pseq.rel <- microbiome::transform(data.5dg, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$Timepoint))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, Timepoint == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_groups_fdg.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.days, main = "Core Taxa (5DG)")
dev.off()

# Venn Diagram (t0)
samples.t0 <- meta.sorted[meta.sorted$Timepoint == "t0" & meta.sorted$AB_Group %in% c("SSG", "5DG"), ]$SampleID
data.t0 <- prune_samples(samples.t0, data.biom)

pseq.rel <- microbiome::transform(data.t0, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$AB_Group))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, AB_Group == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_time_t0.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups[1:2], main = "Core Taxa (t0)")
dev.off()

core_t0_unique_ssg <- as.data.frame(getTaxonomy(list_core$SSG[!(list_core$SSG %in% list_core$'5DG')], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t0_unique_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[!(list_core$'5DG' %in% list_core$SSG)], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t0_shared_ssg_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[list_core$'5DG' %in% list_core$SSG], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))

write.csv(core_t0_unique_ssg, "../output/08_visualization/tab_taxa_core_ssg_t0.csv", quote = FALSE)
write.csv(core_t0_unique_5dg, "../output/08_visualization/tab_taxa_core_5dg_t0.csv", quote = FALSE)
write.csv(core_t0_shared_ssg_5dg, "../output/08_visualization/tab_taxa_core_shared_t0.csv", quote = FALSE)

# Venn Diagram (t1)
samples.t1 <- meta.sorted[meta.sorted$Timepoint == "t1" & meta.sorted$AB_Group %in% c("SSG", "5DG"), ]$SampleID
data.t1 <- prune_samples(samples.t1, data.biom)

pseq.rel <- microbiome::transform(data.t1, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$AB_Group))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, AB_Group == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_time_t1.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups[1:2], main = "Core Taxa (t1)")
dev.off()

core_t1_unique_ssg <- as.data.frame(getTaxonomy(list_core$SSG[!(list_core$SSG %in% list_core$'5DG')], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t1_unique_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[!(list_core$'5DG' %in% list_core$SSG)], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t1_shared_ssg_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[list_core$'5DG' %in% list_core$SSG], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))

write.csv(core_t1_unique_ssg, "../output/08_visualization/tab_taxa_core_ssg_t1.csv", quote = FALSE)
write.csv(core_t1_unique_5dg, "../output/08_visualization/tab_taxa_core_5dg_t1.csv", quote = FALSE)
write.csv(core_t1_shared_ssg_5dg, "../output/08_visualization/tab_taxa_core_shared_t1.csv", quote = FALSE)

# Venn Diagram (t2)
samples.t2 <- meta.sorted[meta.sorted$Timepoint == "t2" & meta.sorted$AB_Group %in% c("SSG", "5DG"), ]$SampleID
data.t2 <- prune_samples(samples.t2, data.biom)

pseq.rel <- microbiome::transform(data.t2, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$AB_Group))

list_core <- c()

for (n in disease_states){
  ps.sub <- subset_samples(pseq.rel, AB_Group == n)
  core_m <- core_members(ps.sub, detection = 0, prevalence = 0.95)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

png("../output/08_visualization/venn_taxa_time_t2.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups[1:2], main = "Core Taxa (t2)")
dev.off()

core_t2_unique_ssg <- as.data.frame(getTaxonomy(list_core$SSG[!(list_core$SSG %in% list_core$'5DG')], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t2_unique_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[!(list_core$'5DG' %in% list_core$SSG)], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))
core_t2_shared_ssg_5dg <- as.data.frame(getTaxonomy(list_core$'5DG'[list_core$'5DG' %in% list_core$SSG], sqlFile = "../tmp/accessionTaxa.sql", desiredTaxa = c("family", "genus", "species")))

write.csv(core_t2_unique_ssg, "../output/08_visualization/tab_taxa_core_ssg_t2.csv", quote = FALSE)
write.csv(core_t2_unique_5dg, "../output/08_visualization/tab_taxa_core_5dg_t2.csv", quote = FALSE)
write.csv(core_t2_shared_ssg_5dg, "../output/08_visualization/tab_taxa_core_shared_t2.csv", quote = FALSE)