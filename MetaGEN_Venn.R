# --------------------------------------------------------------------------------------------------------
# Title: MetaGEN_Venn.R
# Author: Silver A. Wolf
# Last Modified: Thu, 23.06.2022
# Version: 0.0.1
# --------------------------------------------------------------------------------------------------------

# Libraries
library("eulerr")
library("microbiome")

data.biom <- import_biom("output/02_taxonomic_profiling/kraken_biom/bracken_update.biom", parseFunction = parse_taxonomy_default)
data.alpha <- microbiome::alpha(data.biom)
meta.raw <- read.csv("metadata/Horses_Overview.csv", sep = "\t", na.strings = "XXX")
meta.sorted = meta.raw[match(rownames(data.alpha), meta.raw$SampleID),]
rownames(meta.sorted) <- meta.sorted$SampleID
data.biom@sam_data <- sample_data(meta.sorted)

colours.groups = c("SSG" = "#00ff7f",
                   "5DG" = "#ffa500",
                   "SWITCHED" = "#00bfff",
                   "REF" = "#606060"
                   )

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

png("output/08_visualization/venn_taxa_groups.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "Groups")
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

png("output/08_visualization/venn_taxa_ssg.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "SSG")
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

png("output/08_visualization/venn_taxa_5dg.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "5DG")
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

png("output/08_visualization/venn_taxa_time_t0.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "t0")
dev.off()

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

png("output/08_visualization/venn_taxa_time_t1.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "t1")
dev.off()

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

png("output/08_visualization/venn_taxa_time_t2.png", width = 10, height = 10, units = "cm", res = 500)
plot(venn(list_core), fills = colours.groups, main = "t2")
dev.off()